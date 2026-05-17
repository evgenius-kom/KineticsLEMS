"""ICTAC consistency check: pairwise E(α) divergence."""
from __future__ import annotations

import numpy as np

from kinetics_lems.methods import IsoconversionalResult, consistency_check


def _iso(name: str, E_kJ: np.ndarray) -> IsoconversionalResult:
    alpha = np.linspace(0.1, 0.9, E_kJ.size)
    return IsoconversionalResult(
        method=name,
        alpha=alpha,
        Ea_J_per_mol=E_kJ * 1000.0,
        intercept=np.zeros_like(alpha),
        r_squared=np.ones_like(alpha),
    )


def test_identical_methods_have_zero_difference() -> None:
    E = np.full(9, 120.0)
    res = consistency_check({"a": _iso("a", E), "b": _iso("b", E)})
    assert len(res.pairs) == 1
    assert res.pairs[0].max_relative_difference == 0.0
    assert res.warnings == []


def test_pairs_exceeding_threshold_produce_warning() -> None:
    Ea = np.full(9, 100.0)
    Eb = np.full(9, 130.0)  # 30/130 ≈ 0.231 relative diff > 0.10 default
    res = consistency_check({"a": _iso("a", Ea), "b": _iso("b", Eb)})
    assert len(res.warnings) == 1
    assert "a vs b" in res.warnings[0]
    assert res.pairs[0].max_relative_difference > 0.10


def test_locates_worst_alpha() -> None:
    Ea = np.full(9, 100.0)
    Eb = np.full(9, 100.0)
    Eb[5] = 200.0  # spike in the middle
    res = consistency_check({"a": _iso("a", Ea), "b": _iso("b", Eb)})
    p = res.pairs[0]
    # 9 points over α∈[0.1, 0.9]; index 5 → α=0.6
    assert abs(p.alpha_of_max - 0.6) < 1e-6
    assert p.max_relative_difference > 0.4


def test_threshold_override() -> None:
    Ea = np.full(9, 100.0)
    Eb = np.full(9, 108.0)  # ~7.4% diff
    res_strict = consistency_check({"a": _iso("a", Ea), "b": _iso("b", Eb)}, threshold=0.05)
    res_lax = consistency_check({"a": _iso("a", Ea), "b": _iso("b", Eb)}, threshold=0.10)
    assert len(res_strict.warnings) == 1
    assert len(res_lax.warnings) == 0


def test_three_methods_produce_three_pairs() -> None:
    E = np.full(9, 120.0)
    res = consistency_check(
        {"friedman": _iso("friedman", E), "kas": _iso("kas", E), "vyazovkin": _iso("vyazovkin", E)}
    )
    assert len(res.pairs) == 3
    names = {(p.method_a, p.method_b) for p in res.pairs}
    assert names == {("friedman", "kas"), ("friedman", "vyazovkin"), ("kas", "vyazovkin")}


def test_nan_in_one_method_is_ignored() -> None:
    Ea = np.full(9, 120.0)
    Eb = np.full(9, 120.0)
    Eb[0] = np.nan  # tail NaN should be skipped, rest is identical
    res = consistency_check({"a": _iso("a", Ea), "b": _iso("b", Eb)})
    assert res.pairs[0].max_relative_difference == 0.0

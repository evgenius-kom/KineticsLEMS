"""Multi-step E(α) segmentation."""
from __future__ import annotations

import numpy as np
import pytest

from kinetics_lems.conversion import build_runs
from kinetics_lems.methods.common import IsoconversionalResult
from kinetics_lems.methods.multistep import detect_steps
from kinetics_lems.methods.vyazovkin import vyazovkin
from kinetics_lems.synthetic import generate_case


def _flat_iso(Ea_kJ: float, n: int = 19) -> IsoconversionalResult:
    a = np.linspace(0.05, 0.95, n)
    Ea = np.full(n, Ea_kJ * 1000.0)
    return IsoconversionalResult(
        method="synthetic",
        alpha=a,
        Ea_J_per_mol=Ea,
        intercept=np.full(n, np.nan),
        r_squared=np.full(n, np.nan),
    )


def _step_iso(Ea_lo_kJ: float, Ea_hi_kJ: float, split: float = 0.5) -> IsoconversionalResult:
    a = np.linspace(0.05, 0.95, 19)
    Ea = np.where(a < split, Ea_lo_kJ, Ea_hi_kJ) * 1000.0
    return IsoconversionalResult("synthetic", a, Ea, np.full(19, np.nan), np.full(19, np.nan))


def test_flat_E_yields_single_step():
    iso = _flat_iso(120.0)
    res = detect_steps(iso, relative_jump_threshold=0.05)
    assert res.n_steps == 1
    assert res.flatness_score == pytest.approx(0.0, abs=1e-12)
    assert res.steps[0].Ea_kJ_per_mol_median == pytest.approx(120.0)


def test_step_function_yields_two_steps():
    """A square jump from 100 → 180 kJ/mol at α=0.5 gives 2 steps."""
    iso = _step_iso(100.0, 180.0, split=0.5)
    res = detect_steps(iso, relative_jump_threshold=0.10, min_segment_size=2)
    assert res.n_steps == 2, f"Got {res.n_steps} steps: {res.steps}"
    medians = sorted(s.Ea_kJ_per_mol_median for s in res.steps)
    assert medians[0] == pytest.approx(100.0, abs=0.5)
    assert medians[1] == pytest.approx(180.0, abs=0.5)


def test_synthetic_single_step_F1_yields_one_step():
    """End-to-end: clean synthetic F1 should still be detected as single-step."""
    case = generate_case(
        rates_K_per_min=[5.0, 10.0, 20.0],
        Ea_J_per_mol=120_000.0,
        A_per_sec=1.0e10,
        model="F1",
        n_points=2000,
        seed=0,
    )
    runs = build_runs(case)
    alphas = np.linspace(0.1, 0.9, 17)
    iso = vyazovkin(runs, alphas)
    res = detect_steps(iso, relative_jump_threshold=0.05, min_segment_size=3)
    assert res.n_steps == 1, f"Synthetic F1 split into {res.n_steps} segments"
    assert res.flatness_score < 0.05


def test_min_segment_merges_tiny_segments():
    """Tiny noisy segments should be merged into a larger neighbour."""
    a = np.linspace(0.05, 0.95, 19)
    # Mostly 100, but a single outlier point at α≈0.5.
    Ea = np.full_like(a, 100.0)
    Ea[9] = 180.0
    iso = IsoconversionalResult(
        "synthetic", a, Ea * 1000.0,
        np.full_like(a, np.nan), np.full_like(a, np.nan),
    )
    res = detect_steps(iso, relative_jump_threshold=0.10, min_segment_size=3)
    # Outlier alone forms a 1-point segment → must be merged.
    assert res.n_steps <= 2

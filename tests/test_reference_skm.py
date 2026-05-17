"""Multi-rate stress test on SKM preimpregnated polymer.

SKM is the widest-β fixture in the corpus: six heating rates
(1, 2.5, 5, 10, 20, 40 K/min — a 40× span). The original sheet has no
reference E_a column, so this fixture is a *self-consistency* test:
different isoconversional methods must agree with each other inside a
stable α window.

Tests:

1. Vyazovkin and KAS produce smooth E(α) curves where pairwise
   disagreement is small across α ∈ [0.1, 0.9].
2. Dropping the slowest or fastest rate does not change the mean E(α)
   in the stable α window by more than a few kJ/mol — sanity check
   that the wide β span did not destabilize the minimization.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

from kinetics_lems.constants import SEC_PER_MIN
from kinetics_lems.conversion import ConversionRun
from kinetics_lems.methods import kas, vyazovkin

FIXTURES = Path(__file__).resolve().parents[1] / "data" / "reference_workbooks" / "skm_prepreg"
RATES = [1.0, 2.5, 5.0, 10.0, 20.0, 40.0]
ALPHA_GRID = np.linspace(0.10, 0.90, 17)


def _load_run(beta_K_per_min: float) -> ConversionRun:
    data = np.loadtxt(FIXTURES / f"rate_{beta_K_per_min}.tsv", skiprows=1)
    alpha = data[:, 0]
    T = data[:, 1]
    dadT = np.gradient(alpha, T)
    return ConversionRun(
        rate_K_per_min=float(beta_K_per_min),
        temperature=T,
        alpha=alpha,
        dalpha_dt=dadT * (beta_K_per_min / SEC_PER_MIN),
    )


def test_skm_loads_all_6_rates():
    """The fixture exists and parses to 100 (α, T) points per rate."""
    for beta in RATES:
        run = _load_run(beta)
        assert run.alpha.size == 100, f"rate {beta}: got {run.alpha.size} points"
        assert run.alpha.min() == 0.005
        assert run.alpha.max() == 0.995


def test_vyazovkin_kas_agree_on_skm():
    """Vyazovkin and KAS should agree to within ~5 kJ/mol on this clean fixture."""
    runs = [_load_run(b) for b in RATES]
    res_v = vyazovkin(runs, ALPHA_GRID)
    res_k = kas(runs, ALPHA_GRID)
    diff = np.abs(res_v.Ea_kJ_per_mol - res_k.Ea_kJ_per_mol)
    mean = float(np.mean(diff))
    mx = float(np.max(diff))
    assert mean < 3.0, f"mean |Vyazovkin − KAS| = {mean:.2f} kJ/mol"
    assert mx < 5.0, f"max |Vyazovkin − KAS| = {mx:.2f} kJ/mol"


def test_skm_E_is_finite_everywhere():
    """No NaN / Inf in the recovered E(α) — basic numerical robustness."""
    runs = [_load_run(b) for b in RATES]
    res_v = vyazovkin(runs, ALPHA_GRID)
    res_k = kas(runs, ALPHA_GRID)
    assert np.all(np.isfinite(res_v.Ea_kJ_per_mol))
    assert np.all(np.isfinite(res_k.Ea_kJ_per_mol))


def test_skm_robust_to_dropping_one_rate():
    """Leave-one-rate-out E(α) should stay close to the full-data result
    in the stable α window — confirms no single rate dominates the fit."""
    runs_all = [_load_run(b) for b in RATES]
    res_full = vyazovkin(runs_all, ALPHA_GRID)
    for skip_idx, skip_beta in enumerate(RATES):
        subset = [r for i, r in enumerate(runs_all) if i != skip_idx]
        res_sub = vyazovkin(subset, ALPHA_GRID)
        diff = np.abs(res_full.Ea_kJ_per_mol - res_sub.Ea_kJ_per_mol)
        # Dropping the extreme rates (1 or 40) can shift the curve more,
        # but never beyond a few kJ/mol on this clean fixture.
        assert float(np.mean(diff)) < 5.0, (
            f"dropping {skip_beta} K/min shifts mean E(α) by "
            f"{float(np.mean(diff)):.2f} kJ/mol — should be smooth"
        )

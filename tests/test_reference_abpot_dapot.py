"""Regression tests against ABPOT / DAPOT reference α(T) fixtures.

The fixtures under ``data/reference_workbooks/{abpot,dapot}/`` contain
(α, T) tables extracted from an independent Excel-based Vyazovkin
workbook plus a per-α reference E_a column. We feed those tables into
our integral methods and assert agreement within material-specific
tolerances.

Why integral methods only: the fixture is pre-processed α(T), not raw DSC.
Differential methods (Friedman, Vyazovkin-AIC) need fine-grained dα/dt
that this format throws away — those are exercised on synthetic raw waves
in ``test_synthetic_recovery.py``.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from kinetics_lems.constants import SEC_PER_MIN
from kinetics_lems.conversion import ConversionRun
from kinetics_lems.methods import IsoconversionalResult, kas, ofw, vyazovkin

FIXTURES = Path(__file__).resolve().parents[1] / "data" / "reference_workbooks"


def _load_run(folder: Path, beta_K_per_min: float) -> ConversionRun:
    data = np.loadtxt(folder / f"rate_{beta_K_per_min}.tsv", skiprows=1)
    alpha = data[:, 0]
    T = data[:, 1]
    beta_per_sec = beta_K_per_min / SEC_PER_MIN
    dadT = np.gradient(alpha, T)
    return ConversionRun(
        rate_K_per_min=float(beta_K_per_min),
        temperature=T,
        alpha=alpha,
        dalpha_dt=dadT * beta_per_sec,
    )


def _load_reference(folder: Path) -> tuple[np.ndarray, np.ndarray]:
    data = np.loadtxt(folder / "reference_Ea.tsv", skiprows=1)
    return data[:, 0], data[:, 1]


def _stats_in_window(
    res: IsoconversionalResult,
    ref_alpha: np.ndarray,
    ref_Ea: np.ndarray,
    alpha_lo: float,
    alpha_hi: float,
) -> tuple[float, float]:
    """Return (mean_abs_diff, max_abs_diff) across α ∈ [alpha_lo, alpha_hi]."""
    mask_eval = (res.alpha >= alpha_lo) & (res.alpha <= alpha_hi)
    mask_ref = (ref_alpha >= alpha_lo) & (ref_alpha <= alpha_hi)
    # Match on common α points (fixtures and our analysis use the same grid).
    common = np.intersect1d(res.alpha[mask_eval], ref_alpha[mask_ref])
    ours = np.array([res.Ea_kJ_per_mol[res.alpha == a][0] for a in common])
    refs = np.array([ref_Ea[ref_alpha == a][0] for a in common])
    diff = np.abs(ours - refs)
    return float(np.mean(diff)), float(np.max(diff))


# Per-material spec: (folder name, list of heating rates, evaluation window,
# tolerances per method). Tolerances chosen ~30% above currently-observed
# errors so the test is sensitive to genuine drift but not to floating-point
# jitter or a one-iteration change in scipy's bounded minimizer.
_CASES = [
    pytest.param(
        "abpot",
        [0.2, 1.0, 2.5, 5.0],
        (0.10, 0.80),
        {
            "kas": (0.5, 2.0),         # mean, max kJ/mol
            "vyazovkin": (0.5, 2.0),
            "ofw": (1.5, 3.0),
        },
        id="abpot",
    ),
    pytest.param(
        "dapot",
        [1.0, 2.5, 5.0, 10.0],
        (0.10, 0.80),
        {
            "kas": (2.0, 5.0),
            "vyazovkin": (2.0, 5.0),
            "ofw": (3.0, 6.0),
        },
        id="dapot",
    ),
]


@pytest.mark.parametrize("material,rates,window,tolerances", _CASES)
def test_integral_methods_match_reference(material, rates, window, tolerances):
    folder = FIXTURES / material
    runs = [_load_run(folder, b) for b in rates]
    ref_alpha, ref_Ea = _load_reference(folder)
    # Evaluate on the reference's α-grid.
    alphas = ref_alpha

    methods = {"kas": kas, "vyazovkin": vyazovkin, "ofw": ofw}
    failures: list[str] = []
    for name, fn in methods.items():
        res = fn(runs, alphas)
        mean_d, max_d = _stats_in_window(res, ref_alpha, ref_Ea, *window)
        mean_tol, max_tol = tolerances[name]
        if mean_d > mean_tol or max_d > max_tol:
            failures.append(
                f"{name}: mean Δ = {mean_d:.2f} (tol {mean_tol}), "
                f"max Δ = {max_d:.2f} (tol {max_tol})"
            )
    assert not failures, (
        f"{material} drifted from reference in α∈{window}:\n  " + "\n  ".join(failures)
    )

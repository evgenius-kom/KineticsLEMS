"""Validation test on the Epoxy PV15-0,5 reference fixture.

Same structure as ``test_reference_abpot_dapot.py``: feed the (α, T)
tables of all 3 heating rates into our integral methods and assert
agreement with the per-α reference E_a column from an independent
Excel-based Vyazovkin computation.

Epoxy PV15-0,5 is qualitatively different from ABPOT/DAPOT — its
reference E_a is a U-shape (≈52 → 41 → 46 kJ/mol across α), so this
fixture additionally exercises the methods on a non-flat E(α) target.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

from kinetics_lems.constants import SEC_PER_MIN
from kinetics_lems.conversion import ConversionRun
from kinetics_lems.methods import IsoconversionalResult, kas, ofw, vyazovkin

FIXTURES = Path(__file__).resolve().parents[1] / "data" / "reference_workbooks" / "epoxy_pv15"
RATES = [2.5, 5.0, 10.0]
ALPHA_WINDOW = (0.10, 0.80)


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


def _load_reference() -> tuple[np.ndarray, np.ndarray]:
    data = np.loadtxt(FIXTURES / "reference_Ea.tsv", skiprows=1)
    return data[:, 0], data[:, 1]


def _stats(
    res: IsoconversionalResult,
    ref_alpha: np.ndarray,
    ref_Ea: np.ndarray,
) -> tuple[float, float]:
    lo, hi = ALPHA_WINDOW
    common = np.intersect1d(
        res.alpha[(res.alpha >= lo) & (res.alpha <= hi)],
        ref_alpha[(ref_alpha >= lo) & (ref_alpha <= hi)],
    )
    ours = np.array([res.Ea_kJ_per_mol[res.alpha == a][0] for a in common])
    refs = np.array([ref_Ea[ref_alpha == a][0] for a in common])
    diff = np.abs(ours - refs)
    return float(np.mean(diff)), float(np.max(diff))


def test_epoxy_integral_methods_match_reference():
    runs = [_load_run(b) for b in RATES]
    ref_alpha, ref_Ea = _load_reference()
    alphas = ref_alpha

    # Tolerances chosen ~30% above currently-observed errors on this fixture.
    # The non-flat E(α) profile is harder than ABPOT/DAPOT's flat target, so
    # tolerances are wider — especially for OFW which uses Doyle's approx and
    # is biased when E_a varies across α.
    tolerances = {
        "kas":       (1.5, 4.0),
        "vyazovkin": (1.5, 4.0),
        "ofw":       (5.0, 7.0),
    }
    methods = {"kas": kas, "vyazovkin": vyazovkin, "ofw": ofw}
    failures: list[str] = []
    for name, fn in methods.items():
        res = fn(runs, alphas)
        mean_d, max_d = _stats(res, ref_alpha, ref_Ea)
        mean_tol, max_tol = tolerances[name]
        if mean_d > mean_tol or max_d > max_tol:
            failures.append(
                f"{name}: mean Δ = {mean_d:.2f} (tol {mean_tol}), "
                f"max Δ = {max_d:.2f} (tol {max_tol})"
            )
    assert not failures, (
        f"epoxy drifted from reference in α∈{ALPHA_WINDOW}:\n  " + "\n  ".join(failures)
    )


def test_epoxy_E_profile_is_non_flat():
    """Sanity check that this fixture exercises non-flat E(α) territory.

    Epoxy E(α) starts at ~52 kJ/mol, declines to ~41 around α≈0.9, then
    spikes up to ~46 at α=1. The range alone is ~10 kJ/mol, more than
    twice the flat-curve scatter seen on ABPOT/DAPOT."""
    _, ref_Ea = _load_reference()
    e_range = float(np.max(ref_Ea) - np.min(ref_Ea))
    assert e_range > 5.0, f"E_a range only {e_range:.1f} kJ/mol — fixture not as expected"

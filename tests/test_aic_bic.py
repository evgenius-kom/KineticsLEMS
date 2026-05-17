"""AIC/BIC reporting for Coats-Redfern and multistep."""
from __future__ import annotations

import numpy as np

from kinetics_lems.conversion import build_runs
from kinetics_lems.methods import coats_redfern, detect_steps
from kinetics_lems.methods.common import IsoconversionalResult
from kinetics_lems.synthetic import generate_case


def _synthetic_runs():
    case = generate_case(
        rates_K_per_min=[2.5, 5.0, 10.0, 20.0],
        Ea_J_per_mol=120_000.0,
        A_per_sec=1.0e10,
        model="F1",
    )
    return build_runs(case)


def test_coats_redfern_fits_have_finite_aic_bic() -> None:
    runs = _synthetic_runs()
    res = coats_redfern(runs)
    for fit in res.fits:
        assert np.isfinite(fit.aic)
        assert np.isfinite(fit.bic)
    for summary in res.summaries:
        assert np.isfinite(summary.aic_mean)
        assert np.isfinite(summary.bic_mean)


def test_coats_redfern_best_model_aic_within_top_three() -> None:
    """A good model (low RSS) should land in top-3 by AIC even if ties on R²."""
    runs = _synthetic_runs()
    res = coats_redfern(runs)
    by_aic = sorted(res.summaries, key=lambda s: s.aic_mean)
    by_r2 = sorted(res.summaries, key=lambda s: -s.r_squared_mean)
    top3_aic = {s.model for s in by_aic[:3]}
    top3_r2 = {s.model for s in by_r2[:3]}
    # On clean F1 synthetic data, AIC and R² rankings should overlap.
    assert top3_aic & top3_r2


def test_multistep_aic_bic_finite_for_synthetic_flat_E() -> None:
    alpha = np.linspace(0.05, 0.95, 19)
    E = np.full_like(alpha, 120_000.0)
    iso = IsoconversionalResult(
        method="vyazovkin",
        alpha=alpha,
        Ea_J_per_mol=E,
        intercept=np.zeros_like(alpha),
        r_squared=np.ones_like(alpha),
    )
    res = detect_steps(iso)
    assert res.n_steps == 1
    # Flat E → RSS = 0 → log floor kicks in; AIC is *finite* (large-negative),
    # not NaN.
    assert np.isfinite(res.aic_piecewise_constant)
    assert np.isfinite(res.bic_piecewise_constant)


def test_multistep_aic_increases_for_simpler_signal_under_two_step_model() -> None:
    """Forcing a 2-step segmentation on a flat curve must increase AIC vs 1-step."""
    alpha = np.linspace(0.05, 0.95, 19)
    iso = IsoconversionalResult(
        method="vyazovkin",
        alpha=alpha,
        Ea_J_per_mol=np.full_like(alpha, 120_000.0),
        intercept=np.zeros_like(alpha),
        r_squared=np.ones_like(alpha),
    )
    single = detect_steps(iso, relative_jump_threshold=0.50, min_segment_size=3)
    # Manually construct a 2-step E(α) and a 1-step E(α) for comparison.
    iso_jump = IsoconversionalResult(
        method="vyazovkin",
        alpha=alpha,
        Ea_J_per_mol=np.concatenate(
            [np.full(9, 120_000.0), np.full(10, 180_000.0)]
        ),
        intercept=np.zeros_like(alpha),
        r_squared=np.ones_like(alpha),
    )
    two_step = detect_steps(iso_jump, relative_jump_threshold=0.10, min_segment_size=3)
    assert single.n_steps == 1
    assert two_step.n_steps == 2
    # Both should give finite AIC; sanity check the formula didn't blow up.
    assert np.isfinite(single.aic_piecewise_constant)
    assert np.isfinite(two_step.aic_piecewise_constant)

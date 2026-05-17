"""Model-free isothermal prediction + arbitrary T(t) triplet prediction."""
from __future__ import annotations

import numpy as np

from kinetics_lems.constants import R_GAS, SEC_PER_MIN
from kinetics_lems.conversion import build_runs, temperature_at_conversion
from kinetics_lems.methods import (
    predict_alpha_of_t,
    predict_arbitrary_program_modelfree,
    predict_isothermal_modelfree,
    predict_under_program,
    vyazovkin,
)
from kinetics_lems.synthetic import generate_case


def _setup():
    case = generate_case(
        rates_K_per_min=[2.5, 5.0, 10.0, 20.0],
        Ea_J_per_mol=120_000.0,
        A_per_sec=1.0e10,
        model="F1",
    )
    runs = build_runs(case)
    alphas = np.linspace(0.1, 0.9, 17)
    iso = vyazovkin(runs, alphas)
    return runs, alphas, iso


def test_modelfree_isothermal_prediction_returns_monotone_t() -> None:
    runs, alphas, iso = _setup()
    T_alpha = temperature_at_conversion(runs[0], alphas)
    pred = predict_isothermal_modelfree(
        T_K=300.0 + 273.15,
        Ea_J_per_mol=iso.Ea_J_per_mol,
        T_alpha_K=T_alpha,
        beta_K_per_min=runs[0].rate_K_per_min,
        alpha_grid=alphas,
    )
    assert pred.time_sec.size == alphas.size
    assert np.all(np.diff(pred.time_sec) >= -1e-9)


def test_modelfree_isothermal_matches_f1_triplet_within_factor() -> None:
    """On clean F1 data, the model-free path should be within an order of
    magnitude of the triplet-based prediction (they disagree only by the
    accuracy of the f(α) shape, which is 'exact' for F1)."""
    runs, alphas, iso = _setup()
    T_alpha = temperature_at_conversion(runs[0], alphas)
    T_iso = 500.0  # K — well inside the kinetic regime, predictions non-trivial

    modelfree = predict_isothermal_modelfree(
        T_K=T_iso,
        Ea_J_per_mol=iso.Ea_J_per_mol,
        T_alpha_K=T_alpha,
        beta_K_per_min=runs[0].rate_K_per_min,
        alpha_grid=alphas,
    )
    triplet = predict_alpha_of_t(
        T_K=T_iso,
        Ea_J_per_mol=iso.Ea_J_per_mol,
        alpha_grid=alphas,
        A_per_sec=1.0e10,
        model="F1",
    )
    # Compare times at α=0.5.
    t_mf = modelfree.time_to_alpha(0.5)
    t_tr = triplet.time_to_alpha(0.5)
    assert 0.1 < t_mf / t_tr < 10.0


def test_modelfree_arbitrary_program_isothermal_matches_isothermal_helper() -> None:
    runs, alphas, iso = _setup()
    T_alpha = temperature_at_conversion(runs[0], alphas)
    # Build a constant-T program over 10 hours.
    time_s = np.linspace(0.0, 10 * 3600.0, 4001)
    T_program = np.full_like(time_s, 500.0)

    iso_pred = predict_isothermal_modelfree(
        T_K=500.0,
        Ea_J_per_mol=iso.Ea_J_per_mol,
        T_alpha_K=T_alpha,
        beta_K_per_min=runs[0].rate_K_per_min,
        alpha_grid=alphas,
    )
    arb_pred = predict_arbitrary_program_modelfree(
        time_s=time_s,
        T_K_program=T_program,
        Ea_J_per_mol=iso.Ea_J_per_mol,
        T_alpha_K=T_alpha,
        beta_K_per_min=runs[0].rate_K_per_min,
        alpha_grid=alphas,
    )
    # The two should give nearly identical t(α) for the α values that
    # were reached within the 10 h window.
    finite = np.isfinite(arb_pred.time_sec)
    if finite.any():
        np.testing.assert_allclose(
            arb_pred.time_sec[finite], iso_pred.time_sec[finite], rtol=0.05
        )


def test_predict_under_program_returns_alpha_in_unit_interval() -> None:
    runs, alphas, iso = _setup()
    # Ramp 25→200 °C over 1 hour.
    time_s = np.linspace(0.0, 3600.0, 601)
    T_K = 298.15 + (475.0 - 298.15) * (time_s / time_s[-1])
    pred = predict_under_program(
        time_s=time_s,
        T_K=T_K,
        Ea_J_per_mol=iso.Ea_J_per_mol,
        alpha_grid=alphas,
        A_per_sec=1.0e10,
        model="F1",
    )
    assert pred.alpha.shape == time_s.shape
    assert np.all((pred.alpha >= 0.0) & (pred.alpha <= 1.0))
    assert np.all(np.diff(pred.alpha) >= -1e-9)


def test_predict_under_program_linear_heating_matches_synthetic_alpha() -> None:
    """When T(t) is exactly the same linear ramp as the experiment, the
    predicted α(T) curve should match the experimental α(T) within a few %."""
    runs, alphas, iso = _setup()
    run = runs[1]  # 5 K/min
    beta_K_per_sec = run.rate_K_per_min / SEC_PER_MIN
    # Use the same temperature grid; build a time grid from it.
    T_K = run.temperature
    time_s = (T_K - T_K[0]) / beta_K_per_sec
    pred = predict_under_program(
        time_s=time_s,
        T_K=T_K,
        Ea_J_per_mol=iso.Ea_J_per_mol,
        alpha_grid=alphas,
        A_per_sec=1.0e10,
        model="F1",
    )
    # Compare experimental α(T) vs predicted α(T) at α = 0.5.
    T_exp_half = float(np.interp(0.5, run.alpha, T_K))
    T_pred_half = float(np.interp(0.5, pred.alpha, T_K))
    rel_err = abs(T_exp_half - T_pred_half) / T_exp_half
    # E(α) reconstructed from Vyazovkin is not exact; ±2% on T is plenty.
    assert rel_err < 0.02


def _F1_analytical_time(alpha: float, T_K: float, A: float, E_J: float) -> float:
    """Closed-form: -ln(1-α) = A · exp(-E/RT) · t  →  t = -ln(1-α) / k."""
    k = A * np.exp(-E_J / (R_GAS * T_K))
    return -np.log(max(1.0 - alpha, 1e-12)) / k

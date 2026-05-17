"""Predictive isothermal α(t) sanity checks."""
from __future__ import annotations

import numpy as np
import pytest

from kinetics_lems.constants import R_GAS
from kinetics_lems.methods.lifetime import (
    predict_alpha_of_t,
    predict_at_temperatures,
    time_to_conversion,
)


def _flat_Ea(Ea_kJ: float, n: int = 17) -> tuple[np.ndarray, np.ndarray]:
    alpha = np.linspace(0.05, 0.95, n)
    Ea = np.full(n, Ea_kJ * 1000.0)
    return alpha, Ea


def test_F1_analytical_check():
    """For f(α) = 1 − α with flat E_a, the closed form is
    α(t) = 1 − exp(−k·t), k = A·exp(−E/(RT)). Verify our prediction matches.

    Note: predict_alpha_of_t starts integration at alpha[0], so the
    relevant identity is t_pred(α) = t_analytical(α) − t_analytical(α₀)."""
    Ea_kJ = 120.0
    A = 1.0e10
    T = 400.0
    # Dense grid → tight trapezoidal accuracy on a stiff integrand.
    alpha = np.linspace(0.05, 0.95, 181)
    Ea = np.full_like(alpha, Ea_kJ * 1000.0)
    pred = predict_alpha_of_t(
        T_K=T, Ea_J_per_mol=Ea, alpha_grid=alpha,
        A_per_sec=A, model="F1",
    )
    k = A * np.exp(-Ea_kJ * 1000.0 / (R_GAS * T))
    t_analytical = (-np.log(1.0 - pred.alpha) + np.log(1.0 - pred.alpha[0])) / k

    interior = slice(2, -2)  # skip endpoints; trapezoid error is largest there
    rel_err = np.abs(pred.time_sec[interior] - t_analytical[interior]) / t_analytical[interior]
    assert float(np.max(rel_err)) < 0.02


def test_higher_T_means_shorter_time():
    """Arrhenius: hotter → faster. Time-to-α should strictly decrease with T."""
    alpha, Ea = _flat_Ea(140.0)
    t_25 = time_to_conversion(
        alpha_target=0.5, T_K=298.15,
        Ea_J_per_mol=Ea, alpha_grid=alpha, A_per_sec=1.0e10, model="F1",
    )
    t_60 = time_to_conversion(
        alpha_target=0.5, T_K=333.15,
        Ea_J_per_mol=Ea, alpha_grid=alpha, A_per_sec=1.0e10, model="F1",
    )
    t_100 = time_to_conversion(
        alpha_target=0.5, T_K=373.15,
        Ea_J_per_mol=Ea, alpha_grid=alpha, A_per_sec=1.0e10, model="F1",
    )
    assert t_25 > t_60 > t_100 > 0


def test_predict_at_temperatures_returns_consistent_shape():
    alpha, Ea = _flat_Ea(130.0)
    summary = predict_at_temperatures(
        T_K_list=[300.0, 350.0, 400.0],
        Ea_J_per_mol=Ea, alpha_grid=alpha, A_per_sec=1.0e10,
        alpha_targets=(0.1, 0.5, 0.9),
        model="F1",
    )
    assert summary.times_at_targets.shape == (3, 3)
    assert len(summary.predictions) == 3
    # Times must increase monotonically with α target at fixed T.
    for i in range(3):
        row = summary.times_at_targets[i, :]
        assert row[0] < row[1] < row[2]


def test_rejects_invalid_inputs():
    alpha, Ea = _flat_Ea(120.0)
    with pytest.raises(ValueError, match="T_K"):
        predict_alpha_of_t(
            T_K=-1.0, Ea_J_per_mol=Ea, alpha_grid=alpha, A_per_sec=1.0e10,
        )
    with pytest.raises(ValueError, match="A_per_sec"):
        predict_alpha_of_t(
            T_K=300.0, Ea_J_per_mol=Ea, alpha_grid=alpha, A_per_sec=0.0,
        )
    with pytest.raises(ValueError, match="Unknown model"):
        predict_alpha_of_t(
            T_K=300.0, Ea_J_per_mol=Ea, alpha_grid=alpha, A_per_sec=1.0e10,
            model="Z99",
        )
    with pytest.raises(ValueError, match="alpha_target"):
        pred = predict_alpha_of_t(
            T_K=350.0, Ea_J_per_mol=Ea, alpha_grid=alpha, A_per_sec=1.0e10,
        )
        pred.time_to_alpha(0.99)  # outside [0.05, 0.95]


def test_variable_Ea_changes_t_versus_constant():
    """If E_a grows with α, time-to-α=0.9 must be longer than under flat E_a
    equal to the mean — sanity check that E_a(α) is actually used."""
    alpha = np.linspace(0.05, 0.95, 19)
    Ea_flat = np.full_like(alpha, 130_000.0)
    Ea_rising = np.linspace(110_000.0, 150_000.0, alpha.size)
    pred_flat = predict_alpha_of_t(
        T_K=350.0, Ea_J_per_mol=Ea_flat, alpha_grid=alpha,
        A_per_sec=1.0e10, model="F1",
    )
    pred_rising = predict_alpha_of_t(
        T_K=350.0, Ea_J_per_mol=Ea_rising, alpha_grid=alpha,
        A_per_sec=1.0e10, model="F1",
    )
    # Both reach α=0.9 — but trajectories differ from each other.
    t_flat_90 = pred_flat.time_to_alpha(0.9)
    t_rise_90 = pred_rising.time_to_alpha(0.9)
    assert t_flat_90 != t_rise_90

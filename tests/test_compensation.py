"""Compensation effect ln A ≈ a + b·E."""
from __future__ import annotations

import numpy as np
import pytest

from kinetics_lems.conversion import build_runs
from kinetics_lems.methods import (
    coats_redfern,
    compensation_from_coats_redfern,
    compensation_from_isoconversional,
)
from kinetics_lems.synthetic import generate_case


def test_compensation_from_synthetic_clean_data_has_high_r2() -> None:
    case = generate_case(
        rates_K_per_min=[2.5, 5.0, 10.0, 20.0],
        Ea_J_per_mol=120_000.0,
        A_per_sec=1.0e10,
        model="F1",
    )
    runs = build_runs(case)
    cr = coats_redfern(runs)
    fit = compensation_from_coats_redfern(cr)
    assert fit.n_points >= 4
    # The compensation line should be nearly perfect across CR models on
    # clean synthetic data (the (E, ln A) pairs collapse onto a line).
    assert fit.r_squared > 0.9


def test_compensation_from_isoconversional_basic() -> None:
    E = np.array([100.0, 120.0, 140.0, 160.0])
    # Perfect compensation: ln A = 1.0 + 0.1 · E
    A = np.exp(1.0 + 0.1 * E)
    fit = compensation_from_isoconversional(E, A)
    assert abs(fit.slope - 0.1) < 1e-9
    assert abs(fit.intercept - 1.0) < 1e-9
    assert fit.r_squared > 0.99999


def test_compensation_filters_non_finite_inputs() -> None:
    E = np.array([100.0, 120.0, np.nan, 160.0])
    A = np.array([1e8, 1e9, 1e10, 1e11])
    fit = compensation_from_isoconversional(E, A)
    assert fit.n_points == 3


def test_compensation_raises_without_enough_points() -> None:
    with pytest.raises(ValueError):
        compensation_from_isoconversional(np.array([100.0]), np.array([1e8]))

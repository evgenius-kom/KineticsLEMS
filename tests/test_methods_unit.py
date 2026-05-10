"""Unit-level tests for individual method primitives."""
from __future__ import annotations

import numpy as np
import pytest

from kinetics_lems.conversion import build_runs
from kinetics_lems.methods.common import linear_regression
from kinetics_lems.methods.friedman import friedman
from kinetics_lems.synthetic import generate_case


def test_linear_regression_recovers_known_line():
    x = np.linspace(0, 10, 50)
    slope_true, intercept_true = -2.5, 7.0
    y = slope_true * x + intercept_true
    slope, intercept, r2 = linear_regression(x, y)
    assert pytest.approx(slope, abs=1e-9) == slope_true
    assert pytest.approx(intercept, abs=1e-9) == intercept_true
    assert pytest.approx(r2, abs=1e-12) == 1.0


def test_build_runs_alpha_is_monotone_and_in_range():
    case = generate_case([5.0, 10.0], n_points=1000, seed=0)
    runs = build_runs(case)
    assert [r.rate_K_per_min for r in runs] == [5.0, 10.0]
    for r in runs:
        assert r.alpha[0] >= 0.0 and r.alpha[-1] <= 1.0
        assert np.all(np.diff(np.maximum.accumulate(r.alpha)) >= 0)


def test_friedman_with_two_runs_only():
    """Need at least 2 heating rates for any isoconversional method."""
    case = generate_case([5.0, 10.0], Ea_J_per_mol=100_000.0, n_points=1000, seed=0)
    runs = build_runs(case)
    alphas = np.linspace(0.1, 0.9, 9)
    res = friedman(runs, alphas)
    assert pytest.approx(np.nanmean(res.Ea_kJ_per_mol), abs=2.0) == 100.0

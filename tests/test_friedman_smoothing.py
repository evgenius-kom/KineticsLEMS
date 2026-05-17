"""Savitzky-Golay smoothing option in Friedman."""
from __future__ import annotations

import numpy as np
import pytest

from kinetics_lems.conversion import build_runs
from kinetics_lems.methods import friedman
from kinetics_lems.synthetic import generate_case


def _noisy_case(noise: float, seed: int = 42):
    case = generate_case(
        rates_K_per_min=[2.5, 5.0, 10.0, 20.0],
        Ea_J_per_mol=120_000.0,
        A_per_sec=1.0e10,
        model="F1",
        noise_std=noise,
        seed=seed,
    )
    return build_runs(case)


def test_friedman_default_path_unchanged_without_smoothing() -> None:
    runs = _noisy_case(0.0)
    alphas = np.linspace(0.1, 0.9, 17)
    raw = friedman(runs, alphas)
    again = friedman(runs, alphas, smooth_window=None)
    np.testing.assert_array_equal(raw.Ea_J_per_mol, again.Ea_J_per_mol)


def test_smoothing_reduces_noise_on_noisy_synthetic_data() -> None:
    runs = _noisy_case(0.05, seed=123)
    alphas = np.linspace(0.2, 0.8, 13)
    raw = friedman(runs, alphas)
    smoothed = friedman(runs, alphas, smooth_window=15, smooth_poly=3)

    # Both should be near the ground truth 120 kJ/mol; smoothed version
    # should have lower variance across α.
    truth = 120_000.0
    finite_raw = raw.Ea_J_per_mol[np.isfinite(raw.Ea_J_per_mol)]
    finite_smooth = smoothed.Ea_J_per_mol[np.isfinite(smoothed.Ea_J_per_mol)]
    assert abs(np.mean(finite_raw) - truth) < 30_000.0
    assert abs(np.mean(finite_smooth) - truth) < 30_000.0
    assert np.std(finite_smooth) <= np.std(finite_raw)


def test_invalid_window_size_raises() -> None:
    runs = _noisy_case(0.0)
    alphas = np.linspace(0.1, 0.9, 5)
    with pytest.raises(ValueError):
        friedman(runs, alphas, smooth_window=4)  # even
    with pytest.raises(ValueError):
        friedman(runs, alphas, smooth_window=3, smooth_poly=3)  # too small for poly


def test_smoothing_window_larger_than_run_is_shrunk_gracefully() -> None:
    runs = _noisy_case(0.0)
    alphas = np.linspace(0.1, 0.9, 5)
    # Window enormous → internal fallback to longest valid odd window.
    result = friedman(runs, alphas, smooth_window=10_001, smooth_poly=3)
    assert result.Ea_J_per_mol.shape == alphas.shape

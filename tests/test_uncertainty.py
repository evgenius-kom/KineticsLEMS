"""Jackknife-by-run uncertainty for E(α)."""
from __future__ import annotations

import numpy as np
import pytest

from kinetics_lems.conversion import build_runs
from kinetics_lems.methods.uncertainty import jackknife_isoconversional
from kinetics_lems.methods.vyazovkin import vyazovkin
from kinetics_lems.synthetic import generate_case


def test_returns_nan_se_with_only_two_runs():
    """Need ≥ 3 runs for jackknife. With 2, every LOO subset has 1 run."""
    case = generate_case(
        rates_K_per_min=[5.0, 10.0],
        Ea_J_per_mol=120_000.0,
        A_per_sec=1.0e10,
        model="F1",
        n_points=1000,
        seed=0,
    )
    runs = build_runs(case)
    alphas = np.linspace(0.1, 0.9, 9)
    result = jackknife_isoconversional(runs, alphas, vyazovkin)
    assert result.n_runs == 2
    assert np.all(np.isnan(result.Ea_kJ_per_mol_se))


def test_jackknife_se_small_on_noiseless_data():
    """With 4 noiseless rates, leave-one-out estimators should agree → SE ≈ 0."""
    case = generate_case(
        rates_K_per_min=[2.5, 5.0, 10.0, 20.0],
        Ea_J_per_mol=130_000.0,
        A_per_sec=1.0e10,
        model="F1",
        n_points=2000,
        seed=0,
    )
    runs = build_runs(case)
    alphas = np.linspace(0.1, 0.9, 9)
    result = jackknife_isoconversional(runs, alphas, vyazovkin)
    assert result.n_runs == 4
    finite = np.isfinite(result.Ea_kJ_per_mol_se)
    assert finite.all()
    assert float(np.mean(result.Ea_kJ_per_mol_se[finite])) < 0.5
    assert float(np.mean(result.Ea_kJ_per_mol_mean[finite])) == pytest.approx(130.0, abs=0.5)


def test_jackknife_se_grows_with_noise():
    """Adding noise should inflate the SE — sanity check the mechanism."""
    case_clean = generate_case(
        rates_K_per_min=[2.5, 5.0, 10.0, 20.0],
        Ea_J_per_mol=130_000.0,
        A_per_sec=1.0e10,
        model="F1",
        n_points=2000,
        noise_std=0.0,
        seed=0,
    )
    case_noisy = generate_case(
        rates_K_per_min=[2.5, 5.0, 10.0, 20.0],
        Ea_J_per_mol=130_000.0,
        A_per_sec=1.0e10,
        model="F1",
        n_points=2000,
        noise_std=0.05,
        seed=1,
    )
    alphas = np.linspace(0.2, 0.8, 7)
    clean = jackknife_isoconversional(build_runs(case_clean), alphas, vyazovkin)
    noisy = jackknife_isoconversional(build_runs(case_noisy), alphas, vyazovkin)
    finite_clean = np.isfinite(clean.Ea_kJ_per_mol_se)
    finite_noisy = np.isfinite(noisy.Ea_kJ_per_mol_se)
    mean_clean = float(np.mean(clean.Ea_kJ_per_mol_se[finite_clean]))
    mean_noisy = float(np.mean(noisy.Ea_kJ_per_mol_se[finite_noisy]))
    assert mean_noisy >= mean_clean, (
        f"noisy SE ({mean_noisy:.3f}) should be ≥ clean SE ({mean_clean:.3f})"
    )


def test_ci_bounds_bracket_mean():
    case = generate_case(
        rates_K_per_min=[2.5, 5.0, 10.0, 20.0],
        Ea_J_per_mol=130_000.0,
        A_per_sec=1.0e10,
        model="F1",
        n_points=2000,
        noise_std=0.03,
        seed=2,
    )
    runs = build_runs(case)
    alphas = np.linspace(0.2, 0.8, 7)
    result = jackknife_isoconversional(runs, alphas, vyazovkin)
    finite = np.isfinite(result.Ea_kJ_per_mol_se)
    lo = result.Ea_kJ_per_mol_ci95_low[finite]
    hi = result.Ea_kJ_per_mol_ci95_high[finite]
    mean = result.Ea_kJ_per_mol_mean[finite]
    assert np.all(lo <= mean)
    assert np.all(mean <= hi)

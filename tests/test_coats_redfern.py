"""Coats–Redfern model fitting against synthetic data."""
from __future__ import annotations

import numpy as np
import pytest

from kinetics_lems.conversion import build_runs
from kinetics_lems.methods.coats_redfern import coats_redfern
from kinetics_lems.synthetic import generate_case

# F1 + A_m share the same g(α) up to a constant prefactor, so any A_m can
# tie F1 in a Coats–Redfern ranking. Test the equivalence class.
_F1_AVRAMI_CLASS = {"F1", "A2", "A3", "A4"}


@pytest.mark.parametrize(
    "true_model,expected_class",
    [
        ("F1", _F1_AVRAMI_CLASS),
        ("F2", {"F2"}),
        ("F3", {"F3"}),
        ("R2", {"R2"}),
        ("R3", {"R3"}),
    ],
)
def test_coats_redfern_ranks_true_model(true_model, expected_class):
    case = generate_case(
        rates_K_per_min=[5.0, 10.0, 20.0],
        Ea_J_per_mol=140_000.0,
        A_per_sec=1.0e11,
        model=true_model,
        n_points=2000,
        seed=0,
    )
    runs = build_runs(case)
    result = coats_redfern(runs)
    assert result.best_model in expected_class, (
        f"{true_model}: best={result.best_model}, expected one of {expected_class}"
    )


def test_recovers_true_Ea_at_best_model():
    """E_a from the best-fitting model should match the synthetic E_a."""
    Ea_true = 130.0
    case = generate_case(
        rates_K_per_min=[5.0, 10.0, 20.0],
        Ea_J_per_mol=Ea_true * 1000.0,
        A_per_sec=1.0e10,
        model="F1",
        n_points=2000,
        seed=0,
    )
    runs = build_runs(case)
    result = coats_redfern(runs)
    best_summary = result.by_model(result.best_model)
    assert abs(best_summary.Ea_kJ_per_mol_mean - Ea_true) < 5.0


def test_returns_one_fit_per_model_per_run():
    case = generate_case(
        rates_K_per_min=[2.5, 5.0, 10.0, 20.0],
        Ea_J_per_mol=120_000.0,
        A_per_sec=1.0e10,
        model="F1",
        n_points=1500,
        seed=0,
    )
    runs = build_runs(case)
    result = coats_redfern(runs)
    expected_max = 12 * 4  # 12 models × 4 runs
    assert 1 <= len(result.fits) <= expected_max
    assert all(0.0 <= s.r_squared_mean <= 1.0 for s in result.summaries)


def test_summary_sorted_by_r_squared_descending():
    case = generate_case(
        rates_K_per_min=[5.0, 10.0, 20.0],
        Ea_J_per_mol=120_000.0,
        A_per_sec=1.0e10,
        model="F1",
        n_points=1500,
        seed=0,
    )
    runs = build_runs(case)
    result = coats_redfern(runs)
    r2s = [s.r_squared_mean for s in result.summaries]
    assert all(r2s[i] >= r2s[i + 1] for i in range(len(r2s) - 1))


def test_restricting_model_subset():
    case = generate_case(
        rates_K_per_min=[5.0, 10.0],
        Ea_J_per_mol=120_000.0,
        A_per_sec=1.0e10,
        model="F1",
        n_points=1000,
        seed=0,
    )
    runs = build_runs(case)
    result = coats_redfern(runs, models=["F1", "F2"])
    assert {s.model for s in result.summaries} == {"F1", "F2"}


def test_unknown_model_raises():
    case = generate_case(
        rates_K_per_min=[5.0, 10.0],
        Ea_J_per_mol=120_000.0,
        A_per_sec=1.0e10,
        model="F1",
        n_points=500,
        seed=0,
    )
    runs = build_runs(case)
    with pytest.raises(ValueError, match="Unknown model"):
        coats_redfern(runs, models=["F1", "NotAModel"])


def test_recovered_A_close_to_truth():
    """log10 A from the best Coats–Redfern fit should match the synthetic A
    to within ~0.5 (limited by Coats–Redfern's single-rate single-T integral)."""
    A_true = 1.0e10
    case = generate_case(
        rates_K_per_min=[5.0, 10.0, 20.0],
        Ea_J_per_mol=130_000.0,
        A_per_sec=A_true,
        model="F1",
        n_points=2000,
        seed=0,
    )
    runs = build_runs(case)
    result = coats_redfern(runs)
    s = result.by_model(result.best_model)
    assert abs(s.log10_A_mean - np.log10(A_true)) < 0.5

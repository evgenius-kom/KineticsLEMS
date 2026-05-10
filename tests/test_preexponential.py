"""Pre-exponential A recovery from synthetic data."""
from __future__ import annotations

import numpy as np
import pytest

from kinetics_lems.conversion import build_runs
from kinetics_lems.methods.master_plot import rank_models
from kinetics_lems.methods.preexponential import compute_A
from kinetics_lems.methods.vyazovkin import vyazovkin
from kinetics_lems.synthetic import generate_case


@pytest.mark.parametrize(
    "true_model,A_true",
    [("F1", 1.0e10), ("F2", 1.0e12), ("R3", 1.0e11)],
)
def test_recover_A_under_correct_model(true_model, A_true):
    """When the correct f(α) is provided, log10 A should match the truth to <0.2."""
    case = generate_case(
        rates_K_per_min=[5.0, 10.0, 20.0],
        Ea_J_per_mol=130_000.0,
        A_per_sec=A_true,
        model=true_model,
        n_points=2500,
        seed=0,
    )
    runs = build_runs(case)
    alphas = np.linspace(0.1, 0.9, 17)
    vya = vyazovkin(runs, alphas)
    pre = compute_A(runs, alphas, vya.Ea_J_per_mol, model=true_model)
    assert abs(pre.log10_A_median - np.log10(A_true)) < 0.2, (
        f"{true_model}: log10 A = {pre.log10_A_median:.3f}, "
        f"expected {np.log10(A_true):.3f}"
    )
    assert pre.log10_A_mad < 0.15  # tight spread on noiseless synthetic


def test_wrong_model_gives_visibly_biased_A():
    """If we pick the wrong f(α), recovered A should differ noticeably from
    the truth — this is what tells the user the model is wrong."""
    case = generate_case(
        rates_K_per_min=[5.0, 10.0, 20.0],
        Ea_J_per_mol=130_000.0,
        A_per_sec=1.0e10,
        model="R3",
        n_points=2500,
        seed=0,
    )
    runs = build_runs(case)
    alphas = np.linspace(0.1, 0.9, 17)
    vya = vyazovkin(runs, alphas)
    pre_correct = compute_A(runs, alphas, vya.Ea_J_per_mol, model="R3")
    pre_wrong = compute_A(runs, alphas, vya.Ea_J_per_mol, model="F2")
    # Wrong model should have either a biased log10 A or a wider spread.
    correct_err = abs(pre_correct.log10_A_median - 10.0)
    wrong_err = abs(pre_wrong.log10_A_median - 10.0)
    assert wrong_err > correct_err + 0.5 or pre_wrong.log10_A_mad > pre_correct.log10_A_mad + 0.5


def test_auto_picks_model_from_master_plot():
    """The runner-style chain: rank Z(α), then compute A under best model."""
    case = generate_case(
        rates_K_per_min=[5.0, 10.0, 20.0],
        Ea_J_per_mol=130_000.0,
        A_per_sec=1.0e10,
        model="R3",
        n_points=2500,
        seed=0,
    )
    runs = build_runs(case)
    alphas = np.linspace(0.05, 0.95, 19)
    ranking = rank_models(runs, alphas)
    assert ranking.best_model == "R3"

    vya_alphas = np.linspace(0.1, 0.9, 17)
    vya = vyazovkin(runs, vya_alphas)
    pre = compute_A(runs, vya_alphas, vya.Ea_J_per_mol, model=ranking.best_model)
    assert abs(pre.log10_A_median - 10.0) < 0.3

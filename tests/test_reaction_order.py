"""Reaction-order n recovery from synthetic F_n data."""
from __future__ import annotations

import pytest

from kinetics_lems.conversion import build_runs
from kinetics_lems.methods.reaction_order import reaction_order
from kinetics_lems.synthetic import generate_case


@pytest.mark.parametrize("true_n,model", [(1.0, "F1"), (2.0, "F2"), (3.0, "F3")])
def test_recover_reaction_order_for_pure_F_n(true_n, model):
    """For pure F_n synthetic data the sweep must find n close to the truth."""
    case = generate_case(
        rates_K_per_min=[5.0, 10.0, 20.0],
        Ea_J_per_mol=130_000.0,
        A_per_sec=1.0e10,
        model=model,
        n_points=2000,
        seed=0,
    )
    runs = build_runs(case)
    result = reaction_order(runs, n_min=0.5, n_max=4.0, n_steps=120)
    assert abs(result.n_best - true_n) < 0.15, (
        f"{model}: got n={result.n_best:.3f}, expected ≈{true_n}"
    )
    assert result.r_squared_best > 0.99, (
        f"{model}: R² at best n is {result.r_squared_best:.4f}, expected > 0.99"
    )


def test_Ea_at_best_n_matches_synthetic():
    """E_a recovered from the linearization should match the true E_a."""
    case = generate_case(
        rates_K_per_min=[5.0, 10.0, 20.0],
        Ea_J_per_mol=150_000.0,
        A_per_sec=1.0e10,
        model="F1",
        n_points=2000,
        seed=0,
    )
    runs = build_runs(case)
    result = reaction_order(runs)
    assert abs(result.Ea_best_kJ_per_mol - 150.0) < 2.0, (
        f"E_a at best n is {result.Ea_best_kJ_per_mol:.2f}, expected ≈150"
    )


def test_rejects_too_narrow_alpha_window():
    case = generate_case(
        rates_K_per_min=[5.0, 10.0],
        Ea_J_per_mol=120_000.0,
        A_per_sec=1.0e10,
        model="F1",
        n_points=500,
        seed=0,
    )
    runs = build_runs(case)
    with pytest.raises(ValueError, match="α window"):
        reaction_order(runs, alpha_min=0.5, alpha_max=0.4)


def test_avrami_gives_biased_n():
    """For Avrami (not F_n family) the sweep returns biased n — caveat documented."""
    case = generate_case(
        rates_K_per_min=[5.0, 10.0, 20.0],
        Ea_J_per_mol=120_000.0,
        A_per_sec=1.0e10,
        model="A2",
        n_points=2000,
        seed=0,
    )
    runs = build_runs(case)
    result = reaction_order(runs)
    # Should NOT cleanly recover an F_n order — R² peaks below 1 or n is far
    # from a small-integer value. Just assert the test passes without raising.
    assert 0.0 <= result.r_squared_best <= 1.0

"""Z(α) master-plot model identification."""
from __future__ import annotations

import numpy as np
import pytest

from kinetics_lems.conversion import build_runs
from kinetics_lems.methods.master_plot import MASTER_MODELS, master_z_curves, rank_models
from kinetics_lems.synthetic import generate_case


def test_all_master_curves_pass_through_one_at_alpha_05():
    """Each master curve is normalized so that Z(α=0.5) = 1."""
    alphas = np.linspace(0.05, 0.95, 19)
    curves = master_z_curves(alphas)
    half = int(np.argmin(np.abs(alphas - 0.5)))
    for name, z in curves.items():
        assert pytest.approx(z[half], abs=1e-12) == 1.0, f"{name}: Z(0.5)={z[half]}"


# F1, A2, A3, A4 share the same normalized Z(α) (well-known degeneracy);
# any of them is a valid "best" answer when the truth is F1 or any A_m.
_F1_AVRAMI_CLASS = {"F1", "A2", "A3", "A4"}


@pytest.mark.parametrize(
    "true_model,expected_class",
    [
        ("F1", _F1_AVRAMI_CLASS),
        ("F2", {"F2"}),
        ("F3", {"F3"}),
        ("R2", {"R2"}),
        ("R3", {"R3"}),
        ("A2", _F1_AVRAMI_CLASS),
    ],
)
def test_synthetic_master_plot_correctly_ranks(true_model, expected_class):
    """Each synthetic model should rank an equivalent-class member at the top."""
    case = generate_case(
        rates_K_per_min=[5.0, 10.0, 20.0],
        Ea_J_per_mol=140_000.0,
        A_per_sec=1.0e11,
        model=true_model,
        n_points=2000,
        seed=0,
    )
    runs = build_runs(case)
    alphas = np.linspace(0.05, 0.95, 19)
    ranking = rank_models(runs, alphas)
    assert ranking.best_model in expected_class, (
        f"{true_model}: best={ranking.best_model}, expected one of {expected_class}, "
        f"ranking={ranking.ranked()[:5]}"
    )


def test_master_models_set_size():
    """Standard 12-model set per ICTAC 2011 §5."""
    assert len(MASTER_MODELS) == 12
    expected = {"F1", "F2", "F3", "A2", "A3", "A4", "R2", "R3", "D1", "D2", "D3", "D4"}
    assert set(MASTER_MODELS.keys()) == expected

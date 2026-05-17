"""Sestak-Berggren and Prout-Tompkins fits."""
from __future__ import annotations

import numpy as np

from kinetics_lems.conversion import build_runs
from kinetics_lems.methods import fit_prout_tompkins, fit_sestak_berggren
from kinetics_lems.synthetic import generate_case


def _runs(model: str):
    case = generate_case(
        rates_K_per_min=[2.5, 5.0, 10.0, 20.0],
        Ea_J_per_mol=120_000.0,
        A_per_sec=1.0e10,
        model=model,
    )
    return build_runs(case)


def test_prout_tompkins_fit_recovers_F1_shape() -> None:
    """Synthetic F1 has f(α) = (1−α)^1; PT should fit m≈0, n≈1."""
    runs = _runs("F1")
    alphas = np.linspace(0.05, 0.95, 19)
    fit = fit_prout_tompkins(runs, alphas)
    # PT degenerates to F1 when m → 0 and n → 1.
    assert fit.parameters["n"] > 0.5
    assert fit.r_squared > 0.95
    assert fit.rms < 0.5


def test_sestak_berggren_fit_high_r_squared() -> None:
    runs = _runs("F1")
    alphas = np.linspace(0.05, 0.95, 19)
    fit = fit_sestak_berggren(runs, alphas)
    assert fit.r_squared > 0.95


def test_empirical_fit_model_is_pluggable_in_z_curve_pipeline() -> None:
    """The returned ReactionModel must expose f, g, z without raising."""
    runs = _runs("F1")
    alphas = np.linspace(0.05, 0.95, 19)
    fit = fit_prout_tompkins(runs, alphas)
    alpha = np.linspace(0.1, 0.9, 81)
    f_vals = fit.model.f(alpha)
    g_vals = fit.model.g(alpha)
    z_vals = fit.model.z(alpha)
    assert np.all(np.isfinite(f_vals))
    assert np.all(np.isfinite(g_vals))
    assert np.all(np.isfinite(z_vals))


def test_runner_emits_empirical_fits_when_enabled() -> None:
    from pathlib import Path

    from kinetics_lems.config import EmpiricalModelsConfig, load_config
    from kinetics_lems.io import load_case
    from kinetics_lems.runner import run_analysis

    case_path = (
        Path(__file__).resolve().parent.parent / "examples" / "synthetic" / "F1_120kJ"
    )
    if not case_path.is_dir():
        import pytest

        pytest.skip("synthetic case missing")
    base = load_config()
    # Swap in an enabled empirical_models config.
    cfg = type(base)(
        conversion=base.conversion,
        enabled_methods=base.enabled_methods,
        friedman=base.friedman,
        vyazovkin=base.vyazovkin,
        vyazovkin_aic=base.vyazovkin_aic,
        preexponential=base.preexponential,
        multistep=base.multistep,
        reaction_order=base.reaction_order,
        coats_redfern=base.coats_redfern,
        uncertainty=base.uncertainty,
        consistency=base.consistency,
        empirical_models=EmpiricalModelsConfig(
            enable_prout_tompkins=True, enable_sestak_berggren=True
        ),
        lifetime=base.lifetime,
        output=base.output,
    )
    results = run_analysis(load_case(case_path), cfg)
    assert "prout_tompkins" in results.empirical_fits
    assert "sestak_berggren" in results.empirical_fits

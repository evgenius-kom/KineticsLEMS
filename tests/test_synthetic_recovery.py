"""End-to-end test: synthetic data → all methods recover the true E_a."""
from __future__ import annotations

import numpy as np
import pytest

from kinetics_lems.config import load_config
from kinetics_lems.runner import run_analysis
from kinetics_lems.synthetic import generate_case


@pytest.mark.parametrize("Ea_kJ", [80.0, 120.0, 180.0])
def test_all_methods_recover_first_order(Ea_kJ: float):
    case = generate_case(
        rates_K_per_min=[2.5, 5.0, 10.0, 20.0],
        Ea_J_per_mol=Ea_kJ * 1000.0,
        A_per_sec=1.0e10,
        model="F1",
        n_points=2000,
        seed=0,
    )
    cfg = load_config()
    res = run_analysis(case, cfg)

    # Tolerances: differential & advanced methods are exact in noise-free synthetic;
    # OFW relies on Doyle's approximation so allow a few percent.
    tolerances = {
        "friedman": 1.0,        # kJ/mol
        "kas": 1.5,
        "ofw": 5.0,
        "vyazovkin": 1.0,
        "vyazovkin_aic": 1.0,
    }
    for name, tol in tolerances.items():
        Ea_recovered = float(np.nanmean(res.isoconversional[name].Ea_kJ_per_mol))
        assert abs(Ea_recovered - Ea_kJ) < tol, (
            f"{name}: expected {Ea_kJ}, got {Ea_recovered:.2f} (tol {tol})"
        )

    assert res.kissinger is not None
    assert abs(res.kissinger.Ea_kJ_per_mol - Ea_kJ) < 2.0


def test_ea_is_constant_for_single_step_reaction():
    """For a true single-step reaction E(α) should be flat."""
    case = generate_case(
        rates_K_per_min=[2.0, 5.0, 10.0, 20.0],
        Ea_J_per_mol=100_000.0,
        A_per_sec=1.0e9,
        model="F1",
        n_points=2000,
        seed=0,
    )
    cfg = load_config()
    res = run_analysis(case, cfg)
    for name in ("friedman", "vyazovkin", "vyazovkin_aic"):
        Ea = res.isoconversional[name].Ea_kJ_per_mol
        Ea = Ea[~np.isnan(Ea)]
        # Standard deviation across α should be much smaller than the mean.
        assert Ea.std() < 1.0, f"{name}: E(α) not flat, std={Ea.std():.3f}"


def test_noise_robustness():
    """With moderate noise (3%) we still recover E_a within a few %."""
    case = generate_case(
        rates_K_per_min=[2.5, 5.0, 10.0, 20.0],
        Ea_J_per_mol=120_000.0,
        A_per_sec=1.0e10,
        model="F1",
        n_points=3000,
        noise_std=0.03,
        seed=42,
    )
    cfg = load_config()
    res = run_analysis(case, cfg)
    # Differential methods amplify noise; keep tolerance loose for Friedman.
    assert abs(np.nanmean(res.isoconversional["friedman"].Ea_kJ_per_mol) - 120.0) < 6.0
    assert abs(np.nanmean(res.isoconversional["vyazovkin"].Ea_kJ_per_mol) - 120.0) < 4.0
    assert abs(np.nanmean(res.isoconversional["vyazovkin_aic"].Ea_kJ_per_mol) - 120.0) < 4.0


def test_different_reaction_models_recover_ea():
    """Methods should give correct E_a regardless of f(α)."""
    cfg = load_config()
    for model in ("F1", "F2", "R3"):
        case = generate_case(
            rates_K_per_min=[5.0, 10.0, 20.0],
            Ea_J_per_mol=150_000.0,
            A_per_sec=1.0e12,
            model=model,
            n_points=3000,
            seed=0,
        )
        res = run_analysis(case, cfg)
        Ea_friedman = float(np.nanmean(res.isoconversional["friedman"].Ea_kJ_per_mol))
        Ea_vyaz = float(np.nanmean(res.isoconversional["vyazovkin"].Ea_kJ_per_mol))
        assert abs(Ea_friedman - 150.0) < 1.5, f"{model}: Friedman={Ea_friedman}"
        assert abs(Ea_vyaz - 150.0) < 1.5, f"{model}: Vyazovkin={Ea_vyaz}"

"""New synthetic datasets — sanity checks + pipeline smoke tests."""
from __future__ import annotations

import numpy as np

from kinetics_lems.conversion import build_runs
from kinetics_lems.methods import assess_endpoints, vyazovkin
from kinetics_lems.synthetic import (
    generate_arbitrary_program_case,
    generate_daem_gaussian_case,
    generate_two_parallel_case,
)


def test_two_parallel_case_produces_multi_step_E_alpha() -> None:
    case = generate_two_parallel_case(
        rates_K_per_min=[2.5, 5.0, 10.0, 20.0],
        E1_J_per_mol=100_000.0,
        A1_per_sec=1.0e9,
        weight1=0.5,
        E2_J_per_mol=180_000.0,
        A2_per_sec=1.0e12,
    )
    runs = build_runs(case)
    alphas = np.linspace(0.1, 0.9, 17)
    iso = vyazovkin(runs, alphas)
    # Two channels of different E → E(α) must vary noticeably across α.
    finite = iso.Ea_J_per_mol[np.isfinite(iso.Ea_J_per_mol)]
    spread = finite.max() - finite.min()
    assert spread > 30_000.0, f"expected E spread > 30 kJ/mol, got {spread/1000:.1f}"
    rel = assess_endpoints(iso)
    assert any("varies" in w for w in rel.warnings)


def test_daem_gaussian_case_E_alpha_is_smooth_and_centered() -> None:
    case = generate_daem_gaussian_case(
        rates_K_per_min=[2.5, 5.0, 10.0, 20.0],
        E_mean_J_per_mol=150_000.0,
        E_sigma_J_per_mol=10_000.0,
        A_per_sec=1.0e11,
    )
    runs = build_runs(case)
    alphas = np.linspace(0.1, 0.9, 17)
    iso = vyazovkin(runs, alphas)
    finite = iso.Ea_J_per_mol[np.isfinite(iso.Ea_J_per_mol)]
    # Vyazovkin sees the average E across the Gaussian; should be within
    # ±20% of Ē=150 kJ/mol for a 10 kJ/mol σ.
    assert 120_000.0 < np.mean(finite) < 180_000.0


def test_arbitrary_program_case_carries_three_column_wave() -> None:
    case = generate_arbitrary_program_case(program="modulated")
    wave = next(iter(case.waves.values()))
    assert wave.has_recorded_time
    assert wave.t_seconds is not None and wave.t_seconds[0] == 0.0
    runs = build_runs(case)
    # α(T) and dα/dt(α) are non-trivial.
    assert runs[0].alpha[-1] > 0.5
    assert runs[0].dalpha_dt.max() > 0

"""Endpoint reliability and single-step diagnostic."""
from __future__ import annotations

import numpy as np

from kinetics_lems.methods import IsoconversionalResult, assess_endpoints


def _iso(alpha: np.ndarray, E_kJ: np.ndarray) -> IsoconversionalResult:
    return IsoconversionalResult(
        method="vyazovkin",
        alpha=alpha,
        Ea_J_per_mol=E_kJ * 1000.0,
        intercept=np.zeros_like(alpha),
        r_squared=np.ones_like(alpha),
    )


def test_flat_curve_no_warnings() -> None:
    alpha = np.linspace(0.05, 0.95, 19)
    E = np.full_like(alpha, 120.0)
    rel = assess_endpoints(_iso(alpha, E))
    assert rel.warnings == []
    assert rel.flatness_in_core < 1e-9


def test_strong_variation_in_core_triggers_multistep_warning() -> None:
    alpha = np.linspace(0.05, 0.95, 19)
    E = 100.0 + 50.0 * alpha  # 50 kJ/mol spread over 0..1 → ~50% in the core
    rel = assess_endpoints(_iso(alpha, E))
    assert any("varies" in w for w in rel.warnings)
    assert rel.flatness_in_core > 0.10


def test_low_alpha_tail_nan_triggers_tail_warning() -> None:
    alpha = np.linspace(0.05, 0.95, 19)
    E = np.full_like(alpha, 120.0)
    E[alpha < 0.1] = np.nan  # all low-α tail invalid
    rel = assess_endpoints(_iso(alpha, E), alpha_low=0.1, alpha_high=0.9)
    assert any("low-α tail" in w for w in rel.warnings)


def test_high_alpha_tail_nan_triggers_tail_warning() -> None:
    alpha = np.linspace(0.05, 0.95, 19)
    E = np.full_like(alpha, 120.0)
    E[alpha > 0.9] = np.nan
    rel = assess_endpoints(_iso(alpha, E), alpha_low=0.1, alpha_high=0.9)
    assert any("high-α tail" in w for w in rel.warnings)


def test_runner_emits_endpoint_reliability_for_every_iso_method() -> None:
    from pathlib import Path

    from kinetics_lems.config import load_config
    from kinetics_lems.io import load_case
    from kinetics_lems.runner import run_analysis

    case_path = (
        Path(__file__).resolve().parent.parent / "examples" / "synthetic" / "F1_120kJ"
    )
    if not case_path.is_dir():
        import pytest

        pytest.skip("synthetic case missing")
    config = load_config()
    results = run_analysis(load_case(case_path), config)
    assert set(results.endpoint_reliability) == set(results.isoconversional)

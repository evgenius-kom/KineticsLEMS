"""Config validation + numerical-robustness checks."""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from kinetics_lems.config import ConversionConfig, load_config
from kinetics_lems.conversion import ConversionRun
from kinetics_lems.methods.friedman import friedman

# ---------- Config validation ----------

def test_default_config_loads():
    cfg = load_config()
    assert cfg.conversion.step > 0
    assert "vyazovkin" in cfg.enabled_methods


@pytest.mark.parametrize(
    "bad_kwargs",
    [
        {"min": -0.1, "max": 0.9},   # min < 0
        {"min": 0.5, "max": 0.4},    # min > max
        {"min": 0.05, "max": 1.0},   # max == 1
        {"min": 0.05, "max": 0.95, "step": 0.0},  # zero step
    ],
)
def test_bad_conversion_config_rejected(bad_kwargs):
    cfg = ConversionConfig(**bad_kwargs)
    with pytest.raises(ValueError):
        cfg.grid()


def test_unknown_method_in_toml_rejected(tmp_path: Path):
    bad = tmp_path / "bad.toml"
    bad.write_text(
        '[methods]\nenabled = ["friedman", "magic_method"]\n', encoding="utf-8"
    )
    with pytest.raises(ValueError, match="Unknown method"):
        load_config(bad)


# ---------- Friedman robustness against non-positive dα/dt ----------

def test_friedman_returns_nan_on_zero_rate(monkeypatch):
    """When dα/dt is 0 (or negative) at some α, ln() is undefined.

    Friedman should return NaN at that α and continue without raising.
    """
    # Two synthetic runs whose dα/dt has a zero entry at α = 0.5.
    T1 = np.array([300.0, 350.0, 400.0])
    T2 = np.array([320.0, 370.0, 420.0])
    alpha = np.array([0.0, 0.5, 1.0])
    rate1 = np.array([1e-3, 0.0, 1e-3])  # zero at α=0.5
    rate2 = np.array([2e-3, 1e-3, 2e-3])
    runs = [
        ConversionRun(rate_K_per_min=5.0, temperature=T1, alpha=alpha, dalpha_dt=rate1),
        ConversionRun(rate_K_per_min=10.0, temperature=T2, alpha=alpha, dalpha_dt=rate2),
    ]

    res = friedman(runs, np.array([0.25, 0.5, 0.75]))
    # α=0.5 hits the zero rate in run 1 → NaN there, finite at others.
    assert np.isnan(res.Ea_J_per_mol[1])
    assert np.isfinite(res.Ea_J_per_mol[0])
    assert np.isfinite(res.Ea_J_per_mol[2])

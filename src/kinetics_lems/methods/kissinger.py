"""Classical Kissinger (1957) — single Ea from peak temperatures.

Plot ln(β / T_p²) vs 1/T_p across runs (β in K/s, T_p = peak temperature).
Slope = -E_a / R; intercept ≈ ln(A·R/E_a) (assuming first-order).
"""
from __future__ import annotations

import numpy as np

from ..constants import R_GAS, SEC_PER_MIN
from ..conversion import ConversionRun
from .common import KissingerResult, linear_regression


def kissinger(runs: list[ConversionRun]) -> KissingerResult:
    if len(runs) < 2:
        raise ValueError("Kissinger needs at least 2 heating rates")

    Tp = np.array([_peak_temperature(r) for r in runs])
    betas_per_sec = np.array([r.rate_K_per_min / SEC_PER_MIN for r in runs])

    x = 1.0 / Tp
    y = np.log(betas_per_sec / Tp**2)
    slope, intercept, r2 = linear_regression(x, y)

    Ea = -slope * R_GAS
    A = np.exp(intercept) * Ea / R_GAS  # assumes first-order, n=1
    return KissingerResult(
        Ea_J_per_mol=Ea,
        pre_exponential=float(A),
        r_squared=r2,
        Tp_K=Tp,
        rates_K_per_min=np.array([r.rate_K_per_min for r in runs]),
    )


def _peak_temperature(run: ConversionRun) -> float:
    """Temperature at maximum dα/dt (corresponds to the DSC/DTG peak)."""
    return float(run.temperature[int(np.argmax(run.dalpha_dt))])


__all__ = ["kissinger"]

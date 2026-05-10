"""Kissinger–Akahira–Sunose (KAS) — integral isoconversional method.

For each fixed α plot ln(β / T²) vs 1/T(α) across runs.
β is in K/s. Slope = -E_a/R.
"""
from __future__ import annotations

import numpy as np

from ..constants import R_GAS, SEC_PER_MIN
from ..conversion import ConversionRun, temperature_at_conversion
from .common import IsoconversionalResult, linear_regression


def kas(runs: list[ConversionRun], alphas: np.ndarray) -> IsoconversionalResult:
    if len(runs) < 2:
        raise ValueError("KAS needs at least 2 heating rates")

    betas_per_sec = np.array([r.rate_K_per_min / SEC_PER_MIN for r in runs])
    T_table = np.vstack([temperature_at_conversion(r, alphas) for r in runs])  # (n_runs, n_α)

    Ea = np.empty_like(alphas, dtype=float)
    intercepts = np.empty_like(alphas, dtype=float)
    r2 = np.empty_like(alphas, dtype=float)

    for i in range(alphas.size):
        T = T_table[:, i]
        x = 1.0 / T
        y = np.log(betas_per_sec / T**2)
        slope, intercept, r = linear_regression(x, y)
        Ea[i] = -slope * R_GAS
        intercepts[i] = intercept
        r2[i] = r

    return IsoconversionalResult("KAS", alphas, Ea, intercepts, r2)


__all__ = ["kas"]

"""Ozawa–Flynn–Wall (OFW) — integral isoconversional method using Doyle's approximation.

For each fixed α plot ln(β) vs 1/T(α) across runs.
β is in K/s. Slope = -1.052 · E_a / R.
"""
from __future__ import annotations

import numpy as np

from ..constants import DOYLE_CORRECTION, R_GAS, SEC_PER_MIN
from ..conversion import ConversionRun, temperature_at_conversion
from .common import IsoconversionalResult, linear_regression


def ofw(runs: list[ConversionRun], alphas: np.ndarray) -> IsoconversionalResult:
    if len(runs) < 2:
        raise ValueError("OFW needs at least 2 heating rates")

    betas_per_sec = np.array([r.rate_K_per_min / SEC_PER_MIN for r in runs])
    T_table = np.vstack([temperature_at_conversion(r, alphas) for r in runs])

    Ea = np.empty_like(alphas, dtype=float)
    intercepts = np.empty_like(alphas, dtype=float)
    r2 = np.empty_like(alphas, dtype=float)

    for i in range(alphas.size):
        x = 1.0 / T_table[:, i]
        y = np.log(betas_per_sec)
        slope, intercept, r = linear_regression(x, y)
        Ea[i] = -slope * R_GAS / DOYLE_CORRECTION
        intercepts[i] = intercept
        r2[i] = r

    return IsoconversionalResult("OFW", alphas, Ea, intercepts, r2)


__all__ = ["ofw"]

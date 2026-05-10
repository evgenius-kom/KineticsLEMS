"""Friedman (1964) — differential isoconversional method.

For each fixed α plot ln(dα/dt) vs 1/T(α) across runs (varying β).
The slope gives -E_a/R; the intercept estimates ln{A·f(α)}.
"""
from __future__ import annotations

import numpy as np

from ..constants import R_GAS
from ..conversion import ConversionRun, dalpha_dt_at_conversion, temperature_at_conversion
from .common import IsoconversionalResult, linear_regression


def friedman(runs: list[ConversionRun], alphas: np.ndarray) -> IsoconversionalResult:
    if len(runs) < 2:
        raise ValueError("Friedman needs at least 2 heating rates")

    Ea = np.empty_like(alphas, dtype=float)
    intercepts = np.empty_like(alphas, dtype=float)
    r2 = np.empty_like(alphas, dtype=float)

    T_table = np.vstack([temperature_at_conversion(r, alphas) for r in runs])  # shape (n_runs, n_α)
    rate_table = np.vstack([dalpha_dt_at_conversion(r, alphas) for r in runs])

    # ln(dα/dt) is undefined for non-positive rates. On real DSC data the
    # tails near α → 0 and α → 1 can dip to zero or below from noise;
    # let those points become NaN and propagate without spamming warnings.
    with np.errstate(divide="ignore", invalid="ignore"):
        for i in range(alphas.size):
            x = 1.0 / T_table[:, i]
            y = np.log(rate_table[:, i])
            if not np.all(np.isfinite(y)):
                Ea[i] = np.nan
                intercepts[i] = np.nan
                r2[i] = np.nan
                continue
            slope, intercept, r = linear_regression(x, y)
            Ea[i] = -slope * R_GAS
            intercepts[i] = intercept
            r2[i] = r

    return IsoconversionalResult("Friedman", alphas, Ea, intercepts, r2)


__all__ = ["friedman"]

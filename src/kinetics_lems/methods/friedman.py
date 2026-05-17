"""Friedman (1964) — differential isoconversional method.

For each fixed α plot ln(dα/dt) vs 1/T(α) across runs (varying β).
The slope gives -E_a/R; the intercept estimates ln{A·f(α)}.

Optional pre-smoothing
----------------------
Differentiation amplifies noise; Friedman is the only method we have that
operates on dα/dt directly (the integral methods integrate noise out).
When ``smooth_window`` is set, dα/dt is filtered with a Savitzky-Golay
smoother *per heating-rate run*, before α-resampling. This dramatically
improves Friedman robustness on noisy DSC traces (ICTAC 2020 §3.4
mentions noise pre-treatment as a routine pre-step).

The smoother is **opt-in** — by default Friedman runs on the raw rates so
that synthetic and reference results stay reproducible to past versions.
"""
from __future__ import annotations

import numpy as np
from scipy.signal import savgol_filter

from ..constants import R_GAS
from ..conversion import ConversionRun, dalpha_dt_at_conversion, temperature_at_conversion
from .common import IsoconversionalResult, linear_regression


def friedman(
    runs: list[ConversionRun],
    alphas: np.ndarray,
    *,
    smooth_window: int | None = None,
    smooth_poly: int = 3,
) -> IsoconversionalResult:
    """Differential isoconversional E(α).

    Parameters
    ----------
    runs, alphas :
        Conversion curves and target α grid.
    smooth_window :
        Length of the Savitzky-Golay window applied to dα/dt **before**
        α-resampling, per run. Must be odd and ≥ ``smooth_poly + 2``.
        ``None`` (default) disables smoothing.
    smooth_poly :
        Polynomial order of the SavGol fit. Default 3 — preserves peak shape
        well; 2 is more aggressive, 4 keeps more high-frequency content.
    """
    if len(runs) < 2:
        raise ValueError("Friedman needs at least 2 heating rates")
    if smooth_window is not None:
        if smooth_window % 2 == 0 or smooth_window < smooth_poly + 2:
            raise ValueError(
                f"smooth_window must be odd and >= smooth_poly+2 "
                f"(got window={smooth_window}, poly={smooth_poly})"
            )

    Ea = np.empty_like(alphas, dtype=float)
    intercepts = np.empty_like(alphas, dtype=float)
    r2 = np.empty_like(alphas, dtype=float)

    if smooth_window is None:
        smoothed_runs = runs
    else:
        smoothed_runs = [_smooth_run(r, smooth_window, smooth_poly) for r in runs]

    T_table = np.vstack(
        [temperature_at_conversion(r, alphas) for r in smoothed_runs]
    )  # shape (n_runs, n_α)
    rate_table = np.vstack(
        [dalpha_dt_at_conversion(r, alphas) for r in smoothed_runs]
    )

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


def _smooth_run(run: ConversionRun, window: int, poly: int) -> ConversionRun:
    """Apply a Savitzky-Golay smoother to dα/dt; keep α and T as-is."""
    if window > run.dalpha_dt.size:
        # Shrink to the largest valid odd window for short runs.
        window = max(poly + 3, run.dalpha_dt.size - (run.dalpha_dt.size + 1) % 2)
        if window % 2 == 0:
            window -= 1
        if window < poly + 2:
            return run  # too few points to smooth; pass through
    smoothed = savgol_filter(run.dalpha_dt, window_length=window, polyorder=poly)
    return ConversionRun(
        rate_K_per_min=run.rate_K_per_min,
        temperature=run.temperature,
        alpha=run.alpha,
        dalpha_dt=smoothed,
    )


__all__ = ["friedman"]

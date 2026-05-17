"""Compute α(T) from a DSC peak (or analogous curve) and extract T(α) per heating rate.

Workflow (see docs/ALGORITHMS.md §A):
    1. Subtract a linear baseline from the peak.
    2. Cumulatively integrate y vs x (trapezoid) — gives total heat / mass change.
    3. Normalize by the total integral to obtain α ∈ [0, 1].
    4. Interpolate T(α) at user-requested conversions.

3-column support
----------------
When a wave carries ``t_seconds`` (3-column input read by
:func:`kinetics_lems.io.wave_reader.read_wave`), the integration is done
against the recorded time, *not* against T under a linear-heating
assumption. This makes the α(T) and dα/dt curves correct for arbitrary
T(t) programs (modulated DSC, T-jump, fast-cycling FSC).
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .models import CaseData, Wave


@dataclass(frozen=True)
class ConversionRun:
    """α(T) curve plus dα/dt for a single experiment (one heating rate)."""

    rate_K_per_min: float
    """Heating rate β in K/min as supplied by the user."""

    temperature: np.ndarray  # K
    alpha: np.ndarray  # dimensionless, monotone non-decreasing in [0, 1]
    dalpha_dt: np.ndarray  # 1/s — derivative of α with respect to time


def build_runs(case: CaseData) -> list[ConversionRun]:
    runs: list[ConversionRun] = []
    for rate_K_per_min, wave in case.waves.items():
        runs.append(_build_run(rate_K_per_min, wave))
    runs.sort(key=lambda r: r.rate_K_per_min)
    return runs


def _build_run(rate_K_per_min: float, peak: Wave) -> ConversionRun:
    if rate_K_per_min <= 0:
        raise ValueError(f"Heating rate must be positive, got {rate_K_per_min}")

    peak_corrected = peak.subtract_baseline()

    if peak_corrected.has_recorded_time:
        return _build_run_from_time(rate_K_per_min, peak_corrected)
    return _build_run_linear_heating(rate_K_per_min, peak_corrected)


def _build_run_linear_heating(rate_K_per_min: float, peak: Wave) -> ConversionRun:
    """Classical 2-column path: integrate vs T, assume β = const."""
    T = peak.x.astype(float)
    y = peak.y.astype(float)

    if not np.all(np.diff(T) > 0):
        order = np.argsort(T)
        T = T[order]
        y = y[order]

    integral = _cumulative_trapezoid(y, T)
    total = integral[-1]
    if total <= 0:
        # Endothermic / mass-loss — flip sign so α grows from 0 to 1.
        integral = -integral
        total = integral[-1]
    if total <= 0:
        raise ValueError("Degenerate peak: integral is zero — cannot normalize α")
    alpha = np.clip(integral / total, 0.0, 1.0)

    # dα/dt: chain rule. dα/dT = y / total (since integral'(T) = y), and dT/dt = β (K/s).
    beta_K_per_sec = rate_K_per_min / 60.0
    dalpha_dt = (y / total) * beta_K_per_sec

    return ConversionRun(
        rate_K_per_min=float(rate_K_per_min),
        temperature=T,
        alpha=alpha,
        dalpha_dt=dalpha_dt,
    )


def _build_run_from_time(rate_K_per_min: float, peak: Wave) -> ConversionRun:
    """3-column path: integrate vs t, use the recorded T(t) directly.

    Treats ``y`` as proportional to the *rate* of the kinetic signal
    (heat flow, mass-loss rate, dα/dt). The integral ``∫ y dt`` is the
    total extent and α(t) = ∫₀^t y dt' / ∫ y dt is monotonic on [0, 1].
    α(T) is then obtained by interpolation along the recorded T(t).

    ``rate_K_per_min`` is retained on the result as a *nominal* β (for
    sorting and label purposes) but plays no role in this math —
    isoconversional methods that depend on β explicitly (KAS, OFW,
    classical Vyazovkin) still use the rate from settings.json, which is
    acceptable when T(t) is *almost* linear with small modulation. For
    strongly non-linear programs, prefer :func:`vyazovkin_aic`.
    """
    if peak.t_seconds is None:  # pragma: no cover — guarded by has_recorded_time
        raise AssertionError("3-column path entered without recorded time")
    t = np.asarray(peak.t_seconds, dtype=float)
    T = peak.x.astype(float)
    y = peak.y.astype(float)

    integral = _cumulative_trapezoid(y, t)
    total = integral[-1]
    if total <= 0:
        integral = -integral
        total = integral[-1]
    if total <= 0:
        raise ValueError("Degenerate 3-column wave: integral of y is zero")
    alpha = np.clip(integral / total, 0.0, 1.0)
    dalpha_dt = y / total  # 1/s; absolute calibration of y cancels

    return ConversionRun(
        rate_K_per_min=float(rate_K_per_min),
        temperature=T,
        alpha=alpha,
        dalpha_dt=dalpha_dt,
    )


def temperature_at_conversion(run: ConversionRun, alphas: np.ndarray) -> np.ndarray:
    """Return T at each requested α via linear interpolation on the monotone branch."""
    a = run.alpha
    T = run.temperature
    # Restrict to the monotone non-decreasing prefix to keep np.interp well-defined.
    mono = np.maximum.accumulate(a)
    return np.interp(alphas, mono, T)


def dalpha_dt_at_conversion(run: ConversionRun, alphas: np.ndarray) -> np.ndarray:
    """Linear interpolation of dα/dt at each requested α."""
    mono = np.maximum.accumulate(run.alpha)
    return np.interp(alphas, mono, run.dalpha_dt)


def _cumulative_trapezoid(y: np.ndarray, x: np.ndarray) -> np.ndarray:
    """Cumulative trapezoid integral with the same length as y, starting at 0."""
    out = np.zeros_like(y, dtype=float)
    out[1:] = np.cumsum(0.5 * (y[1:] + y[:-1]) * np.diff(x))
    return out


__all__ = [
    "ConversionRun",
    "build_runs",
    "temperature_at_conversion",
    "dalpha_dt_at_conversion",
]

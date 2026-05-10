"""Compute α(T) from a DSC peak (or analogous curve) and extract T(α) per heating rate.

Workflow (see docs/ALGORITHMS.md §A):
    1. Subtract a linear baseline from the peak.
    2. Cumulatively integrate y vs x (trapezoid) — gives total heat / mass change.
    3. Normalize by the total integral to obtain α ∈ [0, 1].
    4. Interpolate T(α) at user-requested conversions.
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
    T = peak_corrected.x.astype(float)
    y = peak_corrected.y.astype(float)

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

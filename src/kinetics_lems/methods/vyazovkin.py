"""Vyazovkin isoconversional methods.

Two flavors:

* :func:`vyazovkin` — classical 1996/1997 formulation, valid for linear-heating
  experiments only. Per α, find E that minimizes

      Φ(E) = sum_{i ≠ j}  I(E, T_α,i) · β_j  /  ( I(E, T_α,j) · β_i )

  with the temperature-integral approximation

      I(E, T)  ≈  (E / R) · p(x),     x = E / (R · T),
      p(x)     =  exp(-x) / x · (x² + 10x + 18) / (x³ + 12x² + 36x + 24)   (Senum–Yang).

* :func:`vyazovkin_aic` — advanced isoconversional (Vyazovkin 2001), evaluates
  the integral numerically over a small window α ∈ [α − Δα, α + Δα] using the
  actual T(t) of each run. Works for arbitrary heating programs and is more
  robust to E(α) variability.
"""
from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass

import numpy as np
from scipy.optimize import minimize_scalar

from ..constants import R_GAS, SEC_PER_MIN
from ..conversion import ConversionRun, temperature_at_conversion
from .common import IsoconversionalResult


def _minimize_bounded(phi: Callable[[float], float], bounds_J: tuple[float, float]) -> float:
    """Minimize ``phi(E)`` on a bounded 1-D interval and return the argmin in J/mol."""
    res = minimize_scalar(phi, bounds=bounds_J, method="bounded", options={"xatol": 1.0})
    return float(res.x)  # pyright: ignore[reportAttributeAccessIssue]

# ---------- Senum–Yang p(x) approximation ----------
# Accurate to better than 1e-5 for x > 20, which covers any realistic
# E_a / (R·T) in thermal-analysis kinetics (typically 20–100).

def _p_senum_yang(x: float) -> float:
    num = x * x + 10.0 * x + 18.0
    den = x * x * x + 12.0 * x * x + 36.0 * x + 24.0
    return float(np.exp(-x) / x * num / den)


# ---------- Classical Vyazovkin (linear heating only) ----------

def vyazovkin(
    runs: list[ConversionRun],
    alphas: np.ndarray,
    Ea_bracket_kJ: tuple[float, float] = (1.0, 600.0),
) -> IsoconversionalResult:
    """Classical Vyazovkin (linear heating). Returns E(α) in J/mol."""
    if len(runs) < 2:
        raise ValueError("Vyazovkin needs at least 2 heating rates")

    betas_per_sec = np.array([r.rate_K_per_min / SEC_PER_MIN for r in runs])
    T_table = np.vstack([temperature_at_conversion(r, alphas) for r in runs])  # (n_runs, n_α)

    Ea = np.empty_like(alphas, dtype=float)
    for i in range(alphas.size):
        T_alpha = T_table[:, i]
        Ea[i] = _minimize_phi_classical(T_alpha, betas_per_sec, Ea_bracket_kJ)

    intercept = np.full_like(alphas, np.nan, dtype=float)
    r2 = np.full_like(alphas, np.nan, dtype=float)
    return IsoconversionalResult("Vyazovkin", alphas, Ea, intercept, r2)


def _minimize_phi_classical(
    T_alpha: np.ndarray,
    betas_per_sec: np.ndarray,
    Ea_bracket_kJ: tuple[float, float],
) -> float:
    n = T_alpha.size

    def phi(E_J: float) -> float:
        # I(E, T) under linear heating ≈ (E / (R·β)) · p(x); the β cancels naturally
        # in the ratio below, so we keep it explicit per the docx formulation:
        # term = I_i · β_j  /  (I_j · β_i)  with I_k = (E/R) · p(x_k).
        x = E_J / (R_GAS * T_alpha)
        # The (E/R) prefactor of I(E, T_α) = (E/R)·p(x) cancels in every
        # numerator/denominator pair, so we work with p(x) directly.
        p = np.array([_p_senum_yang(xi) for xi in x])
        total = 0.0
        for i in range(n):
            for j in range(n):
                if i == j:
                    continue
                total += (p[i] * betas_per_sec[j]) / (p[j] * betas_per_sec[i])
        return total

    bounds = (Ea_bracket_kJ[0] * 1000.0, Ea_bracket_kJ[1] * 1000.0)
    return _minimize_bounded(phi, bounds)


# ---------- Advanced Isoconversional (AIC, Vyazovkin 2001) ----------

@dataclass(frozen=True)
class _RunForAIC:
    rate_K_per_min: float
    t: np.ndarray  # seconds
    T: np.ndarray  # K
    alpha: np.ndarray


def vyazovkin_aic(
    runs: list[ConversionRun],
    alphas: np.ndarray,
    delta_alpha: float = 0.02,
    Ea_bracket_kJ: tuple[float, float] = (1.0, 600.0),
) -> IsoconversionalResult:
    """Advanced isoconversional (Vyazovkin 2001).

    For each α ∈ {α₁, …, α_m} (with α_k ∈ [Δα, 1−Δα]) compute

        J(E, [α-Δα, α+Δα])  =  ∫_{t(α-Δα)}^{t(α+Δα)} exp(-E / (R·T(t))) dt,

    then minimize  Φ(E) = Σ_{i ≠ j} J_i / J_j.

    No assumption of linear heating: each run's true T(t) is used.
    """
    if len(runs) < 2:
        raise ValueError("Vyazovkin AIC needs at least 2 runs")
    if not (0.0 < delta_alpha < 0.5):
        raise ValueError("delta_alpha must be in (0, 0.5)")

    aic_runs = [_to_aic_run(r) for r in runs]

    Ea = np.empty_like(alphas, dtype=float)
    for k, a in enumerate(alphas):
        if a - delta_alpha < 0.0 or a + delta_alpha > 1.0:
            Ea[k] = np.nan
            continue
        windows = [_extract_window(run, a - delta_alpha, a + delta_alpha) for run in aic_runs]
        Ea[k] = _minimize_phi_aic(windows, Ea_bracket_kJ)

    intercept = np.full_like(alphas, np.nan, dtype=float)
    r2 = np.full_like(alphas, np.nan, dtype=float)
    return IsoconversionalResult("Vyazovkin-AIC", alphas, Ea, intercept, r2)


def _to_aic_run(run: ConversionRun) -> _RunForAIC:
    """Build a (t, T, α) triple. Time is reconstructed from T and the heating rate."""
    beta_per_sec = run.rate_K_per_min / SEC_PER_MIN
    t = (run.temperature - run.temperature[0]) / beta_per_sec
    return _RunForAIC(
        rate_K_per_min=run.rate_K_per_min,
        t=t,
        T=run.temperature,
        alpha=np.maximum.accumulate(run.alpha),
    )


def _extract_window(run: _RunForAIC, a_lo: float, a_hi: float) -> tuple[np.ndarray, np.ndarray]:
    """Return (t, T) sampled densely on [α_lo, α_hi] for one run."""
    a = run.alpha
    t = run.t
    T = run.T
    t_lo = float(np.interp(a_lo, a, t))
    t_hi = float(np.interp(a_hi, a, t))
    if t_hi <= t_lo:
        # Pathological: extremely sharp transition; fall back to two-point window.
        return np.array([t_lo, max(t_lo + 1e-9, t_hi)]), np.array(
            [float(np.interp(a_lo, a, T)), float(np.interp(a_hi, a, T))]
        )
    n_points = 64
    t_grid = np.linspace(t_lo, t_hi, n_points)
    T_grid = np.interp(t_grid, t, T)
    return t_grid, T_grid


def _minimize_phi_aic(
    windows: list[tuple[np.ndarray, np.ndarray]],
    Ea_bracket_kJ: tuple[float, float],
) -> float:
    n = len(windows)

    def J(E_J: float, t_grid: np.ndarray, T_grid: np.ndarray) -> float:
        return float(np.trapezoid(np.exp(-E_J / (R_GAS * T_grid)), t_grid))

    def phi(E_J: float) -> float:
        Js = np.array([J(E_J, *w) for w in windows])
        # Numerical floor — windows can shrink to ~0 for nearly degenerate runs.
        Js = np.maximum(Js, 1e-300)
        total = 0.0
        for i in range(n):
            for j in range(n):
                if i == j:
                    continue
                total += Js[i] / Js[j]
        return total

    bounds = (Ea_bracket_kJ[0] * 1000.0, Ea_bracket_kJ[1] * 1000.0)
    return _minimize_bounded(phi, bounds)


__all__ = ["vyazovkin", "vyazovkin_aic"]

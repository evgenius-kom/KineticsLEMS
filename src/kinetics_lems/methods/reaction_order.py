"""Reaction-order n from a linearization sweep.

For a generic n-th-order rate law f(α) = (1 − α)ⁿ, equation (1) of
:mod:`~kinetics_lems.methods.friedman` rewritten in the T-domain gives

    ln(dα/dT) + ln(β) − n · ln(1 − α)  =  ln A  −  E_a / (R·T)        (1)

i.e. once the (constant) β term is moved to the LHS, the LHS is linear
in 1/T with a rate-independent intercept ln A and slope −E_a/R for the
correct n. Pooling several heating rates is then valid — the regression
absorbs no run-dependent intercept.

Sweeping n over a grid and picking the value that maximizes pooled R²
gives an empirical n that is exact for pure F_n kinetics and a useful
proxy otherwise.

Caveat: meaningful only when the underlying kinetics actually belongs
to the F_n family. For Avrami, contracting-geometry, or diffusion
kinetics the recovered n is biased — confirm with the master-plot
ranking before quoting a number.

Reference: ICTAC 2011 §4; spreadsheet "Проверка порядка" in
`theory/Сводная таблица.xlsx`.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..constants import R_GAS
from ..conversion import ConversionRun
from .common import linear_regression


@dataclass(frozen=True)
class ReactionOrderResult:
    """Best-fit n in f(α) = (1 − α)ⁿ from the linearization sweep."""

    n_grid: np.ndarray
    """Candidate orders tried."""

    r_squared: np.ndarray
    """Pooled R² for each candidate n, aligned with ``n_grid``."""

    Ea_kJ_per_mol: np.ndarray
    """E_a recovered at each candidate n (from the slope of the best line)."""

    n_best: float
    """n with the highest R²."""

    Ea_best_kJ_per_mol: float
    """E_a at ``n_best``."""

    r_squared_best: float
    """R² at ``n_best``."""


def reaction_order(
    runs: list[ConversionRun],
    *,
    alpha_min: float = 0.1,
    alpha_max: float = 0.9,
    n_min: float = 0.1,
    n_max: float = 4.0,
    n_steps: int = 80,
) -> ReactionOrderResult:
    """Sweep n and find the value that best linearizes equation (1).

    Parameters
    ----------
    runs:
        Conversion runs at different heating rates. Pooled across runs
        — the linearization absorbs the (A/β) constant into the
        intercept, so different β just contribute different intercept
        clouds without biasing the slope.
    alpha_min, alpha_max:
        Trim α to a stable window — near 0 and 1 dα/dt and ln(1 − α)
        diverge and dominate the fit.
    n_min, n_max, n_steps:
        Linear grid of candidate n values.
    """
    if len(runs) < 1:
        raise ValueError("Need at least one run")
    if not 0.0 < alpha_min < alpha_max < 1.0:
        raise ValueError(f"Bad α window: [{alpha_min}, {alpha_max}]")
    if n_steps < 2:
        raise ValueError("n_steps must be >= 2")

    inv_T_parts: list[np.ndarray] = []
    ln_lhs_base_parts: list[np.ndarray] = []  # ln(dα/dT) + ln(β) per point
    ln_one_minus_alpha_parts: list[np.ndarray] = []

    for run in runs:
        a = run.alpha
        mask = (a >= alpha_min) & (a <= alpha_max) & (run.dalpha_dt > 0)
        if mask.sum() < 2:
            continue
        T = run.temperature[mask]
        with np.errstate(divide="ignore", invalid="ignore"):
            # ln(dα/dT) + ln(β) = ln(dα/dt), which removes the run-dependent
            # ln(A/β) intercept and lets us pool across heating rates safely.
            ln_lhs_base_parts.append(np.log(run.dalpha_dt[mask]))
            ln_one_minus_alpha_parts.append(np.log(1.0 - a[mask]))
        inv_T_parts.append(1.0 / T)

    if not inv_T_parts:
        raise ValueError("No usable points in the requested α window")

    inv_T = np.concatenate(inv_T_parts)
    ln_lhs_base = np.concatenate(ln_lhs_base_parts)
    ln_one_minus_alpha = np.concatenate(ln_one_minus_alpha_parts)

    finite = np.isfinite(inv_T) & np.isfinite(ln_lhs_base) & np.isfinite(ln_one_minus_alpha)
    inv_T = inv_T[finite]
    ln_lhs_base = ln_lhs_base[finite]
    ln_one_minus_alpha = ln_one_minus_alpha[finite]
    if inv_T.size < 2:
        raise ValueError("Too few finite points to fit")

    n_grid = np.linspace(n_min, n_max, n_steps)
    r2 = np.empty_like(n_grid)
    Ea = np.empty_like(n_grid)
    for i, n in enumerate(n_grid):
        y = ln_lhs_base - n * ln_one_minus_alpha
        slope, _, r = linear_regression(inv_T, y)
        r2[i] = r
        Ea[i] = -slope * R_GAS / 1000.0  # → kJ/mol

    best_idx = int(np.nanargmax(r2))
    return ReactionOrderResult(
        n_grid=n_grid,
        r_squared=r2,
        Ea_kJ_per_mol=Ea,
        n_best=float(n_grid[best_idx]),
        Ea_best_kJ_per_mol=float(Ea[best_idx]),
        r_squared_best=float(r2[best_idx]),
    )


__all__ = ["ReactionOrderResult", "reaction_order"]

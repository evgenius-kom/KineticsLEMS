"""Shared types and helpers for isoconversional methods."""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class IsoconversionalResult:
    """Activation energy as a function of conversion."""

    method: str
    alpha: np.ndarray  # conversions used (dimensionless)
    Ea_J_per_mol: np.ndarray  # E_a(α) in J/mol
    intercept: np.ndarray  # raw intercept of the linear fit (interpretation depends on method)
    r_squared: np.ndarray  # per-α R² of the linear fit

    @property
    def Ea_kJ_per_mol(self) -> np.ndarray:
        return self.Ea_J_per_mol / 1000.0


@dataclass(frozen=True)
class KissingerResult:
    """Single activation energy from peak temperatures."""

    Ea_J_per_mol: float
    pre_exponential: float  # A in 1/s, model-free crude estimate (assumes 1st-order, n=1)
    r_squared: float
    Tp_K: np.ndarray  # peak temperatures used
    rates_K_per_min: np.ndarray

    @property
    def Ea_kJ_per_mol(self) -> float:
        return self.Ea_J_per_mol / 1000.0


def linear_regression(x: np.ndarray, y: np.ndarray) -> tuple[float, float, float]:
    """Return (slope, intercept, R²) of a least-squares fit y ≈ slope·x + intercept."""
    if x.size != y.size or x.size < 2:
        raise ValueError("Need at least 2 paired points for linear regression")
    slope, intercept = np.polyfit(x, y, 1)
    y_pred = slope * x + intercept
    ss_res = float(np.sum((y - y_pred) ** 2))
    ss_tot = float(np.sum((y - np.mean(y)) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")
    return float(slope), float(intercept), r2


__all__ = ["IsoconversionalResult", "KissingerResult", "linear_regression"]

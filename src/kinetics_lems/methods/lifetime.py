"""Predictive isothermal kinetics — α(t) and time-to-α at fixed temperature.

Given the kinetic triplet (E_a(α), A, f(α)) recovered by the
isoconversional + master-plot + pre-exponential pipeline, the Arrhenius
rate equation

    dα/dt  =  A · f(α) · exp(-E_a(α) / (R · T))                  (1)

can be integrated at a user-specified storage / operating temperature T_iso
to predict the time evolution α(t).

This is the headline product for industrial kinetics consultancy:
"how long until 5 % of the material has reacted at 25 °C", or
"what is the safe storage temperature for 1-year shelf life".

Two variants:

* :func:`predict_alpha_of_t` — integrates (1) forward in α from α0 to
  α_target, returning the t(α) curve directly (cheaper, deterministic).
* :func:`time_to_conversion` — convenience wrapper returning a single
  scalar t for a given α_target.

Reference: Vyazovkin (2000), *Thermochim. Acta* 355, 155 — model-based
prediction; ICTAC 2011 §7 — prediction recommendations.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..constants import R_GAS
from .master_plot import MASTER_MODELS, ReactionModel


@dataclass(frozen=True)
class LifetimePrediction:
    """Predicted α(t) at one isothermal temperature."""

    T_K: float
    alpha: np.ndarray
    time_sec: np.ndarray
    model_name: str
    A_per_sec: float

    def time_to_alpha(self, alpha_target: float) -> float:
        """Time (s) needed to reach ``alpha_target``. Linear interpolation."""
        if alpha_target < self.alpha[0] or alpha_target > self.alpha[-1]:
            raise ValueError(
                f"alpha_target={alpha_target} outside range "
                f"[{self.alpha[0]:.3f}, {self.alpha[-1]:.3f}]"
            )
        return float(np.interp(alpha_target, self.alpha, self.time_sec))


def predict_alpha_of_t(
    *,
    T_K: float,
    Ea_J_per_mol: np.ndarray,
    alpha_grid: np.ndarray,
    A_per_sec: float,
    model: str | ReactionModel = "F1",
    alpha_start: float = 0.0,
    eps_f: float = 1e-12,
) -> LifetimePrediction:
    """Integrate dα/dt = A·f(α)·exp(-E_a(α)/RT) at constant T.

    The integral can be rearranged as

        t(α) = ∫_{α_start}^{α}  dα' / [A · f(α') · exp(-E_a(α')/RT)]   (2)

    so a single trapezoid sweep over the supplied α grid gives the
    full α → t mapping with no ODE solver needed.

    Parameters
    ----------
    T_K:
        Isothermal temperature (K).
    Ea_J_per_mol:
        Activation energy as a function of α, sampled at ``alpha_grid``.
    alpha_grid:
        Conversion grid; must be strictly increasing on (0, 1).
    A_per_sec:
        Pre-exponential factor in 1/s (typically from
        :func:`~kinetics_lems.methods.preexponential.compute_A`).
    model:
        Reaction model name (``"F1"`` etc.) or a :class:`ReactionModel`.
    alpha_start:
        α at t = 0. Must be inside or below ``alpha_grid``.
    eps_f:
        Floor for f(α) to avoid division by zero at α → 1.
    """
    if T_K <= 0:
        raise ValueError(f"T_K must be positive, got {T_K}")
    if A_per_sec <= 0:
        raise ValueError(f"A_per_sec must be positive, got {A_per_sec}")
    if Ea_J_per_mol.shape != alpha_grid.shape:
        raise ValueError("Ea_J_per_mol must have the same shape as alpha_grid")
    if not np.all(np.diff(alpha_grid) > 0):
        raise ValueError("alpha_grid must be strictly increasing")
    if alpha_start >= alpha_grid[-1]:
        raise ValueError(
            f"alpha_start={alpha_start} must be below the end of alpha_grid"
        )

    if isinstance(model, str):
        if model not in MASTER_MODELS:
            raise ValueError(f"Unknown model '{model}'. Allowed: {sorted(MASTER_MODELS)}")
        rm = MASTER_MODELS[model]
        model_name = model
    else:
        rm = model
        model_name = model.name

    # Restrict to α >= alpha_start.
    mask = alpha_grid >= alpha_start
    if not mask.any():
        raise ValueError("alpha_start exceeds the supplied alpha_grid")
    a = alpha_grid[mask]
    Ea = Ea_J_per_mol[mask]

    f_alpha = np.maximum(rm.f(a), eps_f)
    # Per-α reaction rate: dα/dt = A · f(α) · exp(-E/RT)
    arrhenius = np.exp(-Ea / (R_GAS * T_K))
    rate = A_per_sec * f_alpha * arrhenius
    rate = np.where(rate > 0, rate, np.nan)

    # t(α) = ∫ dα'/rate(α'). Trapezoid cumulative integral.
    integrand = 1.0 / rate
    t_cum = np.zeros_like(a)
    t_cum[1:] = np.cumsum(0.5 * (integrand[1:] + integrand[:-1]) * np.diff(a))

    return LifetimePrediction(
        T_K=float(T_K),
        alpha=a,
        time_sec=t_cum,
        model_name=model_name,
        A_per_sec=float(A_per_sec),
    )


def time_to_conversion(
    *,
    alpha_target: float,
    T_K: float,
    Ea_J_per_mol: np.ndarray,
    alpha_grid: np.ndarray,
    A_per_sec: float,
    model: str | ReactionModel = "F1",
    alpha_start: float = 0.0,
) -> float:
    """Convenience wrapper: time (s) to reach ``alpha_target`` at ``T_K``.

    Equivalent to ``predict_alpha_of_t(...).time_to_alpha(alpha_target)``.
    """
    prediction = predict_alpha_of_t(
        T_K=T_K,
        Ea_J_per_mol=Ea_J_per_mol,
        alpha_grid=alpha_grid,
        A_per_sec=A_per_sec,
        model=model,
        alpha_start=alpha_start,
    )
    return prediction.time_to_alpha(alpha_target)


@dataclass(frozen=True)
class LifetimeSummary:
    """α(t) predictions at several isothermal temperatures."""

    predictions: list[LifetimePrediction]
    alpha_targets: tuple[float, ...]
    """The α values for which time-to-α was tabulated."""

    times_at_targets: np.ndarray
    """Shape (n_T, n_targets); seconds to reach each target at each T."""

    def temperatures_K(self) -> np.ndarray:
        return np.array([p.T_K for p in self.predictions])


def predict_at_temperatures(
    *,
    T_K_list: list[float],
    Ea_J_per_mol: np.ndarray,
    alpha_grid: np.ndarray,
    A_per_sec: float,
    model: str | ReactionModel = "F1",
    alpha_targets: tuple[float, ...] = (0.05, 0.10, 0.50, 0.90),
    alpha_start: float = 0.0,
) -> LifetimeSummary:
    """Run :func:`predict_alpha_of_t` at several T and tabulate time-to-α."""
    predictions: list[LifetimePrediction] = []
    times = np.full((len(T_K_list), len(alpha_targets)), np.nan)
    for i, T in enumerate(T_K_list):
        pred = predict_alpha_of_t(
            T_K=T,
            Ea_J_per_mol=Ea_J_per_mol,
            alpha_grid=alpha_grid,
            A_per_sec=A_per_sec,
            model=model,
            alpha_start=alpha_start,
        )
        predictions.append(pred)
        for j, target in enumerate(alpha_targets):
            if pred.alpha[0] <= target <= pred.alpha[-1]:
                times[i, j] = pred.time_to_alpha(target)
    return LifetimeSummary(
        predictions=predictions,
        alpha_targets=tuple(alpha_targets),
        times_at_targets=times,
    )


__all__ = [
    "LifetimePrediction",
    "LifetimeSummary",
    "predict_alpha_of_t",
    "predict_at_temperatures",
    "time_to_conversion",
]

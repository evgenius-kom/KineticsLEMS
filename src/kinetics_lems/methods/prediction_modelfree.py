"""Model-free isothermal prediction (Vyazovkin 2000) — no f(α) needed.

Standard prediction in :mod:`kinetics_lems.methods.lifetime` requires a
kinetic triplet (E(α), A, f(α)) — picking the wrong f(α) can dominate the
prediction error. Vyazovkin (2000) showed that under the isoconversional
assumption, t(α; T_iso) can be computed from E(α) **alone**:

    t(α; T_iso) = (∫₀^{T_ref} exp(−E_α / R T(t')) dt')  ·  exp(E_α / R T_iso)
                                                                       (1)

where the reference integral on the right is taken from one of the
*experimental* heating-rate programs at the same α. The ratio of
temperature integrals folds the unknown f(α) out: as long as α is the
same, ``A · f(α)`` cancels.

Practical form (linear heating, constant β, α reached at T_α):

    t(α; T_iso) = (I(E_α, T_α) / β_exp) · exp(E_α / R T_iso)            (2)

with the Senum-Yang p(x) approximation for the temperature integral.

Reference:
    Vyazovkin (2000), *Thermochim. Acta* 355, 155–163, eq. (16).
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..constants import R_GAS, SEC_PER_MIN


@dataclass(frozen=True)
class ModelFreePrediction:
    """t(α) and α(t) at one isothermal temperature, from E(α) only."""

    T_K: float
    alpha: np.ndarray
    time_sec: np.ndarray
    """t for each requested α; monotone non-decreasing."""

    method: str = "vyazovkin_modelfree"

    def time_to_alpha(self, target: float) -> float:
        if target < self.alpha[0] or target > self.alpha[-1]:
            raise ValueError(
                f"alpha_target={target} outside range "
                f"[{self.alpha[0]:.3f}, {self.alpha[-1]:.3f}]"
            )
        return float(np.interp(target, self.alpha, self.time_sec))


def _p_senum_yang(x: np.ndarray) -> np.ndarray:
    """Senum-Yang rational approximation of p(x). Same as vyazovkin.py."""
    num = x * x + 10.0 * x + 18.0
    den = x * x * x + 12.0 * x * x + 36.0 * x + 24.0
    return np.exp(-x) / x * num / den


def predict_isothermal_modelfree(
    *,
    T_K: float,
    Ea_J_per_mol: np.ndarray,
    T_alpha_K: np.ndarray,
    beta_K_per_min: float,
    alpha_grid: np.ndarray,
) -> ModelFreePrediction:
    """Compute t(α) at constant ``T_K`` from E(α) and a reference β.

    Parameters
    ----------
    T_K :
        Target isothermal temperature.
    Ea_J_per_mol :
        E_α sampled on ``alpha_grid``.
    T_alpha_K :
        Temperatures at which each α was reached during the reference
        experiment (same shape as ``alpha_grid``). Comes from
        :func:`kinetics_lems.conversion.temperature_at_conversion` on one
        of the experimental runs.
    beta_K_per_min :
        Heating rate of the reference experiment, K/min.
    alpha_grid :
        Conversion grid; must match ``Ea`` and ``T_alpha_K``.

    Returns
    -------
    ModelFreePrediction
        t(α) at the requested isothermal temperature.
    """
    if T_K <= 0:
        raise ValueError(f"T_K must be positive, got {T_K}")
    if beta_K_per_min <= 0:
        raise ValueError(f"beta_K_per_min must be positive, got {beta_K_per_min}")
    if not (Ea_J_per_mol.shape == T_alpha_K.shape == alpha_grid.shape):
        raise ValueError("Ea, T_alpha and alpha_grid must share shape")

    beta_K_per_sec = beta_K_per_min / SEC_PER_MIN
    # I(E, T) ≈ (E / R) · p(x), x = E / (R T)
    x_ref = Ea_J_per_mol / (R_GAS * T_alpha_K)
    I_ref = (Ea_J_per_mol / R_GAS) * _p_senum_yang(x_ref)
    # t(α) = I_ref / β · exp(E / R T_iso)
    arrhenius = np.exp(Ea_J_per_mol / (R_GAS * T_K))
    t = I_ref / beta_K_per_sec * arrhenius
    # Numerical artefacts: t must be ≥ 0 and monotone non-decreasing on α.
    t = np.maximum(t, 0.0)
    t = np.maximum.accumulate(t)
    return ModelFreePrediction(T_K=float(T_K), alpha=alpha_grid.copy(), time_sec=t)


def predict_arbitrary_program_modelfree(
    *,
    time_s: np.ndarray,
    T_K_program: np.ndarray,
    Ea_J_per_mol: np.ndarray,
    T_alpha_K: np.ndarray,
    beta_K_per_min: float,
    alpha_grid: np.ndarray,
) -> ModelFreePrediction:
    """Predict α(t) under arbitrary T(t) without assuming f(α).

    Uses the same isoconversional identity as the isothermal version:
    for any program T(t), the time the system spends "at α" is

        τ(α; program) = ∫_{0}^{t} exp(−E_α / R T(t')) dt'                (3)

    and the prediction is the t at which τ equals the *reference*
    integral I_ref(α) (eq. 2 in Vyazovkin 2000).

    Parameters
    ----------
    time_s, T_K_program :
        The target T(t) program, both 1-D arrays of equal length.
    Ea_J_per_mol, T_alpha_K, beta_K_per_min, alpha_grid :
        Same meaning as in :func:`predict_isothermal_modelfree`.
    """
    if time_s.shape != T_K_program.shape:
        raise ValueError("time_s and T_K_program must share shape")
    if not (Ea_J_per_mol.shape == T_alpha_K.shape == alpha_grid.shape):
        raise ValueError("Ea, T_alpha and alpha_grid must share shape")
    if not np.all(np.diff(time_s) > 0):
        raise ValueError("time_s must be strictly increasing")

    beta_K_per_sec = beta_K_per_min / SEC_PER_MIN
    x_ref = Ea_J_per_mol / (R_GAS * T_alpha_K)
    I_ref = (Ea_J_per_mol / R_GAS) * _p_senum_yang(x_ref) / beta_K_per_sec

    t_predicted = np.full_like(alpha_grid, np.nan, dtype=float)
    for i, E in enumerate(Ea_J_per_mol):
        if not np.isfinite(E):
            continue
        integrand = np.exp(-E / (R_GAS * T_K_program))
        tau = np.zeros_like(time_s)
        tau[1:] = np.cumsum(0.5 * (integrand[1:] + integrand[:-1]) * np.diff(time_s))
        target = float(I_ref[i])
        if tau[-1] < target:
            t_predicted[i] = np.nan  # α never reached within program
            continue
        t_predicted[i] = float(np.interp(target, tau, time_s))

    # Trim NaN endpoint and enforce monotonicity.
    finite = np.isfinite(t_predicted)
    if finite.any():
        max_t_so_far = -np.inf
        for i in range(t_predicted.size):
            if np.isfinite(t_predicted[i]):
                if t_predicted[i] < max_t_so_far:
                    t_predicted[i] = max_t_so_far
                else:
                    max_t_so_far = t_predicted[i]

    return ModelFreePrediction(
        T_K=float("nan"),  # not isothermal
        alpha=alpha_grid.copy(),
        time_sec=t_predicted,
    )


__all__ = [
    "ModelFreePrediction",
    "predict_arbitrary_program_modelfree",
    "predict_isothermal_modelfree",
]

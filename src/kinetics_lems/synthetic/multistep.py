"""Synthetic datasets for multi-step and distributed-reactivity validation.

Used as ground-truth test beds for future DAEM and multi-step fitters
(see :mod:`kinetics_lems.fitting`). Each generator returns a
:class:`CaseData` shaped exactly like a real experiment, so the same
loader/pipeline can consume them.

Generators:
    * :func:`generate_two_parallel_case`  — two independent F1 channels
      with their own (E_i, A_i) and a weighting w_i.
    * :func:`generate_daem_gaussian_case` — Gaussian DAEM (continuum of E
      centered at Ē with stddev σ, shared A).
    * :func:`generate_arbitrary_program_case` — a single F1 reaction
      under a non-linear T(t) program (modulated DSC analogue), written
      out as 3-column files.

All curves are stored as ``dα/dT`` so the existing 2-column file format
is preserved for parallel/DAEM and the 3-column format for the arbitrary
program (see ``io/wave_reader.py``).
"""
from __future__ import annotations

from collections.abc import Sequence

import numpy as np

from ..constants import R_GAS, SEC_PER_MIN
from ..models import CaseData, CaseParams, ExperimentType, Method, Wave


def _f_first_order(a: float) -> float:
    return max(1.0 - a, 0.0)


def _integrate_f1_dynamic(
    T_grid: np.ndarray,
    beta_K_per_sec: float,
    E_J: float,
    A: float,
) -> np.ndarray:
    """Integrate dα/dT = (A/β)·(1−α)·exp(−E/RT) by RK4 on T."""
    alpha = np.zeros_like(T_grid)
    for k in range(1, T_grid.size):
        T0 = T_grid[k - 1]
        T1 = T_grid[k]
        dT = T1 - T0
        a = alpha[k - 1]
        rhs = lambda aa, TT: (A / beta_K_per_sec) * max(1.0 - aa, 0.0) * np.exp(  # noqa: E731
            -E_J / (R_GAS * TT)
        )
        k1 = rhs(a, T0)
        k2 = rhs(a + 0.5 * dT * k1, 0.5 * (T0 + T1))
        k3 = rhs(a + 0.5 * dT * k2, 0.5 * (T0 + T1))
        k4 = rhs(a + dT * k3, T1)
        alpha[k] = min(a + dT * (k1 + 2 * k2 + 2 * k3 + k4) / 6.0, 1.0)
    return alpha


def generate_two_parallel_case(
    rates_K_per_min: Sequence[float],
    *,
    E1_J_per_mol: float = 100_000.0,
    A1_per_sec: float = 1.0e9,
    weight1: float = 0.5,
    E2_J_per_mol: float = 180_000.0,
    A2_per_sec: float = 1.0e12,
    T_start_K: float = 400.0,
    T_stop_K: float = 800.0,
    n_points: int = 1500,
    material: str = "synthetic-2parallel",
) -> CaseData:
    """Two independent first-order parallel channels.

    α_total(T) = w1·α1(T) + (1−w1)·α2(T) where each α_i comes from its
    own (E_i, A_i, F1) integration. ``dα/dT`` is built from the
    sum of the two channel derivatives.
    """
    if not 0.0 < weight1 < 1.0:
        raise ValueError(f"weight1 must lie strictly in (0, 1), got {weight1}")
    weight2 = 1.0 - weight1
    T_grid = np.linspace(T_start_K, T_stop_K, n_points)

    waves: dict[float, Wave] = {}
    for beta_min in rates_K_per_min:
        beta_sec = beta_min / SEC_PER_MIN
        a1 = _integrate_f1_dynamic(T_grid, beta_sec, E1_J_per_mol, A1_per_sec)
        a2 = _integrate_f1_dynamic(T_grid, beta_sec, E2_J_per_mol, A2_per_sec)
        a_total = weight1 * a1 + weight2 * a2
        dadT = np.gradient(a_total, T_grid)
        waves[float(beta_min)] = Wave(x=T_grid.copy(), y=dadT)

    file_map = {f"rate_{i:02d}.txt": float(b) for i, b in enumerate(sorted(rates_K_per_min))}
    params = CaseParams(
        material=material,
        experiment_type=ExperimentType.HEATING,
        method=Method.DSC,
        file_to_condition=file_map,
    )
    return CaseData(params=params, waves=waves)


def generate_daem_gaussian_case(
    rates_K_per_min: Sequence[float],
    *,
    E_mean_J_per_mol: float = 150_000.0,
    E_sigma_J_per_mol: float = 10_000.0,
    A_per_sec: float = 1.0e11,
    n_E_quadrature: int = 41,
    T_start_K: float = 400.0,
    T_stop_K: float = 900.0,
    n_points: int = 1500,
    material: str = "synthetic-daem-gaussian",
) -> CaseData:
    """Gaussian-distributed activation energy, F1 channels, common A.

    α_total(T) = ∫ p(E) · α(T; E, A) dE with p(E) ~ N(Ē, σ).
    Discretised by Gauss-Hermite-like equidistant sampling on
    [Ē−4σ, Ē+4σ].
    """
    Es = np.linspace(
        E_mean_J_per_mol - 4 * E_sigma_J_per_mol,
        E_mean_J_per_mol + 4 * E_sigma_J_per_mol,
        n_E_quadrature,
    )
    weights = np.exp(-((Es - E_mean_J_per_mol) ** 2) / (2 * E_sigma_J_per_mol**2))
    weights = weights / weights.sum()
    T_grid = np.linspace(T_start_K, T_stop_K, n_points)

    waves: dict[float, Wave] = {}
    for beta_min in rates_K_per_min:
        beta_sec = beta_min / SEC_PER_MIN
        a_total = np.zeros_like(T_grid)
        for E_i, w_i in zip(Es, weights, strict=True):
            a_i = _integrate_f1_dynamic(T_grid, beta_sec, float(E_i), A_per_sec)
            a_total += w_i * a_i
        dadT = np.gradient(a_total, T_grid)
        waves[float(beta_min)] = Wave(x=T_grid.copy(), y=dadT)

    file_map = {f"rate_{i:02d}.txt": float(b) for i, b in enumerate(sorted(rates_K_per_min))}
    params = CaseParams(
        material=material,
        experiment_type=ExperimentType.HEATING,
        method=Method.DSC,
        file_to_condition=file_map,
    )
    return CaseData(params=params, waves=waves)


def generate_arbitrary_program_case(
    *,
    E_J_per_mol: float = 120_000.0,
    A_per_sec: float = 1.0e10,
    nominal_rate_K_per_min: float = 10.0,
    program: str = "modulated",
    t_max_s: float = 3600.0,
    n_points: int = 4001,
    material: str = "synthetic-modulated",
) -> CaseData:
    """One run under a non-linear T(t) program.

    ``program="modulated"`` superimposes a 2 K, 60 s sinusoidal modulation on
    a 5 K/min underlying ramp — a clean stand-in for modulated DSC.
    ``program="step"`` uses two isothermal holds connected by a ramp.

    The returned :class:`Wave` carries ``t_seconds``, so reading it through
    the standard pipeline exercises the 3-column branch in
    :mod:`kinetics_lems.conversion`.
    """
    t = np.linspace(0.0, t_max_s, n_points)
    if program == "modulated":
        T = 400.0 + (nominal_rate_K_per_min / SEC_PER_MIN) * t + 2.0 * np.sin(
            2 * np.pi * t / 60.0
        )
    elif program == "step":
        T = np.where(
            t < t_max_s / 3,
            500.0,
            np.where(t < 2 * t_max_s / 3, 500.0 + 0.5 * (t - t_max_s / 3), 600.0),
        )
    else:
        raise ValueError(f"Unknown program {program!r}; expected 'modulated' or 'step'")

    # Forward-integrate dα/dt = A·(1−α)·exp(-E/RT(t)) on the recorded t.
    alpha = np.zeros_like(t)
    for k in range(1, t.size):
        a = alpha[k - 1]
        T_mid = 0.5 * (T[k - 1] + T[k])
        rate = A_per_sec * max(1.0 - a, 0.0) * np.exp(-E_J_per_mol / (R_GAS * T_mid))
        alpha[k] = min(a + rate * (t[k] - t[k - 1]), 1.0)
    dadt = np.gradient(alpha, t)

    wave = Wave(x=T.copy(), y=dadt, t_seconds=t.copy())
    file_map = {"run_arbitrary.txt": nominal_rate_K_per_min}
    params = CaseParams(
        material=material,
        experiment_type=ExperimentType.HEATING,
        method=Method.DSC,
        file_to_condition=file_map,
    )
    return CaseData(params=params, waves={nominal_rate_K_per_min: wave})


__all__ = [
    "generate_arbitrary_program_case",
    "generate_daem_gaussian_case",
    "generate_two_parallel_case",
]

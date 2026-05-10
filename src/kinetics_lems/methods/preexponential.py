"""Pre-exponential factor A from per-α activation energy + a chosen f(α).

Closes the kinetic triplet (E, A, f) once Vyazovkin (or any isoconversional
method) gives E_a(α) and the master plot identifies f(α). For each (α, run)
pair compute

    A(α, run) = (dα/dt)_α,run / [ f(α) · exp(-E_a(α) / (R · T_α,run)) ]

then aggregate across α and runs (median + MAD for robustness).

Reference: Vyazovkin et al. (2011) ICTAC recommendations §6 — "Determining
the pre-exponential factor".
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..constants import R_GAS
from ..conversion import ConversionRun, dalpha_dt_at_conversion, temperature_at_conversion
from .master_plot import MASTER_MODELS, ReactionModel


@dataclass(frozen=True)
class PreexponentialResult:
    """Per-α pre-exponential A under a chosen reaction model."""

    model_name: str
    alpha: np.ndarray
    A_per_sec_per_alpha: np.ndarray   # median A across runs at each α
    A_per_sec_median: float           # global median (over α and runs)
    A_per_sec_mad: float              # median absolute deviation (robust spread)
    log10_A_median: float
    log10_A_mad: float


def compute_A(
    runs: list[ConversionRun],
    alphas: np.ndarray,
    Ea_J_per_mol: np.ndarray,
    *,
    model: str | ReactionModel = "F1",
    eps_f: float = 1e-12,
) -> PreexponentialResult:
    """Compute A(α) under ``model``.

    Parameters
    ----------
    runs:
        ConversionRun objects (must cover the requested α grid).
    alphas:
        α values at which E_a was estimated.
    Ea_J_per_mol:
        Activation energy in J/mol per α (e.g., from Vyazovkin). Must be
        the same length as ``alphas``.
    model:
        Either a model name (e.g. ``"F1"``) or a custom :class:`ReactionModel`.
    eps_f:
        Lower clip for f(α) — prevents division-by-zero at α → 1 for
        first-order-style models.
    """
    if Ea_J_per_mol.shape != alphas.shape:
        raise ValueError("Ea must have the same shape as alphas")

    if isinstance(model, str):
        if model not in MASTER_MODELS:
            raise ValueError(f"Unknown model '{model}'. Allowed: {sorted(MASTER_MODELS)}")
        rm = MASTER_MODELS[model]
        model_name = model
    else:
        rm = model
        model_name = model.name

    f_alpha = np.maximum(rm.f(alphas), eps_f)

    # Build A(α, run) matrix.
    A_table = np.empty((len(runs), alphas.size), dtype=float)
    for i, run in enumerate(runs):
        T = temperature_at_conversion(run, alphas)
        rate = dalpha_dt_at_conversion(run, alphas)
        with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
            arrhenius = np.exp(-Ea_J_per_mol / (R_GAS * T))
            A_table[i, :] = rate / (f_alpha * arrhenius)
    A_table = np.where(np.isfinite(A_table), A_table, np.nan)

    # Per-α: median across runs.
    A_per_alpha = np.nanmedian(A_table, axis=0)

    # Robust scalar summary across α and runs.
    flat = A_table[np.isfinite(A_table)]
    flat = flat[flat > 0]  # log space requires positive values
    if flat.size == 0:
        median = float("nan")
        mad = float("nan")
        log_median = float("nan")
        log_mad = float("nan")
    else:
        log_vals = np.log10(flat)
        log_median = float(np.median(log_vals))
        log_mad = float(np.median(np.abs(log_vals - log_median)))
        median = float(10.0**log_median)
        mad = float(median * (10.0**log_mad - 1.0))  # convert log-MAD to linear

    return PreexponentialResult(
        model_name=model_name,
        alpha=alphas,
        A_per_sec_per_alpha=A_per_alpha,
        A_per_sec_median=median,
        A_per_sec_mad=mad,
        log10_A_median=log_median,
        log10_A_mad=log_mad,
    )


__all__ = ["PreexponentialResult", "compute_A"]

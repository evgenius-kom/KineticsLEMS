"""Coats–Redfern (1964) — single-rate model-fitting linearization.

For each *individual* heating rate β and each candidate reaction model
with integral form g(α), the relation

    ln[g(α) / T²]  ≈  ln(AR / (βE))  −  E / (R·T)             (1)

is linear in 1/T with slope −E/R. Fit (1) for every model and every
run; rank models by R² (averaged across runs).

Unlike isoconversional methods, Coats–Redfern *assumes a model*. Useful
as an independent cross-check on the master-plot Z(α) ranking — if both
agree, model identification is solid; if they disagree, the kinetics
may be multi-step or outside the canonical 12-model set.

Reference: Coats, A. W.; Redfern, J. P. (1964). *Nature* 201, 68.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..constants import R_GAS, SEC_PER_MIN
from ..conversion import ConversionRun
from .common import linear_regression
from .master_plot import MASTER_MODELS


@dataclass(frozen=True)
class CoatsRedfernRunFit:
    """Coats–Redfern fit for one (model, run) pair."""

    model: str
    rate_K_per_min: float
    Ea_kJ_per_mol: float
    A_per_sec: float          # extracted from the intercept assuming linear heating
    r_squared: float
    n_points: int
    aic: float = float("nan")
    """Akaike information criterion of the linear fit (k=2 parameters)."""
    bic: float = float("nan")
    """Bayesian information criterion of the linear fit (k=2 parameters)."""


@dataclass(frozen=True)
class CoatsRedfernModelSummary:
    """Aggregated Coats–Redfern fit across runs for one model."""

    model: str
    Ea_kJ_per_mol_mean: float
    Ea_kJ_per_mol_std: float
    log10_A_mean: float
    log10_A_std: float
    r_squared_mean: float
    """Mean R² across runs — higher = better fit of this model."""

    aic_mean: float
    """Mean AIC across runs — lower = better fit *per fixed model complexity*."""
    bic_mean: float

    n_runs: int


@dataclass(frozen=True)
class CoatsRedfernResult:
    fits: list[CoatsRedfernRunFit]
    summaries: list[CoatsRedfernModelSummary]
    """Sorted descending by mean R² — first entry is the best-fitting model."""

    @property
    def best_model(self) -> str:
        return self.summaries[0].model

    def by_model(self, name: str) -> CoatsRedfernModelSummary:
        for s in self.summaries:
            if s.model == name:
                return s
        raise KeyError(f"Model '{name}' not in result")


def coats_redfern(
    runs: list[ConversionRun],
    *,
    alpha_min: float = 0.1,
    alpha_max: float = 0.9,
    models: list[str] | None = None,
) -> CoatsRedfernResult:
    """Coats–Redfern fit for every model and every run.

    Parameters
    ----------
    runs:
        ConversionRun list.
    alpha_min, alpha_max:
        α window to fit. Excluding the tails removes the
        ``g(α) → 0 or ∞`` problem near the endpoints.
    models:
        Subset of MASTER_MODELS to evaluate. None → all 12.
    """
    if len(runs) < 1:
        raise ValueError("Need at least one run")
    if not 0.0 < alpha_min < alpha_max < 1.0:
        raise ValueError(f"Bad α window: [{alpha_min}, {alpha_max}]")

    model_names = list(models) if models is not None else list(MASTER_MODELS)
    for m in model_names:
        if m not in MASTER_MODELS:
            raise ValueError(f"Unknown model '{m}'")

    fits: list[CoatsRedfernRunFit] = []
    by_model: dict[str, list[CoatsRedfernRunFit]] = {m: [] for m in model_names}

    for run in runs:
        a = run.alpha
        mask = (a >= alpha_min) & (a <= alpha_max)
        if mask.sum() < 3:
            continue
        T = run.temperature[mask]
        alpha = a[mask]
        beta_K_per_sec = run.rate_K_per_min / SEC_PER_MIN
        inv_T = 1.0 / T

        for name in model_names:
            model = MASTER_MODELS[name]
            with np.errstate(divide="ignore", invalid="ignore"):
                g = model.g(alpha)
                y = np.log(g / (T * T))
            finite = np.isfinite(y) & (g > 0)
            if finite.sum() < 3:
                continue
            slope, intercept, r2 = linear_regression(inv_T[finite], y[finite])
            Ea_J = -slope * R_GAS
            if Ea_J <= 0:
                continue
            # intercept ≈ ln(A·R / (β·E))  →  A = β·E·exp(intercept) / R
            A = beta_K_per_sec * Ea_J * np.exp(intercept) / R_GAS
            n_pts = int(finite.sum())
            residuals = y[finite] - (slope * inv_T[finite] + intercept)
            aic, bic = _aic_bic_ols(residuals, k=2)
            fit = CoatsRedfernRunFit(
                model=name,
                rate_K_per_min=run.rate_K_per_min,
                Ea_kJ_per_mol=Ea_J / 1000.0,
                A_per_sec=float(A),
                r_squared=r2,
                n_points=n_pts,
                aic=aic,
                bic=bic,
            )
            fits.append(fit)
            by_model[name].append(fit)

    summaries: list[CoatsRedfernModelSummary] = []
    for name, fits_for_model in by_model.items():
        if not fits_for_model:
            continue
        Eas = np.array([f.Ea_kJ_per_mol for f in fits_for_model])
        As = np.array([f.A_per_sec for f in fits_for_model])
        r2s = np.array([f.r_squared for f in fits_for_model])
        aics = np.array([f.aic for f in fits_for_model])
        bics = np.array([f.bic for f in fits_for_model])
        with np.errstate(divide="ignore", invalid="ignore"):
            log_As = np.log10(np.where(As > 0, As, np.nan))
        summaries.append(
            CoatsRedfernModelSummary(
                model=name,
                Ea_kJ_per_mol_mean=float(np.mean(Eas)),
                Ea_kJ_per_mol_std=float(np.std(Eas, ddof=0)),
                log10_A_mean=float(np.nanmean(log_As)),
                log10_A_std=float(np.nanstd(log_As, ddof=0)),
                r_squared_mean=float(np.mean(r2s)),
                aic_mean=float(np.nanmean(aics)),
                bic_mean=float(np.nanmean(bics)),
                n_runs=len(fits_for_model),
            )
        )

    if not summaries:
        raise ValueError("No successful Coats–Redfern fit on any (model, run) pair")

    summaries.sort(key=lambda s: s.r_squared_mean, reverse=True)
    return CoatsRedfernResult(fits=fits, summaries=summaries)


def _aic_bic_ols(residuals: np.ndarray, *, k: int) -> tuple[float, float]:
    """AIC/BIC for an ordinary-least-squares fit.

    Under Gaussian residuals the maximised log-likelihood reduces to
    ``-n/2 · ln(RSS/n) - const`` (Burnham & Anderson, *Model Selection
    and Multimodel Inference*, 2002, §2.2):

        AIC = n · ln(RSS / n) + 2 · k
        BIC = n · ln(RSS / n) + k · ln(n)

    The additive constant ``n · (1 + ln(2π))`` is the same for every
    model on the same data and is dropped — only AIC *differences* are
    meaningful anyway.
    """
    n = residuals.size
    if n <= k:
        return float("nan"), float("nan")
    rss = float(np.sum(residuals * residuals))
    if rss <= 0:
        # Perfect fit (synthetic / degenerate) — AIC → −∞, return finite floor.
        rss = 1e-300
    log_rss_over_n = float(np.log(rss / n))
    aic = n * log_rss_over_n + 2.0 * k
    bic = n * log_rss_over_n + k * float(np.log(n))
    return aic, bic


__all__ = [
    "CoatsRedfernModelSummary",
    "CoatsRedfernResult",
    "CoatsRedfernRunFit",
    "coats_redfern",
]

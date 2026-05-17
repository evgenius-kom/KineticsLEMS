"""Uncertainty estimates for E(α) via leave-one-run-out jackknife.

Isoconversional methods do not have a closed-form covariance matrix
because they pool data across heating-rate runs, not across independent
identical samples. The pragmatic substitute recommended by ICTAC 2011 §3
is *resampling over runs*: rerun the analysis on every possible
(n − 1)-run subset and take the spread of E(α) across subsets as the
uncertainty.

Implementation choice — jackknife (leave-one-out) rather than bootstrap:

* deterministic (no RNG seed needed for reproducibility);
* needs only n recomputations vs ~10³ for bootstrap;
* the bias-corrected jackknife SE is appropriate for moderately
  smooth statistics like E(α), per Efron & Tibshirani (1993) §11.

For n runs and a base estimator T(runs):

    T_(i)      = estimator evaluated on runs with run i removed
    T_(.)      = mean over i of T_(i)
    SE_jack(α) = sqrt[ (n − 1) / n  ·  Σ_i (T_(i)(α) − T_(.)(α))² ]    (1)

This SE has the right order of magnitude for E(α) when n ≥ 3 runs;
with only 2 runs every leave-one-out estimator drops to 1 run, which
isoconversional methods refuse — so jackknife is reported as NaN.
"""
from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass

import numpy as np

from ..conversion import ConversionRun
from .common import IsoconversionalResult

IsoconversionalEstimator = Callable[[list[ConversionRun], np.ndarray], IsoconversionalResult]


@dataclass(frozen=True)
class UncertaintyResult:
    """Jackknife-by-run uncertainty for one isoconversional method."""

    method: str
    alpha: np.ndarray
    Ea_kJ_per_mol_mean: np.ndarray  # mean across leave-one-out subsets, per α
    Ea_kJ_per_mol_se: np.ndarray    # jackknife standard error per α
    n_runs: int
    """Total number of runs used."""

    @property
    def Ea_kJ_per_mol_ci95_low(self) -> np.ndarray:
        """Approximate 95% CI lower bound (mean − 1.96·SE)."""
        return self.Ea_kJ_per_mol_mean - 1.96 * self.Ea_kJ_per_mol_se

    @property
    def Ea_kJ_per_mol_ci95_high(self) -> np.ndarray:
        return self.Ea_kJ_per_mol_mean + 1.96 * self.Ea_kJ_per_mol_se


def jackknife_isoconversional(
    runs: list[ConversionRun],
    alphas: np.ndarray,
    estimator: IsoconversionalEstimator,
    *,
    method_name: str | None = None,
) -> UncertaintyResult:
    """Leave-one-run-out jackknife of an isoconversional estimator.

    Parameters
    ----------
    runs:
        Full list of n ≥ 3 conversion runs. With n = 2, every leave-one-out
        subset has only 1 run, which all isoconversional methods reject —
        in that case SE is reported as NaN (no uncertainty available).
    alphas:
        Conversion grid passed to the estimator.
    estimator:
        Callable ``(runs, alphas) → IsoconversionalResult``. Typically
        a partial of :func:`vyazovkin` or any other E(α) method.
    method_name:
        Display name; defaults to the estimator's __name__.
    """
    n = len(runs)
    name = method_name or getattr(estimator, "__name__", "estimator")

    if n < 3:
        return UncertaintyResult(
            method=name,
            alpha=alphas,
            Ea_kJ_per_mol_mean=np.full_like(alphas, np.nan, dtype=float),
            Ea_kJ_per_mol_se=np.full_like(alphas, np.nan, dtype=float),
            n_runs=n,
        )

    leave_one_out = np.empty((n, alphas.size), dtype=float)
    for i in range(n):
        subset = [r for j, r in enumerate(runs) if j != i]
        result = estimator(subset, alphas)
        leave_one_out[i, :] = result.Ea_kJ_per_mol

    mean = np.nanmean(leave_one_out, axis=0)
    # Jackknife SE per Efron & Tibshirani (1993) eq. 11.5:
    #   SE = sqrt[ (n-1)/n · Σ_i (θ_(i) − θ_(.))² ]
    diffs = leave_one_out - mean[None, :]
    variance = (n - 1) / n * np.nansum(diffs * diffs, axis=0)
    se = np.sqrt(variance)

    return UncertaintyResult(
        method=name,
        alpha=alphas,
        Ea_kJ_per_mol_mean=mean,
        Ea_kJ_per_mol_se=se,
        n_runs=n,
    )


__all__ = ["UncertaintyResult", "jackknife_isoconversional"]

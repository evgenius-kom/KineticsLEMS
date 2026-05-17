"""Diagnostic checks for isoconversional E(α) curves.

ICTAC 2020 §3 notes that

* α < 0.1 and α > 0.9 are typically unreliable: baseline subtraction
  errors dominate the tails;
* variability of E_α(α) greater than ~10–20 % of the average E_α is
  a strong warning that the process is not single-step.

The functions here surface those two flags in a structured way so they
can be emitted by the CLI and embedded in reports.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .common import IsoconversionalResult


@dataclass(frozen=True)
class EndpointReliability:
    """Per-method assessment of how trustworthy α-endpoints are."""

    method: str
    alpha_low: float
    alpha_high: float
    """α range we treat as the *reliable core* (default 0.1–0.9)."""

    finite_fraction_in_core: float
    finite_fraction_in_tail_low: float
    finite_fraction_in_tail_high: float
    """Fraction of finite (non-NaN) E values inside / below / above the core."""

    flatness_in_core: float
    """(max−min)/median of E_a inside the core; ICTAC §3 multi-step diagnostic."""

    warnings: list[str]


def assess_endpoints(
    iso: IsoconversionalResult,
    *,
    alpha_low: float = 0.1,
    alpha_high: float = 0.9,
    flatness_threshold: float = 0.10,
    finite_fraction_threshold: float = 0.5,
) -> EndpointReliability:
    """Score the reliability of an :class:`IsoconversionalResult`.

    Parameters
    ----------
    alpha_low, alpha_high :
        Boundaries of the *reliable core* α range.
    flatness_threshold :
        Above this relative spread of E_a in the core, emit a "multi-step
        suspected" warning. 0.10 is the ICTAC §3 lower bound;
        0.20 is the upper bound (textbook).
    finite_fraction_threshold :
        Tails with fewer than this fraction of finite points raise a
        "tail mostly invalid" warning (signal-to-noise too low to fit).
    """
    alpha = np.asarray(iso.alpha, dtype=float)
    E = np.asarray(iso.Ea_J_per_mol, dtype=float)

    core_mask = (alpha >= alpha_low) & (alpha <= alpha_high)
    low_mask = alpha < alpha_low
    high_mask = alpha > alpha_high

    def _finite_fraction(mask: np.ndarray) -> float:
        if not mask.any():
            return float("nan")
        return float(np.isfinite(E[mask]).sum() / mask.sum())

    core_finite = _finite_fraction(core_mask)
    tail_low_finite = _finite_fraction(low_mask)
    tail_high_finite = _finite_fraction(high_mask)

    core_E = E[core_mask & np.isfinite(E)]
    if core_E.size >= 2:
        median = float(np.median(core_E))
        flatness = float((core_E.max() - core_E.min()) / max(abs(median), 1.0))
    else:
        flatness = float("nan")

    warnings: list[str] = []
    if np.isfinite(core_finite) and core_finite < 1.0:
        warnings.append(
            f"{iso.method}: only {core_finite * 100:.0f}% of core α∈"
            f"[{alpha_low}, {alpha_high}] yielded a finite E"
        )
    if np.isfinite(tail_low_finite) and tail_low_finite < finite_fraction_threshold:
        warnings.append(
            f"{iso.method}: low-α tail α<{alpha_low} mostly invalid "
            f"({tail_low_finite * 100:.0f}% finite); rates near α=0 too small to fit"
        )
    if np.isfinite(tail_high_finite) and tail_high_finite < finite_fraction_threshold:
        warnings.append(
            f"{iso.method}: high-α tail α>{alpha_high} mostly invalid "
            f"({tail_high_finite * 100:.0f}% finite); baseline drift dominates"
        )
    if np.isfinite(flatness) and flatness > flatness_threshold:
        warnings.append(
            f"{iso.method}: E_a varies by {flatness * 100:.0f}% across α∈"
            f"[{alpha_low}, {alpha_high}] (threshold {flatness_threshold * 100:.0f}%); "
            "single-step assumption suspect"
        )

    return EndpointReliability(
        method=iso.method,
        alpha_low=alpha_low,
        alpha_high=alpha_high,
        finite_fraction_in_core=core_finite,
        finite_fraction_in_tail_low=tail_low_finite,
        finite_fraction_in_tail_high=tail_high_finite,
        flatness_in_core=flatness,
        warnings=warnings,
    )


__all__ = ["EndpointReliability", "assess_endpoints"]

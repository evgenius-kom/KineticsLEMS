"""Basic multi-step detection from E(α).

A non-flat E_a(α) curve is the textbook diagnostic of a multi-step
process — when E changes by more than ~10–20 % across α, the underlying
reaction is not a single elementary step (ICTAC 2011 §6, ICTAC 2020 §3).

This module performs *rule-based* segmentation: it splits the α range
into contiguous segments where E is approximately constant, then refits
each segment with a single-step kinetic triplet (E from per-segment
isoconversional analysis, A from the chosen f(α), step contribution
from the cumulative Δα weight).

Out of scope here: full nonlinear least-squares deconvolution of a sum
of N independent kinetic triplets — see the "model-based" follow-up in
[docs/ROADMAP_PROMPT.md](../../docs/ROADMAP_PROMPT.md).
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .common import IsoconversionalResult


@dataclass(frozen=True)
class KineticStep:
    """One segment of an E(α) curve treated as an elementary step."""

    index: int
    alpha_lo: float
    alpha_hi: float
    Ea_kJ_per_mol_median: float
    Ea_kJ_per_mol_mad: float  # robust spread inside the segment
    contribution: float       # fraction of total reaction (= alpha_hi − alpha_lo)


@dataclass(frozen=True)
class MultiStepResult:
    n_steps: int
    steps: list[KineticStep]
    flatness_score: float
    """E_a spread across α normalized by the median.

    < 0.05  → effectively single-step;  0.05–0.20 → mild multi-step;
    > 0.20 → clear multi-step.
    """


def detect_steps(
    iso: IsoconversionalResult,
    *,
    relative_jump_threshold: float = 0.10,
    min_segment_size: int = 3,
) -> MultiStepResult:
    """Detect contiguous α-segments where E_a is approximately constant.

    Parameters
    ----------
    iso :
        Isoconversional result (preferably Vyazovkin or Vyazovkin-AIC for
        a clean E(α) curve).
    relative_jump_threshold :
        A point starts a new segment when its E differs from the
        running-segment median by more than this fraction.
    min_segment_size :
        Smallest acceptable segment length in α-points; trailing
        too-small segments get merged into the preceding one.
    """
    Ea = iso.Ea_kJ_per_mol
    a = iso.alpha
    valid = np.isfinite(Ea)
    if valid.sum() < 2:
        raise ValueError("Need at least 2 finite E(α) points to segment.")

    Ea = Ea[valid]
    a = a[valid]
    median_E = float(np.median(Ea))
    flatness = float((np.max(Ea) - np.min(Ea)) / max(median_E, 1e-12))

    # Greedy segmentation: extend the current segment while new points
    # stay close to the segment's running median.
    boundaries: list[int] = [0]
    seg_start = 0
    for i in range(1, len(Ea)):
        seg_median = float(np.median(Ea[seg_start : i]))
        if abs(Ea[i] - seg_median) > relative_jump_threshold * abs(seg_median):
            boundaries.append(i)
            seg_start = i
    boundaries.append(len(Ea))

    # Merge segments smaller than min_segment_size into a neighbour.
    while True:
        merged = False
        for k in range(len(boundaries) - 1, 0, -1):
            seg_size = boundaries[k] - boundaries[k - 1]
            if seg_size < min_segment_size and len(boundaries) > 2:
                # Remove the boundary that defines this small segment.
                # If it's the first segment, drop its right boundary; else
                # drop its left boundary so it merges with the segment before.
                drop = k if k == 1 else k - 1
                del boundaries[drop]
                merged = True
                break
        if not merged:
            break

    steps: list[KineticStep] = []
    for k in range(len(boundaries) - 1):
        lo, hi = boundaries[k], boundaries[k + 1]
        segment = Ea[lo:hi]
        seg_a = a[lo:hi]
        seg_med = float(np.median(segment))
        seg_mad = float(np.median(np.abs(segment - seg_med)))
        steps.append(
            KineticStep(
                index=k,
                alpha_lo=float(seg_a[0]),
                alpha_hi=float(seg_a[-1]),
                Ea_kJ_per_mol_median=seg_med,
                Ea_kJ_per_mol_mad=seg_mad,
                contribution=float(seg_a[-1] - seg_a[0]),
            )
        )

    return MultiStepResult(n_steps=len(steps), steps=steps, flatness_score=flatness)


__all__ = ["KineticStep", "MultiStepResult", "detect_steps"]

"""ICTAC consistency check across isoconversional methods.

ICTAC 2011/2020 recommend running at least two isoconversional methods
and reporting agreement: large pairwise disagreement is a red flag for
either bad data (insufficient β range, poor baseline) or a kinetic
pathology (multi-step process, distributed reactivity, diffusion control).

This module computes, for every pair (M_i, M_j) of methods, the maximum
relative difference between their E(α) curves:

    Δ_ij(α) = | E_i(α) − E_j(α) | / max(|E_i(α)|, |E_j(α)|, ε)
    Δ_ij    = max_α Δ_ij(α)

and the worst-case α at which the disagreement is observed.

A summary warning is emitted when any pair exceeds a configurable
threshold (default 10%, in line with ICTAC 2020 §3.2 wording about
"systematic differences").

References:
    ICTAC Kinetics Committee Recommendations, Thermochim. Acta 520 (2011)
    1-19, DOI 10.1016/j.tca.2011.03.034 — §6.
    Vyazovkin et al., Thermochim. Acta 689 (2020) 178597 — §3.
"""
from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations

import numpy as np

from .common import IsoconversionalResult


@dataclass(frozen=True)
class PairwiseDifference:
    method_a: str
    method_b: str
    max_relative_difference: float
    """Worst-case |E_a − E_b| / max(|E_a|, |E_b|) over the shared α grid."""

    alpha_of_max: float
    """α at which the worst-case difference occurs."""


@dataclass(frozen=True)
class ConsistencyResult:
    pairs: list[PairwiseDifference]
    threshold: float
    """Above this relative difference a pair is considered inconsistent."""

    warnings: list[str]
    """Human-readable per-pair warnings; empty if every pair is within threshold."""

    @property
    def max_pair_difference(self) -> float:
        return max((p.max_relative_difference for p in self.pairs), default=0.0)


def consistency_check(
    iso_results: dict[str, IsoconversionalResult],
    *,
    threshold: float = 0.10,
    eps_E: float = 1.0e3,  # 1 kJ/mol, in J/mol
) -> ConsistencyResult:
    """Pairwise agreement of E(α) across isoconversional methods.

    Parameters
    ----------
    iso_results :
        Map of method-name → result; expected to share the same α grid
        (true today because the runner uses a single grid for all methods).
    threshold :
        Relative-difference threshold above which a warning is emitted.
        0.10 is the ICTAC 2020 §3 "rough" value; tightening to 0.05 makes
        sense for high-quality data, 0.20 for noisy/limited-β datasets.
    eps_E :
        Floor on the denominator (in J/mol) so a near-zero E doesn't make
        the ratio explode for noisy α-endpoints.

    Returns
    -------
    ConsistencyResult
        One :class:`PairwiseDifference` for every unordered pair, plus any
        warnings triggered.
    """
    methods = sorted(iso_results.keys())
    pairs: list[PairwiseDifference] = []
    warnings: list[str] = []

    for a_name, b_name in combinations(methods, 2):
        ea = iso_results[a_name]
        eb = iso_results[b_name]

        # Restrict to the *intersection* of finite α values, in case some
        # method (e.g. Friedman) NaN'd a tail point.
        finite = np.isfinite(ea.Ea_J_per_mol) & np.isfinite(eb.Ea_J_per_mol)
        if not finite.any():
            pairs.append(PairwiseDifference(a_name, b_name, float("nan"), float("nan")))
            continue

        alpha = ea.alpha[finite]  # assumes shared α grid
        E_a = ea.Ea_J_per_mol[finite]
        E_b = eb.Ea_J_per_mol[finite]
        denom = np.maximum(np.maximum(np.abs(E_a), np.abs(E_b)), eps_E)
        rel = np.abs(E_a - E_b) / denom

        idx = int(np.argmax(rel))
        pairs.append(
            PairwiseDifference(
                method_a=a_name,
                method_b=b_name,
                max_relative_difference=float(rel[idx]),
                alpha_of_max=float(alpha[idx]),
            )
        )

        if rel[idx] > threshold:
            warnings.append(
                f"{a_name} vs {b_name}: max ΔE/E = {rel[idx] * 100:.1f}% "
                f"at α={alpha[idx]:.2f} (threshold {threshold * 100:.0f}%)"
            )

    return ConsistencyResult(pairs=pairs, threshold=threshold, warnings=warnings)


__all__ = ["ConsistencyResult", "PairwiseDifference", "consistency_check"]

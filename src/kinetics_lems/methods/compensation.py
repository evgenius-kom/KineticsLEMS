"""Kinetic compensation effect — ln A vs E linearity diagnostic.

Across α (isoconversional A) and across reaction-model assumptions
(Coats-Redfern) the pairs ``(E_i, ln A_i)`` often fall on a straight line:

    ln A  ≈  a  +  b · E.                                  (compensation law)

This is *not* a physical law but an artefact of the Arrhenius
linearization at a common temperature: it tells you that ``E`` and
``ln A`` are strongly correlated and *not* independently identifiable.
A high R² of the compensation line is therefore a red flag, not a
validation — the kinetic triplet may be statistically degenerate.

Reference:
    Vyazovkin & Wight (2000), *Thermochim. Acta* 340, 53 — discussion
    of the compensation effect as an identifiability marker;
    Galwey & Brown (2002), *Thermochim. Acta* 386, 91.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .coats_redfern import CoatsRedfernResult
from .common import linear_regression


@dataclass(frozen=True)
class CompensationFit:
    """Linear fit ln A ≈ a + b · E."""

    source: str
    """Where the (E, ln A) points came from — "coats_redfern" or "isoconversional"."""

    E_kJ_per_mol: np.ndarray
    ln_A_per_sec: np.ndarray
    slope: float
    """``b`` — the compensation slope (1 / (kJ/mol))."""

    intercept: float
    """``a`` — the compensation intercept (dimensionless)."""

    r_squared: float
    n_points: int


def compensation_from_coats_redfern(
    cr_result: CoatsRedfernResult,
) -> CompensationFit:
    """Fit ln A vs E across the *per-(model, run) fits* of Coats-Redfern.

    Each individual fit contributes one ``(E, ln A)`` pair, so all 12
    canonical f(α) models and every heating rate participate. A tight
    line here is the classical "compensation effect" signature
    described in Vyazovkin & Wight (2000).
    """
    Es = np.array(
        [f.Ea_kJ_per_mol for f in cr_result.fits if f.A_per_sec > 0], dtype=float
    )
    lnA = np.array(
        [np.log(f.A_per_sec) for f in cr_result.fits if f.A_per_sec > 0], dtype=float
    )
    if Es.size < 2:
        raise ValueError("Need at least 2 (E, A) points to fit compensation")
    slope, intercept, r2 = linear_regression(Es, lnA)
    return CompensationFit(
        source="coats_redfern",
        E_kJ_per_mol=Es,
        ln_A_per_sec=lnA,
        slope=float(slope),
        intercept=float(intercept),
        r_squared=float(r2),
        n_points=int(Es.size),
    )


def compensation_from_isoconversional(
    Ea_kJ_per_mol: np.ndarray,
    A_per_sec: np.ndarray,
) -> CompensationFit:
    """Fit ln A vs E across α from the isoconversional kinetic triplet.

    ``Ea_kJ_per_mol`` and ``A_per_sec`` must come from the same α grid
    (e.g. Vyazovkin + the per-α A returned by ``compute_A``).
    """
    finite = np.isfinite(Ea_kJ_per_mol) & np.isfinite(A_per_sec) & (A_per_sec > 0)
    Es = np.asarray(Ea_kJ_per_mol, dtype=float)[finite]
    As = np.asarray(A_per_sec, dtype=float)[finite]
    if Es.size < 2:
        raise ValueError("Need at least 2 finite (E, A) points to fit compensation")
    lnA = np.log(As)
    slope, intercept, r2 = linear_regression(Es, lnA)
    return CompensationFit(
        source="isoconversional",
        E_kJ_per_mol=Es,
        ln_A_per_sec=lnA,
        slope=float(slope),
        intercept=float(intercept),
        r_squared=float(r2),
        n_points=int(Es.size),
    )


__all__ = [
    "CompensationFit",
    "compensation_from_coats_redfern",
    "compensation_from_isoconversional",
]

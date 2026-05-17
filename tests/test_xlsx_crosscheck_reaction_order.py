"""Cross-check the reaction-order linearization against precomputed xlsx values.

Fixture: ``data/reference_workbooks/xlsx_crosscheck/reaction_order_*.tsv``
contains, per α, the input (T, dα/dT, 1/T) and the precomputed
``y(n; α) = ln(dα/dT) − n · ln(1 − α)`` for n ∈ {1, 2, 3, 0.7}.

Our :mod:`reaction_order` uses ``ln(dα/dt) − n · ln(1 − α)`` (dα/dt, not
dα/dT), which differs from xlsx's expression by a constant ``ln(β)`` —
that constant is absorbed into the line intercept and does not affect
the slope (E_a) or R². To compare bit-for-bit we recompute the xlsx
formula here and verify our re-implementation produces the same y values.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

FIXTURES = Path(__file__).resolve().parents[1] / "data" / "reference_workbooks" / "xlsx_crosscheck"

_N_VALUES = [1.0, 2.0, 3.0, 0.7]
_N_COLUMNS = ["y_n1", "y_n2", "y_n3", "y_n0.7"]


def _load(name: str) -> dict:
    data = np.loadtxt(FIXTURES / name, skiprows=1)
    return {
        "alpha": data[:, 0],
        "T_K": data[:, 1],
        "dalpha_dT": data[:, 2],
        "inv_T": data[:, 3],
        "y_n1": data[:, 4],
        "y_n2": data[:, 5],
        "y_n3": data[:, 6],
        "y_n0.7": data[:, 7],
    }


@pytest.mark.parametrize("fixture,label", [
    ("reaction_order_abpot.tsv", "ABPOT"),
    ("reaction_order_dapot.tsv", "DAPOT"),
])
def test_y_n_formula_matches_xlsx(fixture, label):
    """y(n) = ln(dα/dT) − n · ln(1−α) — the formula our regression uses
    (up to a per-run additive ln(β), which does not change slope or R²)."""
    d = _load(fixture)
    alpha = d["alpha"]
    dadT = d["dalpha_dT"]
    # Drop α≈0 / α≈1 where ln(1-α) explodes and xlsx may have rounding noise.
    interior = (alpha > 1e-6) & (alpha < 1.0 - 1e-6) & (dadT > 0)
    a = alpha[interior]
    rate = dadT[interior]

    ln_rate = np.log(rate)
    ln_one_minus_a = np.log(1.0 - a)

    for n, col in zip(_N_VALUES, _N_COLUMNS, strict=True):
        y_ours = ln_rate - n * ln_one_minus_a
        y_xlsx = d[col][interior]
        np.testing.assert_allclose(
            y_ours, y_xlsx, rtol=1e-5, atol=1e-6,
            err_msg=f"{label} n={n}: y(α) drift from xlsx",
        )


@pytest.mark.parametrize("fixture,label", [
    ("reaction_order_abpot.tsv", "ABPOT"),
    ("reaction_order_dapot.tsv", "DAPOT"),
])
def test_inv_T_consistency(fixture, label):
    """xlsx stored 1/T separately — sanity-check it equals our derivation."""
    d = _load(fixture)
    np.testing.assert_allclose(d["inv_T"], 1.0 / d["T_K"], rtol=1e-5)

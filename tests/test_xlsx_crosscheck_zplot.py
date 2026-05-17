"""Cross-check our master_plot Z(α) against the precomputed xlsx Z column.

Fixture: ``data/reference_workbooks/xlsx_crosscheck/zplot_*_5kpm.tsv``
contains the (α, T, dα/dt) input data and the *normalized* Z(α)/Z(0.5)
column produced by the original Excel formulas.

We build a single-rate :class:`ConversionRun` from the xlsx (α, T, dα/dt)
and run :func:`experimental_z_alpha`; the result must match the xlsx
column to floating-point precision. Any disagreement means our
implementation drifted from the formula that motivated the feature.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from kinetics_lems.conversion import ConversionRun
from kinetics_lems.methods.master_plot import experimental_z_alpha

FIXTURES = Path(__file__).resolve().parents[1] / "data" / "reference_workbooks" / "xlsx_crosscheck"


def _load_zplot(name: str) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    data = np.loadtxt(FIXTURES / name, skiprows=1)
    return data[:, 0], data[:, 1], data[:, 2], data[:, 3]


@pytest.mark.parametrize(
    "fixture,label",
    [
        ("zplot_abpot_5kpm.tsv", "ABPOT"),
        ("zplot_dapot_5kpm.tsv", "DAPOT"),
    ],
)
def test_zplot_matches_xlsx(fixture, label):
    alpha, T, dadt, Z_xlsx = _load_zplot(fixture)
    run = ConversionRun(
        rate_K_per_min=5.0,
        temperature=T,
        alpha=alpha,
        dalpha_dt=dadt,
    )
    z_ours = experimental_z_alpha([run], alpha)
    # Both normalize to Z(0.5) = 1, so the curves should coincide.
    # Tolerance: the xlsx stores values truncated to ~10 sig figs;
    # our floats are ~16 sig figs, so 1e-6 relative is comfortable.
    np.testing.assert_allclose(
        z_ours, Z_xlsx, rtol=1e-6, atol=1e-8,
        err_msg=f"{label}: experimental Z(α) drift from xlsx",
    )


def test_z_norm_equals_one_at_alpha_half():
    """Sanity check on the fixture itself: xlsx normalizes at α≈0.5."""
    for fixture in ("zplot_abpot_5kpm.tsv", "zplot_dapot_5kpm.tsv"):
        alpha, _, _, Z_xlsx = _load_zplot(fixture)
        idx = int(np.argmin(np.abs(alpha - 0.5)))
        assert Z_xlsx[idx] == pytest.approx(1.0, abs=1e-6), (
            f"{fixture}: Z_xlsx near α=0.5 should be 1, got {Z_xlsx[idx]}"
        )

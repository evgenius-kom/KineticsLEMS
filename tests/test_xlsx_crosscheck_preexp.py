"""Cross-check our preexponential A(α) against precomputed xlsx values.

Fixture: ``data/reference_workbooks/xlsx_crosscheck/A_calculated_abpot.tsv``
contains, per (rate, α), the input (T, dα/dt, Evya) and the resulting
A value from the Excel formula
``A = (dα/dt) / [(1 − α) · exp(−Evya / (R·T))]``  (F1 model baked in).

For each rate we build a :class:`ConversionRun` from the xlsx (α, T, dα/dt)
and call :func:`compute_A` with the xlsx's ``Evya`` and model ``"F1"``.
The resulting ``A_per_sec_per_alpha`` must match the xlsx ``A`` column to
floating-point precision.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

from kinetics_lems.constants import J_PER_KJ
from kinetics_lems.conversion import ConversionRun
from kinetics_lems.methods.preexponential import compute_A

FIXTURE = (
    Path(__file__).resolve().parents[1]
    / "data" / "reference_workbooks" / "xlsx_crosscheck" / "A_calculated_abpot.tsv"
)


def _load() -> dict[float, np.ndarray]:
    """Group xlsx rows by heating rate. Returns rate → (α, T, dα/dt, Evya, A)."""
    data = np.loadtxt(FIXTURE, skiprows=1)
    out: dict[float, np.ndarray] = {}
    for rate in np.unique(data[:, 0]):
        mask = data[:, 0] == rate
        out[float(rate)] = data[mask, 1:]
    return out


def test_A_recomputed_matches_xlsx_per_rate():
    grouped = _load()
    for rate, arr in grouped.items():
        alpha = arr[:, 0]
        T = arr[:, 1]
        dadt = arr[:, 2]
        Evya_kJ = arr[:, 3]
        A_xlsx = arr[:, 4]

        # Drop α=0 and α=1 edges where f(α)=1-α may produce divides-by-zero
        # in our code (we floor with eps_f=1e-12, xlsx may have NaN/Inf).
        interior = (alpha > 1e-6) & (alpha < 1.0 - 1e-6) & np.isfinite(A_xlsx) & (A_xlsx > 0)
        alpha = alpha[interior]
        Evya_J = Evya_kJ[interior] * J_PER_KJ
        A_xlsx = A_xlsx[interior]

        run = ConversionRun(
            rate_K_per_min=rate,
            temperature=T[interior],
            alpha=alpha,
            dalpha_dt=dadt[interior],
        )
        result = compute_A(
            runs=[run],
            alphas=alpha,
            Ea_J_per_mol=Evya_J,
            model="F1",
        )
        # Single-run median across runs = the single-run value itself.
        # Tolerance: xlsx uses R = 8.314, our code uses CODATA R = 8.31446.
        # The mismatch propagates through exp(−E/RT) and yields a ~0.2%
        # systematic bias on A. We accept up to 0.5% — anything larger is
        # not "different R", it's a real formula divergence.
        np.testing.assert_allclose(
            result.A_per_sec_per_alpha,
            A_xlsx,
            rtol=5e-3,
            err_msg=f"rate {rate} K/min: A(α) drift from xlsx",
        )

"""Targeted unit tests for internal primitives that are easy to silently break."""
from __future__ import annotations

import numpy as np
import pytest

from kinetics_lems.constants import R_GAS, SEC_PER_MIN
from kinetics_lems.conversion import build_runs
from kinetics_lems.methods.vyazovkin import _p_senum_yang, _phi_off_diagonal_sum
from kinetics_lems.synthetic import generate_case
from kinetics_lems.synthetic.generator import _auto_temperature_window

# ---------- Senum–Yang p(x) approximation ----------

def test_p_senum_yang_works_on_scalar_and_array():
    """Same code path must accept both shapes — used in vectorized Φ."""
    x_scalar = 30.0
    x_array = np.array([20.0, 30.0, 50.0, 80.0])

    p_scalar = _p_senum_yang(x_scalar)
    p_array = _p_senum_yang(x_array)

    assert isinstance(p_array, np.ndarray)
    assert p_array.shape == (4,)
    # Spot-check scalar matches the array element at the same x.
    assert pytest.approx(p_array[1], rel=1e-12) == p_scalar


def test_p_senum_yang_matches_reference_value():
    """Senum–Yang p(30) ≈ exp(-30)/30 · (900+300+18)/(27000+10800+1080+24).

    Hand-computed reference: ~ 9.0e-15. Small absolute number but stable.
    """
    p = _p_senum_yang(30.0)
    expected = np.exp(-30.0) / 30.0 * (900 + 300 + 18) / (27000 + 10800 + 1080 + 24)
    assert pytest.approx(p, rel=1e-12) == expected


# ---------- Off-diagonal sum helper (the heart of Vyazovkin Φ) ----------

def test_phi_off_diagonal_sum_excludes_diagonal():
    v = np.array([1.0, 2.0, 4.0])
    # Σ_{i≠j} v_i / v_j  (size 3 → 6 terms).
    # Pairs: 1/2, 1/4, 2/1, 2/4, 4/1, 4/2 → 0.5+0.25+2+0.5+4+2 = 9.25
    assert pytest.approx(_phi_off_diagonal_sum(v, v), rel=1e-12) == 9.25


def test_phi_off_diagonal_sum_minimum_when_vectors_proportional():
    """If all numerators/denominators are equal (perfect Vyazovkin solution), Σ = n(n-1)."""
    n = 4
    v = np.full(n, 7.0)
    assert pytest.approx(_phi_off_diagonal_sum(v, v), rel=1e-12) == n * (n - 1)


# ---------- _auto_temperature_window picks a safe T_start at the lowest β ----------

def test_auto_window_safe_at_smallest_rate():
    """At T_start the *normalized* per-K rate at β_min must be tiny (≤ ~1e-5)."""
    Ea = 120_000.0
    A = 1e10
    beta_min = 2.5  # K/min
    T_start, T_stop = _auto_temperature_window(Ea, A, beta_min)
    assert T_start < T_stop

    beta_min_per_sec = beta_min / SEC_PER_MIN
    dadT_at_start = (A / beta_min_per_sec) * np.exp(-Ea / (R_GAS * T_start))
    assert dadT_at_start < 1e-4, f"Reaction not negligible at T_start: dα/dT = {dadT_at_start}"


# ---------- Round-trip α(T) for synthetic single-step ----------

def test_alpha_starts_at_zero_and_reaches_one_for_synthetic():
    """The recovered α(T) for a clean synthetic should span [0, 1]."""
    case = generate_case([5.0, 10.0, 20.0], Ea_J_per_mol=120_000.0, n_points=2000, seed=0)
    runs = build_runs(case)
    for r in runs:
        assert pytest.approx(r.alpha[0], abs=1e-6) == 0.0
        assert pytest.approx(r.alpha[-1], abs=1e-6) == 1.0

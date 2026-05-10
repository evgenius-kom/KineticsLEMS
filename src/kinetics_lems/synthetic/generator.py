"""Generate synthetic DSC-like peaks from a known kinetic model.

For a single-step Arrhenius reaction with model f(α) under linear heating β:

    dα/dt = A · f(α) · exp(-E / (R·T)),   T(t) = T₀ + β·t.

We integrate this ODE on a fine time grid for each β and convert dα/dt back
into a "DSC heat-flow" curve y = dα/dT (proportional to heat flow, peak shape
preserved). Writing y vs T gives a wave that the loader can read identically
to a real .txt experimental file.
"""
from __future__ import annotations

import json
from collections.abc import Callable
from pathlib import Path

import numpy as np

from ..constants import R_GAS, SEC_PER_MIN
from ..models import CaseData, CaseParams, ExperimentType, Method, Wave


def _f_first_order(a: float) -> float:
    return max(1.0 - a, 0.0)


REACTION_MODELS: dict[str, Callable[[float], float]] = {
    "F1": _f_first_order,                                # first-order: 1 − α
    "F2": lambda a: max(1.0 - a, 0.0) ** 2,              # second-order
    "A2": lambda a: 2.0 * (1.0 - a) * np.sqrt(max(-np.log(max(1.0 - a, 1e-12)), 0.0)),  # Avrami n=2
    "R2": lambda a: 2.0 * np.sqrt(max(1.0 - a, 0.0)),    # contracting cylinder
    "R3": lambda a: 3.0 * (max(1.0 - a, 0.0)) ** (2.0 / 3.0),  # contracting sphere
}


def generate_case(
    rates_K_per_min: list[float],
    *,
    Ea_J_per_mol: float = 100_000.0,
    A_per_sec: float = 1.0e10,
    model: str = "F1",
    T_start: float | None = None,
    T_stop: float | None = None,
    n_points: int = 1500,
    noise_std: float = 0.0,
    seed: int | None = None,
    material: str = "synthetic",
) -> CaseData:
    """Generate a CaseData with one wave per heating rate.

    Parameters
    ----------
    rates_K_per_min : list[float]
        Heating rates β in K/min.
    Ea_J_per_mol, A_per_sec, model :
        Ground-truth Arrhenius parameters and reaction model.
    T_start, T_stop :
        Temperature window of the simulated experiment in K.
        When ``None``, both are picked automatically from the kinetics so the
        peak fits well inside the window with negligible rate at the edges.
    n_points :
        Number of T samples per wave.
    noise_std :
        Optional Gaussian noise stddev added to dα/dT (relative).
    """
    if model not in REACTION_MODELS:
        raise ValueError(f"Unknown model '{model}'. Available: {list(REACTION_MODELS)}")
    f = REACTION_MODELS[model]
    rng = np.random.default_rng(seed)

    auto_start, auto_stop = _auto_temperature_window(
        Ea_J_per_mol, A_per_sec, min(rates_K_per_min)
    )
    T_start_eff = float(T_start) if T_start is not None else auto_start
    T_stop_eff = float(T_stop) if T_stop is not None else auto_stop

    waves: dict[float, Wave] = {}
    T_grid_template = np.linspace(T_start_eff, T_stop_eff, n_points)
    for beta_min in rates_K_per_min:
        beta_sec = beta_min / SEC_PER_MIN
        T_grid = T_grid_template.copy()
        # Integrate dα/dT = A/β · f(α) · exp(-E/(R·T)) by simple RK4 in T.
        alpha = np.zeros_like(T_grid)
        for k in range(1, T_grid.size):
            T0 = T_grid[k - 1]
            T1 = T_grid[k]
            dT = T1 - T0
            a = alpha[k - 1]
            k1 = _rate_T(a, T0, A_per_sec, beta_sec, Ea_J_per_mol, f)
            k2 = _rate_T(a + 0.5 * dT * k1, 0.5 * (T0 + T1), A_per_sec, beta_sec, Ea_J_per_mol, f)
            k3 = _rate_T(a + 0.5 * dT * k2, 0.5 * (T0 + T1), A_per_sec, beta_sec, Ea_J_per_mol, f)
            k4 = _rate_T(a + dT * k3, T1, A_per_sec, beta_sec, Ea_J_per_mol, f)
            alpha[k] = min(a + dT * (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0, 1.0)

        # "Heat-flow"-shaped wave proportional to dα/dT.
        dadT = np.gradient(alpha, T_grid)
        if noise_std > 0:
            dadT = dadT * (1.0 + rng.normal(0.0, noise_std, size=dadT.shape))
        waves[beta_min] = Wave(T_grid, dadT)

    file_map = {
        f"rate_{i:02d}.txt": float(b) for i, b in enumerate(sorted(rates_K_per_min))
    }
    params = CaseParams(
        material=material,
        experiment_type=ExperimentType.HEATING,
        method=Method.DSC,
        file_to_condition=file_map,
    )
    return CaseData(params=params, waves=waves)


def _auto_temperature_window(
    Ea_J_per_mol: float,
    A_per_sec: float,
    beta_min_K_per_min: float,
) -> tuple[float, float]:
    """Pick a T window where the rate at the edges is negligible.

    Picks T_start so that dα/dT = (A/β)·exp(-E/RT) ≪ 1 even at the *slowest*
    heating rate (β_min) — which is the worst case, since lower β makes
    dα/dT larger at any given T. T_stop is set well above the peak so the
    reaction completes inside the window.
    """
    k_low = 1e-6  # target k(T_start) = k_low · β_min  →  dα/dT ≈ k_low at T_start
    k_high = 1e2  # 1/s; reaction effectively over above this
    beta_min_per_sec = beta_min_K_per_min / SEC_PER_MIN

    def T_for_rate(k: float) -> float:
        return Ea_J_per_mol / (R_GAS * np.log(A_per_sec / k))

    T_start = T_for_rate(k_low * beta_min_per_sec)
    T_stop = T_for_rate(k_high)
    # Small symmetric margin.
    span = T_stop - T_start
    return float(T_start - 0.05 * span), float(T_stop + 0.05 * span)


def _rate_T(
    alpha: float,
    T: float,
    A_per_sec: float,
    beta_sec: float,
    Ea: float,
    f: Callable[[float], float],
) -> float:
    if alpha >= 1.0:
        return 0.0
    return (A_per_sec / beta_sec) * f(alpha) * float(np.exp(-Ea / (R_GAS * T)))


def write_case(case: CaseData, output_dir: Path | str) -> Path:
    """Write a CaseData out as <dir>/settings.json + per-rate .txt files."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Re-key the file-to-condition map to actual files we will write.
    file_to_cond: dict[str, float] = {}
    for filename, condition in case.params.file_to_condition.items():
        wave = case.waves[float(condition)]
        path = output_dir / filename
        np.savetxt(path, np.column_stack([wave.x, wave.y]), fmt="%.6f", delimiter="\t")
        file_to_cond[filename] = float(condition)

    settings = {
        "Settings": [
            {"name": "Material", "value": case.params.material},
            {"name": "ExperimentType", "value": case.params.experiment_type.value},
            {"name": "Method", "value": case.params.method.value},
            {"name": "Conditions", "value": {k: str(v) for k, v in file_to_cond.items()}},
        ]
    }
    (output_dir / "settings.json").write_text(json.dumps(settings, indent=2), encoding="utf-8")
    return output_dir


__all__ = ["generate_case", "write_case", "REACTION_MODELS"]

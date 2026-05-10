"""Criado–Málek master plot Z(α) — model-discrimination via shape comparison.

Z(α) = f(α) · g(α), where f and g are the differential and integral
reaction-model functions. Empirically Z(α) ∝ (dα/dt) · T² (up to constants
that cancel after normalising at α = 0.5).

Plotting experimental Z(α)/Z(0.5) vs the master curves of standard models
identifies the f(α) that best describes the kinetics.

**Known degeneracy.** First-order F1 and any Avrami A_m share the same
*normalized* Z(α) shape because Z = m · (1−α) · [−ln(1−α)] for every A_m,
which differs from F1 only by a constant prefactor that cancels after
normalising at α = 0.5. The Z-plot can rank "F1-class" vs "F_n with n>1"
vs "R_n" vs "D_n", but it cannot separate F1 from A_m by itself —
combine with Friedman intercepts or master-plot of y(α) for that.

Reference: Criado, J. M. *Thermochim. Acta* 24 (1978) 186; Málek, J.
*Thermochim. Acta* 200 (1992) 257; ICTAC 2011 §5.
"""
from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass

import numpy as np

from ..conversion import ConversionRun, dalpha_dt_at_conversion, temperature_at_conversion

# ---------- Master curves ----------

def _f_F(n: float) -> Callable[[np.ndarray], np.ndarray]:
    """f(α) for F_n: (1 − α)^n."""
    return lambda a: np.power(1.0 - a, n)


def _g_F(n: float) -> Callable[[np.ndarray], np.ndarray]:
    """g(α) = ∫₀^α dα'/(1−α')^n.

    For n=1 this is −ln(1−α). For n≠1 the closed form is
    ((1−α)^(1−n) − 1) / (n−1) = (1 − (1−α)^(1−n)) / (1−n);
    both denominators give the same sign as the integrand (positive).
    """
    if abs(n - 1.0) < 1e-12:
        return lambda a: -np.log(np.maximum(1.0 - a, 1e-300))
    exp = 1.0 - n
    return lambda a: (np.power(np.maximum(1.0 - a, 1e-300), exp) - 1.0) / (n - 1.0)


def _f_A(m: float) -> Callable[[np.ndarray], np.ndarray]:
    """Avrami f(α): m·(1−α)·[−ln(1−α)]^((m−1)/m)."""
    return lambda a: m * (1.0 - a) * np.power(
        np.maximum(-np.log(np.maximum(1.0 - a, 1e-300)), 0.0), (m - 1.0) / m
    )


def _g_A(m: float) -> Callable[[np.ndarray], np.ndarray]:
    """Avrami g(α): [−ln(1−α)]^(1/m)."""
    return lambda a: np.power(
        np.maximum(-np.log(np.maximum(1.0 - a, 1e-300)), 0.0), 1.0 / m
    )


def _f_R(n: float) -> Callable[[np.ndarray], np.ndarray]:
    """Contracting-geometry f(α): n·(1−α)^((n−1)/n) for R_n."""
    return lambda a: n * np.power(np.maximum(1.0 - a, 0.0), (n - 1.0) / n)


def _g_R(n: float) -> Callable[[np.ndarray], np.ndarray]:
    """Contracting-geometry g(α): 1 − (1−α)^(1/n) for R_n."""
    return lambda a: 1.0 - np.power(np.maximum(1.0 - a, 0.0), 1.0 / n)


# Diffusion models (Málek 1992 Table 1).
def _f_D1(a: np.ndarray) -> np.ndarray:
    return 1.0 / (2.0 * np.maximum(a, 1e-12))


def _g_D1(a: np.ndarray) -> np.ndarray:
    return a * a


def _f_D2(a: np.ndarray) -> np.ndarray:
    # 1 / [-ln(1-α)] (Valensi 2-D diffusion)
    return 1.0 / np.maximum(-np.log(np.maximum(1.0 - a, 1e-300)), 1e-12)


def _g_D2(a: np.ndarray) -> np.ndarray:
    # (1-α) ln(1-α) + α
    return (1.0 - a) * np.log(np.maximum(1.0 - a, 1e-300)) + a


def _f_D3(a: np.ndarray) -> np.ndarray:
    # 3·(1-α)^(2/3) / (2·[1 - (1-α)^(1/3)])  (Jander 3-D diffusion)
    one_minus = np.maximum(1.0 - a, 0.0)
    denom = 2.0 * (1.0 - np.power(one_minus, 1.0 / 3.0))
    return 1.5 * np.power(one_minus, 2.0 / 3.0) / np.maximum(denom, 1e-12)


def _g_D3(a: np.ndarray) -> np.ndarray:
    return np.square(1.0 - np.power(np.maximum(1.0 - a, 0.0), 1.0 / 3.0))


def _f_D4(a: np.ndarray) -> np.ndarray:
    # 3 / [2·((1-α)^(-1/3) - 1)]  (Ginstling–Brounshtein)
    one_minus = np.maximum(1.0 - a, 1e-12)
    denom = 2.0 * (np.power(one_minus, -1.0 / 3.0) - 1.0)
    return 1.5 / np.maximum(denom, 1e-12)


def _g_D4(a: np.ndarray) -> np.ndarray:
    return 1.0 - 2.0 / 3.0 * a - np.power(np.maximum(1.0 - a, 0.0), 2.0 / 3.0)


@dataclass(frozen=True)
class ReactionModel:
    name: str
    f: Callable[[np.ndarray], np.ndarray]  # f(α)
    g: Callable[[np.ndarray], np.ndarray]  # g(α) (integral form)

    def z(self, alpha: np.ndarray) -> np.ndarray:
        """Theoretical Z(α) = f(α) · g(α), unnormalized."""
        return self.f(alpha) * self.g(alpha)


# 12-model standard set (Málek 1992 / ICTAC 2011 §5).
MASTER_MODELS: dict[str, ReactionModel] = {
    "F1": ReactionModel("F1", _f_F(1.0), _g_F(1.0)),
    "F2": ReactionModel("F2", _f_F(2.0), _g_F(2.0)),
    "F3": ReactionModel("F3", _f_F(3.0), _g_F(3.0)),
    "A2": ReactionModel("A2", _f_A(2.0), _g_A(2.0)),
    "A3": ReactionModel("A3", _f_A(3.0), _g_A(3.0)),
    "A4": ReactionModel("A4", _f_A(4.0), _g_A(4.0)),
    "R2": ReactionModel("R2", _f_R(2.0), _g_R(2.0)),
    "R3": ReactionModel("R3", _f_R(3.0), _g_R(3.0)),
    "D1": ReactionModel("D1", _f_D1, _g_D1),
    "D2": ReactionModel("D2", _f_D2, _g_D2),
    "D3": ReactionModel("D3", _f_D3, _g_D3),
    "D4": ReactionModel("D4", _f_D4, _g_D4),
}


# ---------- Experimental Z(α) ----------

def experimental_z_alpha(
    runs: list[ConversionRun],
    alphas: np.ndarray,
) -> np.ndarray:
    """Z(α) averaged across runs.

    For each run k compute Z_k(α) = (dα/dt)_α · T_α², normalize by Z_k(0.5),
    then average across runs. The 0.5 normalization makes the curves
    comparable to the master curves (every master curve is normalized
    similarly so that Z_master(0.5) = 1).
    """
    z_per_run = []
    half_idx = int(np.argmin(np.abs(alphas - 0.5)))
    for r in runs:
        T = temperature_at_conversion(r, alphas)
        rate = dalpha_dt_at_conversion(r, alphas)
        z = rate * T * T
        z_norm = z / z[half_idx] if z[half_idx] > 0 else z
        z_per_run.append(z_norm)
    return np.mean(np.vstack(z_per_run), axis=0)


def master_z_curves(alphas: np.ndarray) -> dict[str, np.ndarray]:
    """Return the 12 master Z(α) curves, each normalized so that Z(0.5) = 1."""
    half_idx = int(np.argmin(np.abs(alphas - 0.5)))
    out: dict[str, np.ndarray] = {}
    for name, model in MASTER_MODELS.items():
        z = model.z(alphas)
        ref = z[half_idx]
        out[name] = z / ref if ref > 0 else z
    return out


# ---------- Ranking ----------

@dataclass(frozen=True)
class ModelRanking:
    alphas: np.ndarray
    experimental_z: np.ndarray
    master_curves: dict[str, np.ndarray]
    rms_distance: dict[str, float]  # smaller = better fit

    @property
    def best_model(self) -> str:
        return min(self.rms_distance, key=lambda n: self.rms_distance[n])

    def ranked(self) -> list[tuple[str, float]]:
        return sorted(self.rms_distance.items(), key=lambda kv: kv[1])


def rank_models(
    runs: list[ConversionRun],
    alphas: np.ndarray,
) -> ModelRanking:
    """Rank reaction models by RMS distance from experimental Z(α)."""
    if alphas[0] > 0.5 or alphas[-1] < 0.5:
        raise ValueError("Z(α) ranking needs an α grid that includes 0.5")
    z_exp = experimental_z_alpha(runs, alphas)
    z_master = master_z_curves(alphas)
    distances = {
        name: float(np.sqrt(np.mean((z_exp - curve) ** 2)))
        for name, curve in z_master.items()
    }
    return ModelRanking(
        alphas=alphas,
        experimental_z=z_exp,
        master_curves=z_master,
        rms_distance=distances,
    )


__all__ = [
    "MASTER_MODELS",
    "ModelRanking",
    "ReactionModel",
    "experimental_z_alpha",
    "master_z_curves",
    "rank_models",
]

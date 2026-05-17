"""Empirical / autocatalytic reaction models (Sestak-Berggren, Prout-Tompkins).

The two models supported here are 2-3 parameter generalisations of the
canonical f(α) family:

* **Prout-Tompkins**:    f(α) = α^m · (1 − α)^n.
* **Sestak-Berggren**:   f(α) = α^m · (1 − α)^n · [−ln(1 − α)]^p.

They are *empirical fits*, not mechanistic models — the parameters
(m, n, p) do not directly map onto physical processes. Use them when the
12 canonical master-plot models all fit poorly (typical for autocatalytic
cure kinetics and other sigmoidal reactions).

Both are fitted by minimising the RMS distance between the experimental
Z(α) (Criado–Málek normalisation at α = 0.5) and the analytical Z_emp(α)
generated from f(α) and its numerical integral g(α). Returning the
``ReactionModel`` so the fitted shape can be slotted directly into the
existing master-plot ranking, A computation, and lifetime prediction.

Warning: Sestak-Berggren is famously *non-identifiable* — many (m, n, p)
combinations produce nearly identical Z(α). Always report R², parameter
covariance, and prefer fewer parameters when in doubt.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.optimize import least_squares

from ..conversion import ConversionRun
from .master_plot import ReactionModel, experimental_z_alpha


@dataclass(frozen=True)
class EmpiricalFit:
    """Best-fit empirical reaction model + diagnostics."""

    name: str  # "Prout-Tompkins" or "Sestak-Berggren"
    parameters: dict[str, float]
    model: ReactionModel  # the fitted f(α), g(α) — pluggable elsewhere
    rms: float
    """RMS distance from experimental Z(α), comparable to master_plot.rms_distance."""

    r_squared: float


def _safe_pow(base: np.ndarray, exp: float, *, floor: float = 1e-300) -> np.ndarray:
    """``base ** exp`` with a floor on ``base`` so 0^positive stays well-defined."""
    return np.power(np.maximum(base, floor), exp)


def _prout_tompkins_f(m: float, n: float):
    """f(α) = α^m · (1 − α)^n. Vectorised."""

    def f(a: np.ndarray) -> np.ndarray:
        return _safe_pow(a, m) * _safe_pow(1.0 - a, n)

    return f


def _sestak_berggren_f(m: float, n: float, p: float):
    """f(α) = α^m · (1 − α)^n · [−ln(1 − α)]^p."""

    def f(a: np.ndarray) -> np.ndarray:
        one_minus = np.maximum(1.0 - a, 1e-300)
        log_term = np.maximum(-np.log(one_minus), 1e-300)
        return _safe_pow(a, m) * _safe_pow(one_minus, n) * _safe_pow(log_term, p)

    return f


def _g_from_f(f, alphas: np.ndarray) -> np.ndarray:
    """g(α) = ∫₀^α dα'/f(α') by cumulative trapezoid on the supplied grid."""
    fa = f(alphas)
    integrand = 1.0 / np.where(fa > 0, fa, np.nan)
    g = np.zeros_like(alphas)
    g[1:] = np.cumsum(0.5 * (integrand[1:] + integrand[:-1]) * np.diff(alphas))
    return g


def _build_reaction_model(name: str, f) -> ReactionModel:
    """Wrap a numerical f(α) into a ReactionModel with numerical g(α)."""

    def g(alphas: np.ndarray) -> np.ndarray:
        return _g_from_f(f, np.asarray(alphas, dtype=float))

    return ReactionModel(name=name, f=f, g=g)


def _z_normalised(f, alphas: np.ndarray) -> np.ndarray:
    """Z(α) = f(α)·g(α) normalised to Z(0.5) = 1."""
    z = f(alphas) * _g_from_f(f, alphas)
    half_idx = int(np.argmin(np.abs(alphas - 0.5)))
    ref = z[half_idx]
    return z / ref if ref > 0 else z


def _rms_against_experimental(
    z_exp: np.ndarray, z_model: np.ndarray
) -> tuple[float, float]:
    """Return ``(rms_distance, r_squared)`` between experimental and model Z(α)."""
    finite = np.isfinite(z_exp) & np.isfinite(z_model)
    if finite.sum() < 2:
        return float("nan"), float("nan")
    diff = z_exp[finite] - z_model[finite]
    rms = float(np.sqrt(np.mean(diff * diff)))
    ss_res = float(np.sum(diff * diff))
    ss_tot = float(np.sum((z_exp[finite] - np.mean(z_exp[finite])) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")
    return rms, r2


def fit_prout_tompkins(
    runs: list[ConversionRun],
    alphas: np.ndarray,
    *,
    initial: tuple[float, float] = (0.5, 1.0),
    bounds: tuple[tuple[float, float], tuple[float, float]] = ((0.0, 4.0), (0.0, 4.0)),
) -> EmpiricalFit:
    """Fit Prout-Tompkins f(α) = α^m (1−α)^n to Z(α) of ``runs``."""
    if alphas[0] > 0.5 or alphas[-1] < 0.5:
        raise ValueError("Z(α) fit needs an α grid that includes 0.5")
    z_exp = experimental_z_alpha(runs, alphas)

    def residuals(params: np.ndarray) -> np.ndarray:
        m, n = params
        z_model = _z_normalised(_prout_tompkins_f(m, n), alphas)
        diff = z_exp - z_model
        return np.where(np.isfinite(diff), diff, 0.0)

    lb = (bounds[0][0], bounds[1][0])
    ub = (bounds[0][1], bounds[1][1])
    res = least_squares(residuals, x0=initial, bounds=(lb, ub))
    m_opt, n_opt = float(res.x[0]), float(res.x[1])
    f_opt = _prout_tompkins_f(m_opt, n_opt)
    rms, r2 = _rms_against_experimental(z_exp, _z_normalised(f_opt, alphas))
    return EmpiricalFit(
        name="Prout-Tompkins",
        parameters={"m": m_opt, "n": n_opt},
        model=_build_reaction_model(f"PT(m={m_opt:.3f}, n={n_opt:.3f})", f_opt),
        rms=rms,
        r_squared=r2,
    )


def fit_sestak_berggren(
    runs: list[ConversionRun],
    alphas: np.ndarray,
    *,
    initial: tuple[float, float, float] = (0.5, 1.0, 0.0),
    bounds: tuple[
        tuple[float, float], tuple[float, float], tuple[float, float]
    ] = ((0.0, 4.0), (0.0, 4.0), (-2.0, 2.0)),
) -> EmpiricalFit:
    """Fit Sestak-Berggren f(α) = α^m (1−α)^n [−ln(1−α)]^p to Z(α)."""
    if alphas[0] > 0.5 or alphas[-1] < 0.5:
        raise ValueError("Z(α) fit needs an α grid that includes 0.5")
    z_exp = experimental_z_alpha(runs, alphas)

    def residuals(params: np.ndarray) -> np.ndarray:
        m, n, p = params
        z_model = _z_normalised(_sestak_berggren_f(m, n, p), alphas)
        diff = z_exp - z_model
        return np.where(np.isfinite(diff), diff, 0.0)

    lb = (bounds[0][0], bounds[1][0], bounds[2][0])
    ub = (bounds[0][1], bounds[1][1], bounds[2][1])
    res = least_squares(residuals, x0=initial, bounds=(lb, ub))
    m_opt, n_opt, p_opt = float(res.x[0]), float(res.x[1]), float(res.x[2])
    f_opt = _sestak_berggren_f(m_opt, n_opt, p_opt)
    rms, r2 = _rms_against_experimental(z_exp, _z_normalised(f_opt, alphas))
    return EmpiricalFit(
        name="Sestak-Berggren",
        parameters={"m": m_opt, "n": n_opt, "p": p_opt},
        model=_build_reaction_model(
            f"SB(m={m_opt:.3f}, n={n_opt:.3f}, p={p_opt:.3f})", f_opt
        ),
        rms=rms,
        r_squared=r2,
    )


__all__ = [
    "EmpiricalFit",
    "fit_prout_tompkins",
    "fit_sestak_berggren",
]

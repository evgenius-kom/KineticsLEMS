"""Objective function components for model-based fitting.

The composite loss recommended by ICTAC 2020 §3 is

    L(θ) = w_α   · normalized_SSE(α)
         + w_dα  · normalized_SSE(dα/dt)
         + w_peak· peak-position penalty                   (optional)
         + regularization

normalising each component by its own data variance so the weights
``w_*`` carry intuitive meaning (1.0 = equal contribution). The split
matters because α is smooth (1/T-bounded RSS) while dα/dt is noisy and
peak-shaped — pure α-only loss tends to under-fit shoulder structure.

This file is infrastructure: the dataclasses are stable, the call site
:func:`Objective.evaluate` raises ``NotImplementedError`` until a
topology ODE is wired up in ``topology.py``.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np


@dataclass(frozen=True)
class ObjectiveWeights:
    """Relative weights of the loss components (all default to 1.0)."""

    alpha: float = 1.0
    rate: float = 1.0
    peak: float = 0.0


@dataclass(frozen=True)
class ObjectiveSpec:
    """Static description of which loss components are active and how."""

    weights: ObjectiveWeights = field(default_factory=ObjectiveWeights)
    use_log_rate: bool = False
    """Fit ``ln(dα/dt)`` instead of ``dα/dt`` — robust to peak-tail scale gap."""

    regularization_l2: float = 0.0


class Objective:
    """Composite loss evaluator. Stateless wrapper around a spec."""

    def __init__(self, spec: ObjectiveSpec) -> None:
        self.spec = spec

    def evaluate(self, alpha_pred, alpha_obs, rate_pred, rate_obs) -> float:
        """Compute the composite loss for a single run.

        TODO: implement once :func:`fitting.topology.build_ode_system`
        returns a working SINGLE-topology RHS. The math is straightforward
        — see this file's module docstring — but it deserves a clean
        per-run implementation rather than inline math at the optimiser
        boundary.
        """
        raise NotImplementedError(
            "Objective.evaluate is not yet implemented; see "
            "fitting/objective.py module docstring for the formula."
        )

    @staticmethod
    def normalised_sse(pred: np.ndarray, obs: np.ndarray) -> float:
        """Sum of squared residuals normalised by the variance of ``obs``."""
        denom = float(np.var(obs))
        if denom == 0.0:
            return float(np.sum((pred - obs) ** 2))
        return float(np.sum((pred - obs) ** 2) / denom)


__all__ = ["Objective", "ObjectiveSpec", "ObjectiveWeights"]

"""Fitting problem and result dataclasses.

These are *plain* dataclasses, not pydantic models, because they hold numpy
arrays in hot paths. The pydantic-friendly view lives in
:mod:`kinetics_lems.schemas` (``FitResultSchema``); a ``to_schema()``
adapter on :class:`FittingResult` is the eventual round-trip.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from ..conversion import ConversionRun
    from .objective import ObjectiveSpec
    from .topology import Topology


@dataclass(frozen=True)
class ParameterSpec:
    """One free parameter of a multi-step kinetic model.

    The fit is performed in *transformed* space (``log A`` instead of ``A``,
    soft-bounded weights via softmax) to keep optimisers happy; the
    ``bounds`` here are in the *natural* parameter space and the transform
    is applied internally — see :func:`Objective.transform_params`.
    """

    name: str
    initial: float
    bounds: tuple[float, float]
    unit: str = ""
    log_space: bool = False
    """If ``True``, the optimiser sees ``log10(value)``; useful for ``A``."""


@dataclass(frozen=True)
class FittingProblem:
    """A fit specification: data + model topology + objective + bounds.

    ``runs`` is the list of experimental conversion curves at different
    heating-rate programs (the same ``ConversionRun`` already produced by
    :func:`kinetics_lems.conversion.build_runs`).

    ``parameters`` lists the free parameters in fit order; the same order
    is used to unpack the parameter vector inside the topology ODE.

    ``holdout_indices`` lets you reserve experiments for validation; the
    optimiser only sees the complement, the diagnostics report both.
    """

    runs: list[ConversionRun]
    topology: Topology
    reaction_model: str  # "F1", "A2", "Fn", "SB", ...
    parameters: list[ParameterSpec]
    objective_spec: ObjectiveSpec
    holdout_indices: tuple[int, ...] = ()


@dataclass
class FittingResult:
    """Output of a fit.

    Fields are filled progressively by optimisers; missing ones are ``None``.
    The eventual JSON sidecar is obtained via
    ``FitResultSchema.model_validate(result.to_schema_dict())``.
    """

    problem: FittingProblem
    best_params: np.ndarray | None = None
    """Parameter vector in fit (transformed) space."""

    natural_params: np.ndarray | None = None
    """Same parameters in natural / human-readable space."""

    loss_total: float | None = None
    loss_components: dict[str, float] = field(default_factory=dict)
    residuals_per_run: list[np.ndarray] = field(default_factory=list)
    parameter_correlation: np.ndarray | None = None
    aic: float | None = None
    bic: float | None = None
    warnings: list[str] = field(default_factory=list)
    optimizer_info: dict[str, object] = field(default_factory=dict)


__all__ = ["FittingProblem", "FittingResult", "ParameterSpec"]

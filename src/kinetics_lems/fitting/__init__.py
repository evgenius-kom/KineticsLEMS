"""Model-based fitting engine — infrastructure scaffold.

This package lays the minimum scaffold for a future global
multi-experiment model-based fitter. Public surface:

* :class:`~.problem.FittingProblem` — bundles experiments + reaction
  model + objective into one immutable object.
* :class:`~.problem.FittingResult` — holds parameters, residuals,
  diagnostics; eventually round-trips through
  :class:`kinetics_lems.schemas.FitResultSchema`.
* :class:`~.objective.Objective` — pure-function loss components
  (currently NotImplemented but with the contract pinned).
* :class:`~.topology.Topology` — single / parallel / consecutive /
  competitive / mixed; all branches raise ``NotImplementedError``
  until the migration in TODO item 1 is started.

**Nothing in this package is wired into the runner yet.** Importing it
costs nothing, and the dataclasses give downstream code (reports,
schemas, prediction-engine) something to reference instead of
inventing types on the fly.

See ``docs/TODO_FEATURES.md`` item 1 ("Implement the model-based global
fitter") for the migration plan, effort estimate, and validation targets.
"""
from __future__ import annotations

from .objective import Objective, ObjectiveSpec, ObjectiveWeights
from .problem import FittingProblem, FittingResult, ParameterSpec
from .topology import Topology

__all__ = [
    "FittingProblem",
    "FittingResult",
    "Objective",
    "ObjectiveSpec",
    "ObjectiveWeights",
    "ParameterSpec",
    "Topology",
]

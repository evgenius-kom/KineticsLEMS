"""Model-based fitting engine — infrastructure scaffold.

KineticsLEMS today is **purely model-free**: Friedman, KAS, OFW, Vyazovkin,
Coats-Redfern + Z(α) master plot, then a separate A and lifetime
prediction. The next major direction (from
``kinetics_research_implementation_notes.md`` §2.3 / §6 P3) is a
**global multi-experiment model-based fitting engine**: fit a parametric
kinetic model (single step, parallel/consecutive/competitive multi-step,
or DAEM) against α(T) and dα/dt simultaneously across all heating rates.

This package lays the minimum scaffold needed to start filling that in
incrementally without breaking the existing API:

* :class:`~.problem.FittingProblem` — bundles experiments + reaction
  model + objective into one immutable object.
* :class:`~.problem.FittingResult` — holds parameters, residuals,
  diagnostics; eventually round-trips through
  :class:`kinetics_lems.schemas.FitResultSchema`.
* :class:`~.objective.Objective` — pure-function loss components
  (currently NotImplemented but with the contract pinned).
* :class:`~.topology.Topology` — single / parallel / consecutive /
  competitive / mixed; only ``SINGLE`` is implementable today, the
  others raise ``NotImplementedError``.

**Nothing in this package is wired into the runner yet.** Importing it
costs nothing, and the dataclasses give downstream code (reports,
schemas, prediction-engine) something to reference instead of
inventing types on the fly. The migration order:

1. Wire a *toy* SINGLE-step fitter (already feasible) into the runner
   under a feature flag.
2. Add PARALLEL once we have the synthetic 2-parallel dataset (#16).
3. Add CONSECUTIVE / COMPETITIVE.
4. Promote :class:`FittingResult` to schemas-backed serialisation.

References:
* ICTAC 2020 multi-step recommendations, DOI 10.1016/j.tca.2020.178597.
* Vyazovkin 2000, Thermochim. Acta 355, 155.
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

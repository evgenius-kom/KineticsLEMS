"""Smoke tests for the model-based fitting scaffold.

The scaffold is intentionally non-functional today; these tests pin down
the public surface (imports, dataclass shapes, NotImplementedError
contract) so the migration plan stays honest.
"""
from __future__ import annotations

import numpy as np
import pytest

from kinetics_lems.fitting import (
    FittingProblem,
    FittingResult,
    Objective,
    ObjectiveSpec,
    ObjectiveWeights,
    ParameterSpec,
    Topology,
)
from kinetics_lems.fitting.topology import build_ode_system


def test_parameter_spec_records_bounds_and_log_flag() -> None:
    p = ParameterSpec(
        name="logA", initial=10.0, bounds=(1.0, 25.0), unit="ln(1/s)", log_space=True
    )
    assert p.name == "logA"
    assert p.bounds == (1.0, 25.0)
    assert p.log_space is True


def test_topology_enum_has_five_branches() -> None:
    assert set(Topology) == {
        Topology.SINGLE,
        Topology.PARALLEL,
        Topology.CONSECUTIVE,
        Topology.COMPETITIVE,
        Topology.MIXED,
    }


def test_build_ode_system_raises_until_implemented() -> None:
    problem = FittingProblem(
        runs=[],
        topology=Topology.SINGLE,
        reaction_model="F1",
        parameters=[ParameterSpec(name="E", initial=120.0, bounds=(1.0, 600.0))],
        objective_spec=ObjectiveSpec(),
    )
    with pytest.raises(NotImplementedError):
        build_ode_system(problem, np.array([120.0]))


def test_objective_evaluate_raises_with_helpful_message() -> None:
    obj = Objective(ObjectiveSpec(weights=ObjectiveWeights(alpha=1.0, rate=0.5)))
    with pytest.raises(NotImplementedError) as exc:
        obj.evaluate(np.zeros(3), np.zeros(3), np.zeros(3), np.zeros(3))
    assert "objective.py" in str(exc.value)


def test_objective_normalised_sse_handles_zero_variance() -> None:
    val = Objective.normalised_sse(np.array([1.0, 2.0]), np.array([1.0, 1.0]))
    assert val == 1.0  # falls back to raw SSE when variance is zero


def test_fitting_result_starts_empty_and_accepts_assignments() -> None:
    problem = FittingProblem(
        runs=[],
        topology=Topology.SINGLE,
        reaction_model="F1",
        parameters=[],
        objective_spec=ObjectiveSpec(),
    )
    result = FittingResult(problem=problem)
    assert result.loss_total is None
    assert result.warnings == []
    result.warnings.append("toy warning")
    assert result.warnings == ["toy warning"]

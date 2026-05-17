"""Smoke tests for the canonical pydantic schemas.

The schemas are infrastructure-only today; we test that they instantiate,
validate, and round-trip through JSON. The full pipeline does not yet
consume them — that migration is tracked in schemas/__init__.py.
"""
from __future__ import annotations

import json

import pytest
from pydantic import ValidationError

from kinetics_lems.schemas import (
    ExperimentSchema,
    FitResultSchema,
    KineticModelSchema,
    KineticStepSchema,
    ParameterEstimate,
    SignalRole,
    SignalSchema,
    TechniqueEnum,
    TemperatureProgramMode,
    TemperatureProgramSchema,
)


def _experiment() -> ExperimentSchema:
    return ExperimentSchema(
        id="exp_001",
        technique=TechniqueEnum.DSC,
        program=TemperatureProgramSchema(
            mode=TemperatureProgramMode.DYNAMIC,
            heating_rate_K_min=10.0,
        ),
        signals=[
            SignalSchema(
                name="heat_flow",
                unit="mW",
                values=[0.1, 0.2, 0.3],
                role=SignalRole.HEAT_FLOW,
            )
        ],
    )


def test_experiment_minimum_required_fields() -> None:
    exp = _experiment()
    assert exp.id == "exp_001"
    assert exp.technique is TechniqueEnum.DSC
    assert exp.program.heating_rate_K_min == 10.0


def test_experiment_rejects_unknown_field() -> None:
    with pytest.raises(ValidationError):
        ExperimentSchema(
            id="x",
            technique=TechniqueEnum.DSC,
            program=TemperatureProgramSchema(),
            unknown_field=1,  # type: ignore[call-arg]
        )


def test_experiment_round_trips_through_json() -> None:
    exp = _experiment()
    dumped = exp.model_dump_json()
    restored = ExperimentSchema.model_validate_json(dumped)
    assert restored == exp


def test_fit_result_schema_with_parameters_and_diagnostics() -> None:
    fit = FitResultSchema(
        model_id="single_step_F1_v1",
        software_version="0.1.0",
        parameters=[
            ParameterEstimate(name="E", value=120.0, unit="kJ/mol", stderr=4.0),
            ParameterEstimate(name="logA", value=10.0, unit="ln(1/s)", stderr=0.3),
        ],
        diagnostics={"aic": 12.3, "bic": 14.1, "rss": 0.002, "warnings": []},
    )
    data = json.loads(fit.model_dump_json())
    assert data["parameters"][0]["name"] == "E"
    assert data["diagnostics"]["aic"] == 12.3


def test_kinetic_model_schema_parallel_topology_round_trip() -> None:
    model = KineticModelSchema(
        model_type="multi_step",
        topology="parallel",
        steps=[
            KineticStepSchema(
                id="step1",
                reaction_model="F1",
                parameters=[ParameterEstimate(name="E", value=120.0, unit="kJ/mol")],
                weight=0.6,
            ),
            KineticStepSchema(
                id="step2",
                reaction_model="F1",
                parameters=[ParameterEstimate(name="E", value=180.0, unit="kJ/mol")],
                weight=0.4,
            ),
        ],
    )
    restored = KineticModelSchema.model_validate_json(model.model_dump_json())
    assert restored.steps[0].weight == 0.6
    assert restored.steps[1].parameters[0].value == 180.0

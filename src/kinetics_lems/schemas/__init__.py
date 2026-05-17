"""Canonical pydantic schemas for inputs, models, fits, and predictions.

Status: **infrastructure only**. The shapes below are intentionally close to
what the existing dataclasses (``models.py``, ``conversion.py``,
``methods/*``) carry around, but currently *nothing in the main pipeline
consumes them yet* — they exist so that:

1. Vendor adapters (:mod:`kinetics_lems.io.vendors`) have a canonical
   destination type once they grow beyond the legacy ``CaseData`` shape.
2. Fit results (:mod:`kinetics_lems.fitting`) can be serialised
   to/from JSON without bespoke per-method code.
3. Reports can render any pipeline artefact via a uniform interface.

Why not pydantic everywhere right now? The frozen dataclasses are fine
for hot per-α loops (no validation overhead). Pydantic earns its keep at
I/O boundaries (load case → validate, write fit result → schema-stable
JSON). That's the boundary we are setting up.

See ``docs/TODO_FEATURES.md`` item 5 ("Promote schemas out of stub
status") for the migration plan.
"""
from __future__ import annotations

from .canonical import (
    AtmosphereSpec,
    ConversionDefinition,
    ExperimentSchema,
    FitResultSchema,
    KineticModelSchema,
    KineticStepSchema,
    ParameterEstimate,
    PredictionResultSchema,
    PreprocessingSpec,
    SampleSpec,
    SignalRole,
    SignalSchema,
    TechniqueEnum,
    TemperatureProgramMode,
    TemperatureProgramSchema,
)

__all__ = [
    "AtmosphereSpec",
    "ConversionDefinition",
    "ExperimentSchema",
    "FitResultSchema",
    "KineticModelSchema",
    "KineticStepSchema",
    "ParameterEstimate",
    "PredictionResultSchema",
    "PreprocessingSpec",
    "SampleSpec",
    "SignalRole",
    "SignalSchema",
    "TechniqueEnum",
    "TemperatureProgramMode",
    "TemperatureProgramSchema",
]

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

When you promote any of these schemas out of stub status, the migration is:

* replace the corresponding ``@dataclass`` in ``models.py`` /
  ``methods/common.py`` with the pydantic ``BaseModel``;
* add a ``.to_legacy()`` adapter if call sites still expect the dataclass;
* update tests last — type narrowing in pydantic catches more upstream.

Why not pydantic everywhere right now? The frozen dataclasses are fine
for hot per-α loops (no validation overhead) and Type-42 plot generators.
Pydantic earns its keep at I/O boundaries (load case → validate, write
fit result → schema-stable JSON). That's the boundary we are setting up.

References:
    ICTAC 2020 recommendations §1 (data documentation requirements):
    Thermochim. Acta 689 (2020) 178597.
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

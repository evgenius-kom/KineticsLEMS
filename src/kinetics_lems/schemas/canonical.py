"""Canonical pydantic models for experiments, kinetic models, and fit results.

Most fields are ``Optional`` because vendor adapters fill in what they have;
the *required* core is enough to run the existing isoconversional pipeline:

* :class:`ExperimentSchema` — ``id``, ``technique``, one ``signals`` entry,
  and a ``program`` describing T(t).
* :class:`KineticModelSchema` — at least one :class:`KineticStepSchema` with
  ``E`` and ``logA``.
* :class:`FitResultSchema` — ``model_id`` + ``parameters`` list, the
  minimum for a JSON sidecar next to the existing CSV outputs.

These intentionally mirror §5.2/5.3 of ``kinetics_research_implementation_notes.md``.
"""
from __future__ import annotations

from enum import StrEnum
from typing import Any

from pydantic import BaseModel, ConfigDict, Field

# --- Enums ------------------------------------------------------------------


class TechniqueEnum(StrEnum):
    """Recognised thermal-analysis techniques.

    The current pipeline treats DSC/TGA/FSC/POM identically (scalar y vs T);
    this enum is for documentation + future technique-specific logic
    (TGA mass normalisation, MS quantitation, ...).
    """

    DSC = "DSC"
    TGA = "TGA"
    FSC = "FSC"
    POM = "POM"
    STA = "STA"
    DTA = "DTA"
    DIL = "DIL"  # dilatometry
    DEA = "DEA"  # dielectric analysis
    ARC = "ARC"  # accelerating rate calorimetry
    MS = "MS"  # mass spectrometry channel
    FTIR = "FTIR"  # infrared band
    TMA = "TMA"
    CUSTOM = "CUSTOM"


class TemperatureProgramMode(StrEnum):
    DYNAMIC = "dynamic"  # constant β linear ramp
    ISOTHERMAL = "isothermal"
    STEP = "step"
    ARBITRARY = "arbitrary"  # user-supplied T(t)


class SignalRole(StrEnum):
    MASS = "mass"
    MASS_LOSS = "mass_loss"
    HEAT_FLOW = "heat_flow"
    DSC = "dsc"
    DTG = "dtg"
    CONVERSION = "conversion"
    RATE = "rate"
    MS_CHANNEL = "ms_channel"
    FTIR_BAND = "ftir_band"
    CUSTOM = "custom"


# --- Sub-schemas ------------------------------------------------------------


class SampleSpec(BaseModel):
    name: str = "unknown"
    mass_initial_mg: float | None = None
    mass_final_mg: float | None = None
    geometry: str | None = None
    particle_size: str | None = None


class AtmosphereSpec(BaseModel):
    gas: str | None = None
    flow_ml_min: float | None = None
    pressure_kPa: float | None = None


class TemperatureProgramSchema(BaseModel):
    mode: TemperatureProgramMode = TemperatureProgramMode.DYNAMIC
    heating_rate_K_min: float | None = None
    """β for ``DYNAMIC`` mode; ``None`` for non-linear T(t)."""

    time_s: list[float] | None = None
    temperature_K: list[float] | None = None
    """Required for ``ARBITRARY`` / ``STEP`` modes; optional otherwise."""

    isothermal_T_K: float | None = None


class SignalSchema(BaseModel):
    name: str
    unit: str
    values: list[float]
    role: SignalRole = SignalRole.CUSTOM


class ConversionDefinition(BaseModel):
    """How α(T) was computed from the raw signal."""

    method: str  # e.g. "tga_mass_loss", "dsc_area"
    baseline: str | None = None  # "linear" today; "sigmoid" / "spline" TODO
    smoothing: dict[str, Any] | None = None  # SavGol parameters if applied


class PreprocessingSpec(BaseModel):
    """Provenance for everything between raw signal and α(T)."""

    baseline: dict[str, Any] | None = None
    smoothing: dict[str, Any] | None = None
    conversion: ConversionDefinition | None = None


class ExperimentSchema(BaseModel):
    """One run at one heating-rate program."""

    model_config = ConfigDict(extra="forbid")

    id: str
    source_file: str | None = None
    source_vendor: str | None = None
    instrument: str | None = None
    technique: TechniqueEnum
    sample: SampleSpec = Field(default_factory=SampleSpec)
    atmosphere: AtmosphereSpec = Field(default_factory=AtmosphereSpec)
    program: TemperatureProgramSchema
    signals: list[SignalSchema] = Field(default_factory=list)
    preprocessing: PreprocessingSpec | None = None

    # TODO: round-trip adapters to/from legacy CaseData live in
    # schemas/adapters.py (not yet implemented). Open a new file rather than
    # cluttering this one; see schemas/__init__.py docstring for the
    # migration plan.


# --- Kinetic models ---------------------------------------------------------


class ParameterEstimate(BaseModel):
    """One fitted parameter with uncertainty and bounds."""

    name: str
    value: float
    unit: str
    stderr: float | None = None
    bounds: tuple[float, float] | None = None


class KineticStepSchema(BaseModel):
    """One step of a (possibly multi-step) kinetic model."""

    id: str
    reaction_model: str  # "F1", "A2", "SB", "Kamal", ...
    parameters: list[ParameterEstimate]
    weight: float | None = None  # for parallel-step models; sum→1


class KineticModelSchema(BaseModel):
    """Multi-step topology + parameters."""

    model_config = ConfigDict(extra="forbid")

    model_type: str  # "single_step" | "multi_step" | "DAEM" | ...
    topology: str | None = None  # "parallel" | "consecutive" | "competitive" | "mixed"
    steps: list[KineticStepSchema]
    validity: dict[str, Any] | None = None
    """Temperature/heating-rate range, atmosphere — extrapolation guard."""


# --- Fit / prediction outputs ----------------------------------------------


class FitResultSchema(BaseModel):
    """JSON-stable sidecar for whatever fit produced ``model_id``."""

    model_config = ConfigDict(extra="forbid")

    model_id: str
    software_version: str
    data_hashes: list[str] = Field(default_factory=list)
    fit_scope: dict[str, Any] | None = None
    """``{"experiments_used": [...], "holdout_experiments": [...]}``."""

    parameters: list[ParameterEstimate]
    parameter_correlation: list[list[float]] | None = None
    objective: dict[str, float] | None = None  # loss_total, loss_alpha, loss_rate
    diagnostics: dict[str, Any] | None = None  # aic, bic, rss, warnings


class PredictionResultSchema(BaseModel):
    """Output of a simulation under some T(t) profile."""

    model_config = ConfigDict(extra="forbid")

    model_id: str
    program: TemperatureProgramSchema
    time_s: list[float]
    temperature_K: list[float]
    alpha: list[float]
    dalpha_dt: list[float] | None = None
    extrapolation_warnings: list[str] = Field(default_factory=list)

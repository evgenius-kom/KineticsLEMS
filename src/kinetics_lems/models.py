"""Domain types for an experimental case and per-heating-rate runs."""
from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field
from enum import StrEnum

import numpy as np


class ExperimentType(StrEnum):
    ISOTHERMAL = "isothermal"
    HEATING = "heating"
    COOLING = "cooling"


class Method(StrEnum):
    DSC = "DSC"
    TGA = "TGA"
    FSC = "FSC"
    POM = "POM"


@dataclass(frozen=True)
class Wave:
    """A 2D experimental curve (e.g. T vs heat-flow, or t vs mass)."""

    x: np.ndarray
    y: np.ndarray

    def __post_init__(self) -> None:
        if self.x.shape != self.y.shape:
            raise ValueError(f"Wave x/y shape mismatch: {self.x.shape} vs {self.y.shape}")
        if self.x.ndim != 1:
            raise ValueError("Wave is 1D only")
        if self.x.size < 2:
            raise ValueError("Wave needs at least two points")

    def __len__(self) -> int:
        return self.x.size

    def linear_baseline(self) -> Wave:
        """Linear baseline drawn between first and last points."""
        bl_y = np.interp(self.x, [self.x[0], self.x[-1]], [self.y[0], self.y[-1]])
        return Wave(self.x.copy(), bl_y)

    def subtract_baseline(self) -> Wave:
        return Wave(self.x.copy(), self.y - self.linear_baseline().y)


@dataclass(frozen=True)
class CaseParams:
    material: str
    experiment_type: ExperimentType
    method: Method
    file_to_condition: Mapping[str, float]
    """For heating/cooling: heating rate β (K/min). For isothermal: temperature (K)."""


@dataclass
class CaseData:
    params: CaseParams
    waves: dict[float, Wave] = field(default_factory=dict)
    """key = heating rate (K/min) or isothermal temperature (K), value = wave."""

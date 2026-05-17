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
    """A 2D experimental curve.

    Two storage modes:

    * **2-column** (``x`` = T or t, ``y`` = signal, ``t_seconds = None``)
      — the legacy linear-heating case. Math in :mod:`kinetics_lems.conversion`
      assumes ``dT/dt = β = const``.

    * **3-column** (``x`` = T, ``y`` = signal, ``t_seconds`` = time array
      with the same shape as ``x``) — arbitrary T(t) program. The
      conversion engine uses the recorded ``t`` to compute α(t) and
      dα/dt without any linear-heating assumption.

    ``t_seconds`` is optional and additive — older code reading ``wave.x``
    / ``wave.y`` keeps working unchanged.
    """

    x: np.ndarray
    y: np.ndarray
    t_seconds: np.ndarray | None = None
    """Recorded time in seconds (3-column files); None for 2-column."""

    def __post_init__(self) -> None:
        if self.x.shape != self.y.shape:
            raise ValueError(f"Wave x/y shape mismatch: {self.x.shape} vs {self.y.shape}")
        if self.x.ndim != 1:
            raise ValueError("Wave is 1D only")
        if self.x.size < 2:
            raise ValueError("Wave needs at least two points")
        if self.t_seconds is not None:
            if self.t_seconds.shape != self.x.shape:
                raise ValueError(
                    f"Wave t shape mismatch: {self.t_seconds.shape} vs {self.x.shape}"
                )
            if not np.all(np.diff(self.t_seconds) > 0):
                raise ValueError("Wave t must be strictly monotone increasing")

    def __len__(self) -> int:
        return self.x.size

    @property
    def has_recorded_time(self) -> bool:
        return self.t_seconds is not None

    def linear_baseline(self) -> Wave:
        """Linear baseline drawn between first and last points."""
        bl_y = np.interp(self.x, [self.x[0], self.x[-1]], [self.y[0], self.y[-1]])
        return Wave(self.x.copy(), bl_y, self.t_seconds)

    def subtract_baseline(self) -> Wave:
        return Wave(
            self.x.copy(),
            self.y - self.linear_baseline().y,
            self.t_seconds,
        )


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

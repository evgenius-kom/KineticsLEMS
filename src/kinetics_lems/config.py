"""Load and validate the algorithm configuration from a TOML file."""
from __future__ import annotations

import tomllib
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

DEFAULT_CONFIG_PATH = Path(__file__).resolve().parent.parent.parent / "configs" / "default.toml"

ALL_METHODS = ("friedman", "kas", "ofw", "kissinger", "vyazovkin", "vyazovkin_aic")


@dataclass(frozen=True)
class ConversionConfig:
    step: float = 0.05
    min: float = 0.05
    max: float = 0.95

    def grid(self) -> np.ndarray:
        """Return an evenly-spaced α grid on [min, max].

        ``step`` controls the *count* of points: ``n ≈ (max − min) / step + 1``,
        rounded to the nearest integer. The actual spacing
        ``(max − min) / (n − 1)`` may differ slightly when ``(max − min)``
        is not an integer multiple of ``step``.
        """
        if not 0.0 < self.min < self.max < 1.0:
            raise ValueError(
                f"Bad α range: need 0 < min < max < 1, got min={self.min}, max={self.max}"
            )
        if self.step <= 0:
            raise ValueError(f"Bad α step: {self.step}")
        n = int(round((self.max - self.min) / self.step)) + 1
        return np.linspace(self.min, self.max, n)


@dataclass(frozen=True)
class VyazovkinConfig:
    ea_bracket_kJ: tuple[float, float] = (1.0, 600.0)


@dataclass(frozen=True)
class VyazovkinAICConfig:
    delta_alpha: float = 0.02
    ea_bracket_kJ: tuple[float, float] = (1.0, 600.0)


@dataclass(frozen=True)
class OutputConfig:
    directory: str = "out"
    save_plots: bool = True
    save_csv: bool = True
    plot_dpi: int = 120


@dataclass(frozen=True)
class Config:
    conversion: ConversionConfig = field(default_factory=ConversionConfig)
    enabled_methods: tuple[str, ...] = ALL_METHODS
    vyazovkin: VyazovkinConfig = field(default_factory=VyazovkinConfig)
    vyazovkin_aic: VyazovkinAICConfig = field(default_factory=VyazovkinAICConfig)
    output: OutputConfig = field(default_factory=OutputConfig)


def load_config(path: Path | str | None = None) -> Config:
    if path is None:
        path = DEFAULT_CONFIG_PATH
    path = Path(path)
    with path.open("rb") as f:
        raw = tomllib.load(f)
    return _parse(raw)


def _parse(raw: dict) -> Config:
    conv_raw = raw.get("conversion", {})
    conversion = ConversionConfig(
        step=float(conv_raw.get("step", 0.05)),
        min=float(conv_raw.get("min", 0.05)),
        max=float(conv_raw.get("max", 0.95)),
    )

    methods_raw = raw.get("methods", {})
    enabled = tuple(methods_raw.get("enabled", ALL_METHODS))
    for m in enabled:
        if m not in ALL_METHODS:
            raise ValueError(f"Unknown method '{m}'. Allowed: {ALL_METHODS}")

    vya_raw = methods_raw.get("vyazovkin", {})
    vya = VyazovkinConfig(ea_bracket_kJ=_pair(vya_raw.get("ea_bracket_kJ", (1.0, 600.0))))

    aic_raw = methods_raw.get("vyazovkin_aic", {})
    aic = VyazovkinAICConfig(
        delta_alpha=float(aic_raw.get("delta_alpha", 0.02)),
        ea_bracket_kJ=_pair(aic_raw.get("ea_bracket_kJ", (1.0, 600.0))),
    )

    out_raw = raw.get("output", {})
    output = OutputConfig(
        directory=str(out_raw.get("directory", "out")),
        save_plots=bool(out_raw.get("save_plots", True)),
        save_csv=bool(out_raw.get("save_csv", True)),
        plot_dpi=int(out_raw.get("plot_dpi", 120)),
    )

    return Config(
        conversion=conversion,
        enabled_methods=enabled,
        vyazovkin=vya,
        vyazovkin_aic=aic,
        output=output,
    )


def _pair(value) -> tuple[float, float]:
    a, b = value
    return (float(a), float(b))


__all__ = [
    "Config",
    "ConversionConfig",
    "VyazovkinConfig",
    "VyazovkinAICConfig",
    "OutputConfig",
    "ALL_METHODS",
    "DEFAULT_CONFIG_PATH",
    "load_config",
]

"""Load and validate the algorithm configuration from a TOML file."""
from __future__ import annotations

import tomllib
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

DEFAULT_CONFIG_PATH = Path(__file__).resolve().parent.parent.parent / "configs" / "default.toml"

ALL_METHODS = (
    "friedman",
    "kas",
    "ofw",
    "kissinger",
    "vyazovkin",
    "vyazovkin_aic",
    "master_plot",
    "preexponential",
    "multistep",
)


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
class PreexponentialConfig:
    model: str | None = None
    """f(α) model name (e.g. ``"F1"``). ``None`` → take best model from
    Z(α) master-plot ranking; falls back to ``"F1"`` if master_plot disabled."""


@dataclass(frozen=True)
class MultiStepConfig:
    jump_threshold: float = 0.10
    """Relative |ΔE / median(E)| above which a new segment starts."""
    min_segment_size: int = 3
    """Smallest segment in α-points; smaller segments are merged."""


@dataclass(frozen=True)
class OutputConfig:
    directory: str = "out"
    save_plots: bool = True
    save_csv: bool = True
    plot_dpi: int = 300
    plot_formats: tuple[str, ...] = ("png", "pdf")
    """Vector formats first (pdf/svg) for journals; png for screen previews."""
    per_method_panels: bool = False
    """If True, also write one figure per method in addition to the overlay."""


@dataclass(frozen=True)
class Config:
    conversion: ConversionConfig = field(default_factory=ConversionConfig)
    enabled_methods: tuple[str, ...] = ALL_METHODS
    vyazovkin: VyazovkinConfig = field(default_factory=VyazovkinConfig)
    vyazovkin_aic: VyazovkinAICConfig = field(default_factory=VyazovkinAICConfig)
    preexponential: PreexponentialConfig = field(default_factory=PreexponentialConfig)
    multistep: MultiStepConfig = field(default_factory=MultiStepConfig)
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

    preexp_raw = methods_raw.get("preexponential", {})
    raw_model = preexp_raw.get("model")
    preexp = PreexponentialConfig(
        # Empty string from TOML → auto-pick from master_plot.
        model=raw_model if raw_model else None,
    )

    multistep_raw = methods_raw.get("multistep", {})
    multistep = MultiStepConfig(
        jump_threshold=float(multistep_raw.get("jump_threshold", 0.10)),
        min_segment_size=int(multistep_raw.get("min_segment_size", 3)),
    )

    out_raw = raw.get("output", {})
    output = OutputConfig(
        directory=str(out_raw.get("directory", "out")),
        save_plots=bool(out_raw.get("save_plots", True)),
        save_csv=bool(out_raw.get("save_csv", True)),
        plot_dpi=int(out_raw.get("plot_dpi", 300)),
        plot_formats=tuple(out_raw.get("plot_formats", ("png", "pdf"))),
        per_method_panels=bool(out_raw.get("per_method_panels", False)),
    )

    return Config(
        conversion=conversion,
        enabled_methods=enabled,
        vyazovkin=vya,
        vyazovkin_aic=aic,
        preexponential=preexp,
        multistep=multistep,
        output=output,
    )


def _pair(value) -> tuple[float, float]:
    a, b = value
    return (float(a), float(b))


__all__ = [
    "ALL_METHODS",
    "Config",
    "ConversionConfig",
    "DEFAULT_CONFIG_PATH",
    "MultiStepConfig",
    "OutputConfig",
    "PreexponentialConfig",
    "VyazovkinAICConfig",
    "VyazovkinConfig",
    "load_config",
]

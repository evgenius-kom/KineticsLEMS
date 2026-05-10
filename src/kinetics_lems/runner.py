"""Run the configured isoconversional analysis on a CaseData."""
from __future__ import annotations

from dataclasses import dataclass

from .config import Config
from .conversion import build_runs
from .methods import (
    IsoconversionalResult,
    KissingerResult,
    friedman,
    kas,
    kissinger,
    ofw,
    vyazovkin,
    vyazovkin_aic,
)
from .models import CaseData


@dataclass
class AnalysisResults:
    isoconversional: dict[str, IsoconversionalResult]  # method name -> result
    kissinger: KissingerResult | None = None


def run_analysis(case: CaseData, config: Config) -> AnalysisResults:
    runs = build_runs(case)
    if len(runs) < 2:
        raise ValueError(
            f"At least 2 heating rates are required for isoconversional analysis "
            f"(got {len(runs)})"
        )

    alphas = config.conversion.grid()

    iso: dict[str, IsoconversionalResult] = {}
    kiss: KissingerResult | None = None

    for name in config.enabled_methods:
        if name == "friedman":
            iso[name] = friedman(runs, alphas)
        elif name == "kas":
            iso[name] = kas(runs, alphas)
        elif name == "ofw":
            iso[name] = ofw(runs, alphas)
        elif name == "kissinger":
            kiss = kissinger(runs)
        elif name == "vyazovkin":
            iso[name] = vyazovkin(
                runs, alphas, Ea_bracket_kJ=config.vyazovkin.ea_bracket_kJ
            )
        elif name == "vyazovkin_aic":
            iso[name] = vyazovkin_aic(
                runs,
                alphas,
                delta_alpha=config.vyazovkin_aic.delta_alpha,
                Ea_bracket_kJ=config.vyazovkin_aic.ea_bracket_kJ,
            )
        else:  # pragma: no cover — guarded by config validation
            raise ValueError(f"Unknown method '{name}'")

    return AnalysisResults(isoconversional=iso, kissinger=kiss)


__all__ = ["AnalysisResults", "run_analysis"]

"""Run the configured isoconversional analysis on a CaseData."""
from __future__ import annotations

from dataclasses import dataclass

from .config import Config
from .conversion import build_runs
from .methods import (
    IsoconversionalResult,
    KissingerResult,
    ModelRanking,
    MultiStepResult,
    PreexponentialResult,
    compute_A,
    detect_steps,
    friedman,
    kas,
    kissinger,
    ofw,
    rank_models,
    vyazovkin,
    vyazovkin_aic,
)
from .models import CaseData


@dataclass
class AnalysisResults:
    isoconversional: dict[str, IsoconversionalResult]  # method name -> result
    kissinger: KissingerResult | None = None
    model_ranking: ModelRanking | None = None  # Z(α) master-plot ranking
    preexponential: PreexponentialResult | None = None  # A from kinetic triplet
    multistep: MultiStepResult | None = None  # E(α) segmentation


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
    ranking: ModelRanking | None = None
    preexp: PreexponentialResult | None = None
    multistep: MultiStepResult | None = None

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
        elif name == "master_plot":
            ranking = rank_models(runs, alphas)
        elif name == "preexponential":
            # Needs both Vyazovkin (E_a per α) and master_plot (best f(α)).
            ea_source = iso.get("vyazovkin") or iso.get("vyazovkin_aic") or iso.get("kas")
            if ea_source is None:
                raise ValueError(
                    "preexponential needs vyazovkin / vyazovkin_aic / kas to run first"
                )
            chosen = (
                config.preexponential.model
                if config.preexponential.model is not None
                else (ranking.best_model if ranking is not None else "F1")
            )
            preexp = compute_A(
                runs,
                ea_source.alpha,
                ea_source.Ea_J_per_mol,
                model=chosen,
            )
        elif name == "multistep":
            ea_source = iso.get("vyazovkin") or iso.get("vyazovkin_aic") or iso.get("kas")
            if ea_source is None:
                raise ValueError(
                    "multistep needs vyazovkin / vyazovkin_aic / kas to run first"
                )
            multistep = detect_steps(
                ea_source,
                relative_jump_threshold=config.multistep.jump_threshold,
                min_segment_size=config.multistep.min_segment_size,
            )
        else:  # pragma: no cover — guarded by config validation
            raise ValueError(f"Unknown method '{name}'")

    return AnalysisResults(
        isoconversional=iso,
        kissinger=kiss,
        model_ranking=ranking,
        preexponential=preexp,
        multistep=multistep,
    )


__all__ = ["AnalysisResults", "run_analysis"]

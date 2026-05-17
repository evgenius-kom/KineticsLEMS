"""Run the configured isoconversional analysis on a CaseData."""
from __future__ import annotations

from dataclasses import dataclass
from functools import partial

from .config import Config
from .conversion import build_runs
from .methods import (
    CoatsRedfernResult,
    IsoconversionalResult,
    KissingerResult,
    LifetimeSummary,
    ModelRanking,
    MultiStepResult,
    PreexponentialResult,
    ReactionOrderResult,
    UncertaintyResult,
    coats_redfern,
    compute_A,
    detect_steps,
    friedman,
    jackknife_isoconversional,
    kas,
    kissinger,
    ofw,
    predict_at_temperatures,
    rank_models,
    reaction_order,
    vyazovkin,
    vyazovkin_aic,
)
from .models import CaseData

_ISO_ESTIMATORS = {
    "vyazovkin": vyazovkin,
    "vyazovkin_aic": vyazovkin_aic,
    "kas": kas,
    "ofw": ofw,
    "friedman": friedman,
}


@dataclass
class AnalysisResults:
    isoconversional: dict[str, IsoconversionalResult]  # method name -> result
    kissinger: KissingerResult | None = None
    model_ranking: ModelRanking | None = None
    preexponential: PreexponentialResult | None = None
    multistep: MultiStepResult | None = None
    reaction_order: ReactionOrderResult | None = None
    coats_redfern: CoatsRedfernResult | None = None
    uncertainty: UncertaintyResult | None = None
    lifetime: LifetimeSummary | None = None


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
    ro: ReactionOrderResult | None = None
    cr: CoatsRedfernResult | None = None
    unc: UncertaintyResult | None = None
    lt: LifetimeSummary | None = None

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
            ea_source = _pick_ea_source(iso)
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
            ea_source = _pick_ea_source(iso)
            multistep = detect_steps(
                ea_source,
                relative_jump_threshold=config.multistep.jump_threshold,
                min_segment_size=config.multistep.min_segment_size,
            )
        elif name == "reaction_order":
            ro = reaction_order(
                runs,
                alpha_min=config.reaction_order.alpha_min,
                alpha_max=config.reaction_order.alpha_max,
                n_min=config.reaction_order.n_min,
                n_max=config.reaction_order.n_max,
                n_steps=config.reaction_order.n_steps,
            )
        elif name == "coats_redfern":
            cr = coats_redfern(
                runs,
                alpha_min=config.coats_redfern.alpha_min,
                alpha_max=config.coats_redfern.alpha_max,
                models=list(config.coats_redfern.models)
                if config.coats_redfern.models is not None
                else None,
            )
        elif name == "uncertainty":
            target = config.uncertainty.method
            if target not in _ISO_ESTIMATORS:
                raise ValueError(
                    f"uncertainty.method must be one of {sorted(_ISO_ESTIMATORS)}, "
                    f"got {target!r}"
                )
            estimator = _ISO_ESTIMATORS[target]
            if target == "vyazovkin":
                estimator = partial(
                    vyazovkin, Ea_bracket_kJ=config.vyazovkin.ea_bracket_kJ
                )
            elif target == "vyazovkin_aic":
                estimator = partial(
                    vyazovkin_aic,
                    delta_alpha=config.vyazovkin_aic.delta_alpha,
                    Ea_bracket_kJ=config.vyazovkin_aic.ea_bracket_kJ,
                )
            unc = jackknife_isoconversional(
                runs, alphas, estimator, method_name=target
            )
        elif name == "lifetime":
            if preexp is None:
                raise ValueError(
                    "lifetime needs preexponential to run first (provides A)"
                )
            ea_source = _pick_ea_source(iso)
            model = (
                config.lifetime.model
                if config.lifetime.model is not None
                else (ranking.best_model if ranking is not None else "F1")
            )
            T_K_list = [t + 273.15 for t in config.lifetime.temperatures_C]
            lt = predict_at_temperatures(
                T_K_list=T_K_list,
                Ea_J_per_mol=ea_source.Ea_J_per_mol,
                alpha_grid=ea_source.alpha,
                A_per_sec=preexp.A_per_sec_median,
                model=model,
                alpha_targets=config.lifetime.alpha_targets,
            )
        else:  # pragma: no cover — guarded by config validation
            raise ValueError(f"Unknown method '{name}'")

    return AnalysisResults(
        isoconversional=iso,
        kissinger=kiss,
        model_ranking=ranking,
        preexponential=preexp,
        multistep=multistep,
        reaction_order=ro,
        coats_redfern=cr,
        uncertainty=unc,
        lifetime=lt,
    )


def _pick_ea_source(iso: dict[str, IsoconversionalResult]) -> IsoconversionalResult:
    """Pick the cleanest available E(α) source for downstream methods."""
    for name in ("vyazovkin", "vyazovkin_aic", "kas"):
        if name in iso:
            return iso[name]
    raise ValueError(
        "Downstream method needs one of vyazovkin / vyazovkin_aic / kas to run first"
    )


__all__ = ["AnalysisResults", "run_analysis"]

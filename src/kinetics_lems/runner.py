"""Run the configured isoconversional analysis on a CaseData."""
from __future__ import annotations

from dataclasses import dataclass, field
from functools import partial

from .config import Config
from .conversion import build_runs
from .methods import (
    CoatsRedfernResult,
    CompensationFit,
    ConsistencyResult,
    EmpiricalFit,
    EndpointReliability,
    IsoconversionalResult,
    KissingerResult,
    LifetimeSummary,
    ModelRanking,
    MultiStepResult,
    PreexponentialResult,
    ReactionOrderResult,
    UncertaintyResult,
    assess_endpoints,
    coats_redfern,
    compensation_from_coats_redfern,
    compute_A,
    consistency_check,
    detect_steps,
    fit_prout_tompkins,
    fit_sestak_berggren,
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
    consistency: ConsistencyResult | None = None
    compensation: CompensationFit | None = None
    empirical_fits: dict[str, EmpiricalFit] = field(default_factory=dict)
    """Optional Prout-Tompkins / Sestak-Berggren fits; enabled via config."""
    lifetime: LifetimeSummary | None = None
    endpoint_reliability: dict[str, EndpointReliability] = field(default_factory=dict)
    """Per-method endpoint warnings — always populated when at least one
    isoconversional curve is present. Empty dict ↔ no methods enabled."""


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
    cons: ConsistencyResult | None = None
    comp: CompensationFit | None = None
    lt: LifetimeSummary | None = None

    for name in config.enabled_methods:
        if name == "friedman":
            iso[name] = friedman(
                runs,
                alphas,
                smooth_window=config.friedman.smooth_window,
                smooth_poly=config.friedman.smooth_poly,
            )
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
        elif name == "consistency":
            if len(iso) < 2:
                # Need at least two isoconversional curves to compare.
                cons = None
            else:
                cons = consistency_check(iso, threshold=config.consistency.threshold)
        elif name == "compensation":
            # Source is Coats-Redfern by default — gives many (E, A) pairs
            # across all 12 models × N heating rates, which makes the line
            # tight enough to interpret. The isoconversional path is in
            # methods.compensation but needs preexponential to have run.
            if cr is not None and len(cr.fits) >= 2:
                comp = compensation_from_coats_redfern(cr)
            else:
                comp = None
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

    # Optional empirical-model fits (Prout-Tompkins / Sestak-Berggren).
    # Off by default — these are non-identifiable and should only be used
    # when the 12 canonical models give a poor Z(α) match.
    empirical_fits: dict[str, EmpiricalFit] = {}
    if config.empirical_models.enable_prout_tompkins:
        empirical_fits["prout_tompkins"] = fit_prout_tompkins(runs, alphas)
    if config.empirical_models.enable_sestak_berggren:
        empirical_fits["sestak_berggren"] = fit_sestak_berggren(runs, alphas)

    # Endpoint reliability is cheap — compute it for every iso method that ran.
    endpoint_reliability: dict[str, EndpointReliability] = {
        name: assess_endpoints(
            res,
            alpha_low=config.conversion.min,
            alpha_high=config.conversion.max,
        )
        for name, res in iso.items()
    }

    return AnalysisResults(
        isoconversional=iso,
        kissinger=kiss,
        model_ranking=ranking,
        preexponential=preexp,
        multistep=multistep,
        reaction_order=ro,
        coats_redfern=cr,
        uncertainty=unc,
        consistency=cons,
        compensation=comp,
        empirical_fits=empirical_fits,
        lifetime=lt,
        endpoint_reliability=endpoint_reliability,
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

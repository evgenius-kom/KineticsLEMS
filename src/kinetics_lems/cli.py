"""Command-line entry points for KineticsLEMS."""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

from .config import DEFAULT_CONFIG_PATH, load_config
from .io import load_case
from .reporting import plot_ea_per_method, plot_ea_vs_alpha, plot_kissinger, write_csv
from .reporting_coats_redfern import plot_coats_redfern, write_coats_redfern_csv
from .reporting_compensation import plot_compensation, write_compensation_csv
from .reporting_consistency import write_consistency_csv
from .reporting_lifetime import plot_lifetime, write_lifetime_csv
from .reporting_markdown import write_markdown_report
from .reporting_master_plot import plot_master_plot, write_master_plot_csv
from .reporting_multistep import plot_multistep, write_multistep_csv
from .reporting_preexp import plot_preexp, write_preexp_csv
from .reporting_reaction_order import plot_reaction_order, write_reaction_order_csv
from .reporting_uncertainty import plot_uncertainty, write_uncertainty_csv
from .runner import run_analysis
from .synthetic import generate_case, write_case


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="kinetics-lems",
        description="Isoconversional kinetic analysis of thermal-analysis data.",
    )
    sub = parser.add_subparsers(dest="cmd", required=True)

    p_run = sub.add_parser("analyze", help="Run analysis on a case (folder or .zip)")
    p_run.add_argument("case", type=Path, help="Path to case folder or .zip archive")
    p_run.add_argument(
        "--config",
        type=Path,
        default=DEFAULT_CONFIG_PATH,
        help="Path to algorithm config TOML (default: bundled configs/default.toml)",
    )
    p_run.add_argument(
        "--out", type=Path, default=None, help="Output directory (overrides config)"
    )

    p_gen = sub.add_parser("generate", help="Generate a synthetic case for testing/validation")
    p_gen.add_argument("output", type=Path, help="Output directory for the synthetic case")
    p_gen.add_argument(
        "--rates", type=float, nargs="+", default=[2.5, 5.0, 10.0, 20.0], help="Heating rates K/min"
    )
    p_gen.add_argument("--ea", type=float, default=120.0, help="True E_a (kJ/mol)")
    p_gen.add_argument("--A", type=float, default=1.0e10, help="True pre-exponential A (1/s)")
    p_gen.add_argument(
        "--model", type=str, default="F1", help="Reaction model: F1, F2, A2, R2, R3"
    )
    p_gen.add_argument("--noise", type=float, default=0.0, help="Relative noise stddev")
    p_gen.add_argument("--seed", type=int, default=0)
    p_gen.add_argument(
        "--T-start", type=float, default=None, help="Lower T (K); auto if omitted"
    )
    p_gen.add_argument(
        "--T-stop", type=float, default=None, help="Upper T (K); auto if omitted"
    )

    return parser


def main(argv: list[str] | None = None) -> int:
    args = _build_parser().parse_args(argv)
    if args.cmd == "analyze":
        return _cmd_analyze(args)
    if args.cmd == "generate":
        return _cmd_generate(args)
    return 2  # pragma: no cover


def _cmd_analyze(args: argparse.Namespace) -> int:
    config = load_config(args.config)
    out_dir = args.out or Path(config.output.directory)

    case = load_case(args.case)
    results = run_analysis(case, config)

    if config.output.save_csv:
        write_csv(results, out_dir)
        if results.model_ranking is not None:
            write_master_plot_csv(results.model_ranking, out_dir)
        if results.preexponential is not None:
            write_preexp_csv(results.preexponential, out_dir)
        if results.multistep is not None:
            write_multistep_csv(results.multistep, out_dir)
        if results.reaction_order is not None:
            write_reaction_order_csv(results.reaction_order, out_dir)
        if results.coats_redfern is not None:
            write_coats_redfern_csv(results.coats_redfern, out_dir)
        if results.uncertainty is not None:
            write_uncertainty_csv(results.uncertainty, out_dir)
        if results.consistency is not None:
            write_consistency_csv(results.consistency, out_dir)
        if results.compensation is not None:
            write_compensation_csv(results.compensation, out_dir)
        if results.lifetime is not None:
            write_lifetime_csv(results.lifetime, out_dir)
    if config.output.save_plots:
        plot_kwargs = {
            "dpi": config.output.plot_dpi,
            "formats": config.output.plot_formats,
        }
        plot_ea_vs_alpha(results, out_dir, **plot_kwargs)
        plot_kissinger(results, out_dir, **plot_kwargs)
        if results.model_ranking is not None:
            plot_master_plot(results.model_ranking, out_dir, **plot_kwargs)
        if results.preexponential is not None:
            plot_preexp(results.preexponential, out_dir, **plot_kwargs)
        if results.multistep is not None:
            ms_iso = (
                results.isoconversional.get("vyazovkin")
                or results.isoconversional.get("vyazovkin_aic")
                or results.isoconversional.get("kas")
            )
            if ms_iso is not None:
                plot_multistep(ms_iso, results.multistep, out_dir, **plot_kwargs)
        if results.reaction_order is not None:
            plot_reaction_order(results.reaction_order, out_dir, **plot_kwargs)
        if results.coats_redfern is not None:
            plot_coats_redfern(results.coats_redfern, out_dir, **plot_kwargs)
        if results.uncertainty is not None:
            plot_uncertainty(results.uncertainty, out_dir, **plot_kwargs)
        if results.compensation is not None:
            plot_compensation(results.compensation, out_dir, **plot_kwargs)
        if results.lifetime is not None:
            plot_lifetime(results.lifetime, out_dir, **plot_kwargs)
        if config.output.per_method_panels:
            plot_ea_per_method(results, out_dir, **plot_kwargs)

    if config.output.write_markdown_report:
        write_markdown_report(results, case, out_dir)

    print(f"Material: {case.params.material}")
    print(f"Method:   {case.params.method.value}")
    print(f"Type:     {case.params.experiment_type.value}")
    print(f"Rates:    {sorted(case.params.file_to_condition.values())} K/min")
    for name, res in results.isoconversional.items():
        print(f"  {name:14s} mean E_a = {_nanmean(res.Ea_kJ_per_mol):7.2f} kJ/mol")
    if results.kissinger is not None:
        print(
            f"  kissinger      E_a = {results.kissinger.Ea_kJ_per_mol:7.2f} kJ/mol "
            f"(R²={results.kissinger.r_squared:.4f})"
        )
    if results.model_ranking is not None:
        ranked = results.model_ranking.ranked()
        top3 = ", ".join(f"{n} (RMS={d:.3f})" for n, d in ranked[:3])
        print(f"  master_plot    best fit: {results.model_ranking.best_model}  [top-3: {top3}]")
    if results.preexponential is not None:
        pre = results.preexponential
        print(
            f"  preexponential under f(α)={pre.model_name}: "
            f"log10 A = {pre.log10_A_median:.2f} ± {pre.log10_A_mad:.2f}  "
            f"(A ≈ {pre.A_per_sec_median:.3e} 1/s)"
        )
    if results.multistep is not None:
        ms = results.multistep
        print(
            f"  multistep      n={ms.n_steps} step(s), "
            f"E_a flatness {ms.flatness_score*100:.1f}%"
        )
        for s in ms.steps:
            print(
                f"    step {s.index}: α∈[{s.alpha_lo:.2f}, {s.alpha_hi:.2f}], "
                f"E = {s.Ea_kJ_per_mol_median:.1f} ± {s.Ea_kJ_per_mol_mad:.1f} kJ/mol, "
                f"contribution {s.contribution*100:.0f}%"
            )
    if results.reaction_order is not None:
        ro = results.reaction_order
        print(
            f"  reaction_order best n = {ro.n_best:.2f}  "
            f"(E_a = {ro.Ea_best_kJ_per_mol:.1f} kJ/mol, R² = {ro.r_squared_best:.4f})"
        )
    if results.coats_redfern is not None:
        cr = results.coats_redfern
        top = cr.summaries[:3]
        ranked = ", ".join(f"{s.model} (R²={s.r_squared_mean:.3f})" for s in top)
        print(f"  coats_redfern  best model: {cr.best_model}  [top-3: {ranked}]")
    if results.uncertainty is not None:
        unc = results.uncertainty
        finite = np.isfinite(unc.Ea_kJ_per_mol_se)
        if finite.any():
            se_mean = float(np.mean(unc.Ea_kJ_per_mol_se[finite]))
            print(
                f"  uncertainty    {unc.method} jackknife: "
                f"mean SE = {se_mean:.2f} kJ/mol  (n_runs = {unc.n_runs})"
            )
        else:
            print(
                f"  uncertainty    jackknife unavailable "
                f"(need ≥ 3 runs, got {unc.n_runs})"
            )
    if results.endpoint_reliability:
        flagged = [
            (name, rel)
            for name, rel in results.endpoint_reliability.items()
            if rel.warnings
        ]
        if flagged:
            print("  endpoints      α-range/single-step warnings:")
            for _name, rel in flagged:
                for msg in rel.warnings:
                    print(f"    ⚠  {msg}")
    if results.consistency is not None:
        cons = results.consistency
        worst = max(cons.pairs, key=lambda p: p.max_relative_difference, default=None)
        if worst is not None:
            print(
                f"  consistency    {len(cons.pairs)} pair(s), "
                f"max ΔE/E = {worst.max_relative_difference*100:.1f}% "
                f"({worst.method_a} vs {worst.method_b}, α={worst.alpha_of_max:.2f}); "
                f"{len(cons.warnings)} warning(s)"
            )
        for msg in cons.warnings:
            print(f"    ⚠  {msg}")
    if results.compensation is not None:
        comp = results.compensation
        print(
            f"  compensation   {comp.source}: ln A = {comp.slope:.3f}·E + "
            f"{comp.intercept:.2f}, R² = {comp.r_squared:.3f} "
            f"(n = {comp.n_points})"
        )
    for label, emp in results.empirical_fits.items():
        params = ", ".join(f"{k}={v:.3f}" for k, v in emp.parameters.items())
        print(
            f"  {label:14s} {emp.name}: {params}, "
            f"RMS = {emp.rms:.4f}, R² = {emp.r_squared:.4f}"
        )
    if results.lifetime is not None:
        lt = results.lifetime
        print(f"  lifetime       isothermal time-to-α (f(α) = {lt.predictions[0].model_name}):")
        header = "       T (°C)" + "".join(
            f"   t(α={a:.2f})" for a in lt.alpha_targets
        )
        print(header)
        for i, pred in enumerate(lt.predictions):
            row = f"       {pred.T_K - 273.15:6.1f}"
            for j in range(len(lt.alpha_targets)):
                t = lt.times_at_targets[i, j]
                row += f"   {_format_time(t):>10s}"
            print(row)
    print(f"Outputs written to: {out_dir.resolve()}")
    return 0


def _format_time(t: float) -> str:
    if not np.isfinite(t):
        return "—"
    if t < 60.0:
        return f"{t:.1f}s"
    if t < 3600.0:
        return f"{t / 60.0:.1f}m"
    if t < 86400.0:
        return f"{t / 3600.0:.1f}h"
    if t < 365.25 * 86400.0:
        return f"{t / 86400.0:.1f}d"
    return f"{t / (365.25 * 86400.0):.1f}y"


def _cmd_generate(args: argparse.Namespace) -> int:
    case = generate_case(
        rates_K_per_min=args.rates,
        Ea_J_per_mol=args.ea * 1000.0,
        A_per_sec=args.A,
        model=args.model,
        T_start=args.T_start,
        T_stop=args.T_stop,
        noise_std=args.noise,
        seed=args.seed,
        material=f"synthetic-{args.model}",
    )
    write_case(case, args.output)
    print(f"Synthetic case written to {args.output.resolve()}")
    return 0


def _nanmean(arr) -> float:
    arr = np.asarray(arr)
    valid = ~np.isnan(arr)
    return float(np.mean(arr[valid])) if valid.any() else float("nan")


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())

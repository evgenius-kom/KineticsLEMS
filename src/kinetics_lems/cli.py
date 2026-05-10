"""Command-line entry points for KineticsLEMS."""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

from .config import DEFAULT_CONFIG_PATH, load_config
from .io import load_case
from .reporting import plot_ea_vs_alpha, plot_kissinger, write_csv
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
    if config.output.save_plots:
        plot_ea_vs_alpha(results, out_dir, dpi=config.output.plot_dpi)
        plot_kissinger(results, out_dir, dpi=config.output.plot_dpi)

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
    print(f"Outputs written to: {out_dir.resolve()}")
    return 0


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
    import numpy as np

    arr = np.asarray(arr)
    valid = ~np.isnan(arr)
    return float(np.mean(arr[valid])) if valid.any() else float("nan")


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())

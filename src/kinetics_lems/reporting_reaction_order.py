"""CSV + figure exports for the reaction-order n sweep."""
from __future__ import annotations

import csv
from collections.abc import Iterable
from pathlib import Path

import matplotlib.pyplot as plt

from .methods import ReactionOrderResult
from .plotting import DEFAULT_FORMATS, paper_style, save_figure


def write_reaction_order_csv(result: ReactionOrderResult, output_dir: Path) -> Path:
    """R² and recovered E_a as a function of candidate n."""
    output_dir.mkdir(parents=True, exist_ok=True)
    path = output_dir / "reaction_order.csv"
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow([f"# n_best = {result.n_best:.3f}"])
        w.writerow([f"# Ea_best_kJ_per_mol = {result.Ea_best_kJ_per_mol:.3f}"])
        w.writerow([f"# r_squared_best = {result.r_squared_best:.6f}"])
        w.writerow([])
        w.writerow(["n", "r_squared", "Ea_kJ_per_mol"])
        for n, r, e in zip(
            result.n_grid, result.r_squared, result.Ea_kJ_per_mol, strict=True
        ):
            w.writerow([f"{n:.3f}", f"{r:.6f}", f"{e:.4f}"])
    return path


def plot_reaction_order(
    result: ReactionOrderResult,
    output_dir: Path,
    *,
    formats: Iterable[str] = DEFAULT_FORMATS,
    dpi: int = 300,
) -> Path | None:
    """R²(n) curve with the best-n marked."""
    output_dir.mkdir(parents=True, exist_ok=True)
    with paper_style(dpi=dpi):
        fig, (ax_r2, ax_e) = plt.subplots(2, 1, figsize=(5.5, 5.5), sharex=True)
        ax_r2.plot(result.n_grid, result.r_squared, marker="o", markersize=3)
        ax_r2.axvline(
            result.n_best, color="0.4", linestyle="--", linewidth=0.8,
            label=f"n_best = {result.n_best:.2f}",
        )
        ax_r2.set_ylabel(r"$R^2$ of pooled fit")
        ax_r2.legend(frameon=False, loc="best", fontsize=9)
        ax_r2.set_title(
            f"Reaction order sweep — best fit: n = {result.n_best:.2f}, "
            f"$E_a$ = {result.Ea_best_kJ_per_mol:.1f} kJ/mol"
        )

        ax_e.plot(result.n_grid, result.Ea_kJ_per_mol, marker="o", markersize=3)
        ax_e.axvline(result.n_best, color="0.4", linestyle="--", linewidth=0.8)
        ax_e.set_xlabel("Candidate reaction order, n")
        ax_e.set_ylabel(r"$E_{\mathrm{a}}$ (kJ/mol)")
        fig.tight_layout()
        paths = save_figure(fig, output_dir / "reaction_order", formats=formats)
        return paths[0] if paths else None


__all__ = ["plot_reaction_order", "write_reaction_order_csv"]

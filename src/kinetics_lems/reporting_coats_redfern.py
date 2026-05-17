"""CSV + figure exports for Coats–Redfern model fitting."""
from __future__ import annotations

import csv
from collections.abc import Iterable
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from .methods import CoatsRedfernResult
from .plotting import DEFAULT_FORMATS, paper_style, save_figure


def write_coats_redfern_csv(result: CoatsRedfernResult, output_dir: Path) -> list[Path]:
    """Two CSVs: per-(model, run) fits, and the ranked model summary."""
    output_dir.mkdir(parents=True, exist_ok=True)
    written: list[Path] = []

    fits_path = output_dir / "coats_redfern_fits.csv"
    with fits_path.open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow([
            "model", "rate_K_per_min", "Ea_kJ_per_mol",
            "A_per_sec", "r_squared", "n_points",
        ])
        for fit in result.fits:
            w.writerow([
                fit.model,
                f"{fit.rate_K_per_min:.4f}",
                f"{fit.Ea_kJ_per_mol:.4f}",
                f"{fit.A_per_sec:.6e}",
                f"{fit.r_squared:.6f}",
                fit.n_points,
            ])
    written.append(fits_path)

    summary_path = output_dir / "coats_redfern_ranking.csv"
    with summary_path.open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow([
            "model", "Ea_kJ_per_mol_mean", "Ea_kJ_per_mol_std",
            "log10_A_mean", "log10_A_std", "r_squared_mean", "n_runs",
        ])
        for s in result.summaries:
            w.writerow([
                s.model,
                f"{s.Ea_kJ_per_mol_mean:.4f}",
                f"{s.Ea_kJ_per_mol_std:.4f}",
                f"{s.log10_A_mean:.4f}",
                f"{s.log10_A_std:.4f}",
                f"{s.r_squared_mean:.6f}",
                s.n_runs,
            ])
    written.append(summary_path)
    return written


def plot_coats_redfern(
    result: CoatsRedfernResult,
    output_dir: Path,
    *,
    formats: Iterable[str] = DEFAULT_FORMATS,
    dpi: int = 300,
) -> Path | None:
    """Bar chart of mean R² per model (ranked descending) with E_a annotations."""
    output_dir.mkdir(parents=True, exist_ok=True)
    if not result.summaries:
        return None
    names = [s.model for s in result.summaries]
    r2 = np.array([s.r_squared_mean for s in result.summaries])
    Ea = np.array([s.Ea_kJ_per_mol_mean for s in result.summaries])

    with paper_style(dpi=dpi):
        fig, ax = plt.subplots(figsize=(6.0, 4.0))
        bars = ax.barh(names, r2, color="#1f78b4", edgecolor="black", linewidth=0.5)
        for bar, e in zip(bars, Ea, strict=True):
            ax.text(
                bar.get_width(),
                bar.get_y() + bar.get_height() / 2,
                f"  {e:.0f} kJ/mol",
                va="center",
                fontsize=8,
            )
        ax.set_xlabel(r"Mean $R^2$ across runs")
        ax.set_xlim(0.0, 1.05)
        ax.invert_yaxis()  # best at top
        ax.set_title(f"Coats–Redfern ranking — best fit: {result.best_model}")
        fig.tight_layout()
        paths = save_figure(fig, output_dir / "coats_redfern", formats=formats)
        return paths[0] if paths else None


__all__ = ["plot_coats_redfern", "write_coats_redfern_csv"]

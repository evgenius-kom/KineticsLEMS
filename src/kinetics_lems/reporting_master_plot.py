"""CSV + figure exports for the Criado–Málek Z(α) ranking."""
from __future__ import annotations

import csv
from collections.abc import Iterable
from pathlib import Path

import matplotlib.pyplot as plt

from .methods import ModelRanking
from .plotting import DEFAULT_FORMATS, paper_style, save_figure


def write_master_plot_csv(ranking: ModelRanking, output_dir: Path) -> Path:
    """Wide-format CSV with α, experimental Z, and every master Z(α) curve."""
    output_dir.mkdir(parents=True, exist_ok=True)
    path = output_dir / "master_plot_z.csv"
    columns = ["alpha", "experimental"] + list(ranking.master_curves.keys())
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(columns)
        for i, a in enumerate(ranking.alphas):
            row = [f"{a:.4f}", f"{ranking.experimental_z[i]:.6g}"]
            for name in ranking.master_curves:
                row.append(f"{ranking.master_curves[name][i]:.6g}")
            w.writerow(row)

    # Companion ranking summary.
    summary = output_dir / "master_plot_ranking.csv"
    with summary.open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["model", "rms_distance"])
        for name, dist in ranking.ranked():
            w.writerow([name, f"{dist:.6f}"])
    return path


def plot_master_plot(
    ranking: ModelRanking,
    output_dir: Path,
    *,
    top_n: int = 3,
    formats: Iterable[str] = DEFAULT_FORMATS,
    dpi: int = 300,
) -> Path | None:
    """Plot experimental Z(α) overlaid on the top-N closest master curves.

    Other master curves are drawn faintly in gray to give context without
    cluttering the figure.
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    ranked = ranking.ranked()
    top_names = {name for name, _ in ranked[:top_n]}

    with paper_style(dpi=dpi):
        fig, ax = plt.subplots(figsize=(6.0, 4.5))
        # Faint background for the rest.
        for name, curve in ranking.master_curves.items():
            if name in top_names:
                continue
            ax.plot(ranking.alphas, curve, color="#cccccc", linewidth=0.6, zorder=1)
        # Highlighted top-N.
        for name, _ in ranked[:top_n]:
            ax.plot(
                ranking.alphas,
                ranking.master_curves[name],
                linestyle="--",
                linewidth=1.2,
                label=f"{name} (RMS={ranking.rms_distance[name]:.3f})",
                zorder=2,
            )
        # Experimental on top.
        ax.plot(
            ranking.alphas,
            ranking.experimental_z,
            color="black",
            marker="o",
            markersize=4,
            linewidth=1.5,
            label="experimental",
            zorder=3,
        )
        ax.set_xlabel(r"Conversion, $\alpha$")
        ax.set_ylabel(r"$Z(\alpha) / Z(0.5)$")
        ax.set_title(f"Master plot — best fit: {ranking.best_model}")
        ax.set_xlim(0.0, 1.0)
        ax.set_ylim(bottom=0.0)
        ax.legend(frameon=False, loc="best", fontsize=8)
        fig.tight_layout()
        paths = save_figure(fig, output_dir / "master_plot_z", formats=formats)
        return paths[0] if paths else None


__all__ = ["plot_master_plot", "write_master_plot_csv"]

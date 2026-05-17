"""CSV + figure exports for jackknife-by-run E(α) uncertainty."""
from __future__ import annotations

import csv
from collections.abc import Iterable
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from .methods import UncertaintyResult
from .plotting import DEFAULT_FORMATS, paper_style, save_figure


def write_uncertainty_csv(result: UncertaintyResult, output_dir: Path) -> Path:
    output_dir.mkdir(parents=True, exist_ok=True)
    path = output_dir / "uncertainty.csv"
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow([f"# method = {result.method}"])
        w.writerow([f"# n_runs = {result.n_runs}"])
        w.writerow([])
        w.writerow([
            "alpha", "Ea_mean_kJ_per_mol", "Ea_se_kJ_per_mol",
            "Ea_ci95_low_kJ_per_mol", "Ea_ci95_high_kJ_per_mol",
        ])
        for a, m, s, lo, hi in zip(
            result.alpha,
            result.Ea_kJ_per_mol_mean,
            result.Ea_kJ_per_mol_se,
            result.Ea_kJ_per_mol_ci95_low,
            result.Ea_kJ_per_mol_ci95_high,
            strict=True,
        ):
            w.writerow([
                f"{a:.4f}", f"{m:.4f}", f"{s:.4f}", f"{lo:.4f}", f"{hi:.4f}",
            ])
    return path


def plot_uncertainty(
    result: UncertaintyResult,
    output_dir: Path,
    *,
    formats: Iterable[str] = DEFAULT_FORMATS,
    dpi: int = 300,
) -> Path | None:
    """E(α) with shaded 95% jackknife CI."""
    output_dir.mkdir(parents=True, exist_ok=True)
    mean = result.Ea_kJ_per_mol_mean
    valid = np.isfinite(mean) & np.isfinite(result.Ea_kJ_per_mol_se)
    if not valid.any():
        return None
    a = result.alpha[valid]

    with paper_style(dpi=dpi):
        fig, ax = plt.subplots(figsize=(6.0, 4.0))
        ax.plot(a, mean[valid], "o-", markersize=4, linewidth=1.2, label="mean E(α)")
        ax.fill_between(
            a,
            result.Ea_kJ_per_mol_ci95_low[valid],
            result.Ea_kJ_per_mol_ci95_high[valid],
            alpha=0.25,
            label="95% jackknife CI",
        )
        ax.set_xlabel(r"Conversion, $\alpha$")
        ax.set_ylabel(r"$E_{\mathrm{a}}$ (kJ/mol)")
        ax.set_title(
            f"Jackknife uncertainty on {result.method} (n_runs = {result.n_runs})"
        )
        ax.set_xlim(0.0, 1.0)
        ax.legend(frameon=False, loc="best")
        fig.tight_layout()
        paths = save_figure(fig, output_dir / "uncertainty", formats=formats)
        return paths[0] if paths else None


__all__ = ["plot_uncertainty", "write_uncertainty_csv"]

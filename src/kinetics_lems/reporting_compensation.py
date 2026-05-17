"""CSV + plot exports for the compensation-effect diagnostic."""
from __future__ import annotations

import csv
from collections.abc import Iterable
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from .methods import CompensationFit
from .plotting import DEFAULT_FORMATS, paper_style, save_figure


def write_compensation_csv(fit: CompensationFit, output_dir: Path) -> Path:
    output_dir.mkdir(parents=True, exist_ok=True)
    path = output_dir / "compensation.csv"
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow([f"# source = {fit.source}"])
        w.writerow([f"# slope (ln A per kJ/mol) = {fit.slope:.6f}"])
        w.writerow([f"# intercept = {fit.intercept:.6f}"])
        w.writerow([f"# r_squared = {fit.r_squared:.6f}"])
        w.writerow([f"# n_points = {fit.n_points}"])
        w.writerow([])
        w.writerow(["Ea_kJ_per_mol", "ln_A_per_sec"])
        for e, a in zip(fit.E_kJ_per_mol, fit.ln_A_per_sec, strict=True):
            w.writerow([f"{e:.4f}", f"{a:.6f}"])
    return path


def plot_compensation(
    fit: CompensationFit,
    output_dir: Path,
    *,
    formats: Iterable[str] = DEFAULT_FORMATS,
    dpi: int = 300,
) -> Path | None:
    if fit.n_points < 2:
        return None
    output_dir.mkdir(parents=True, exist_ok=True)
    with paper_style(dpi=dpi):
        fig, ax = plt.subplots(figsize=(5.5, 4.0))
        ax.plot(fit.E_kJ_per_mol, fit.ln_A_per_sec, "o", color="#1f78b4", markersize=4)
        xs = np.linspace(fit.E_kJ_per_mol.min(), fit.E_kJ_per_mol.max(), 100)
        ax.plot(
            xs,
            fit.slope * xs + fit.intercept,
            "-",
            color="#e31a1c",
            label=f"fit: ln A = {fit.slope:.3f}·E + {fit.intercept:.2f}  "
                  f"(R² = {fit.r_squared:.3f})",
        )
        ax.set_xlabel(r"$E_{\mathrm{a}}$ (kJ/mol)")
        ax.set_ylabel(r"$\ln A$  (A in 1/s)")
        ax.set_title(f"Compensation effect — source: {fit.source}")
        ax.legend(frameon=False, loc="best", fontsize=8)
        fig.tight_layout()
        paths = save_figure(fig, output_dir / "compensation", formats=formats)
        return paths[0] if paths else None


__all__ = ["plot_compensation", "write_compensation_csv"]

"""CSV + plot exports for analysis results.

Plot styling and multi-format saving live in :mod:`kinetics_lems.plotting`;
this module orchestrates per-method CSVs and the standard set of figures
(E_a vs α overlay, individual per-method panels, Kissinger fit).
"""
from __future__ import annotations

import csv
from collections.abc import Iterable
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from .plotting import DEFAULT_FORMATS, paper_style, save_figure
from .runner import AnalysisResults

# ---------- CSV ----------

def write_csv(results: AnalysisResults, output_dir: Path) -> list[Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    written: list[Path] = []
    for name, res in results.isoconversional.items():
        path = output_dir / f"{name}.csv"
        with path.open("w", newline="", encoding="utf-8") as f:
            w = csv.writer(f)
            w.writerow(["alpha", "Ea_kJ_per_mol", "intercept", "r_squared"])
            for a, e, b, r in zip(
                res.alpha, res.Ea_kJ_per_mol, res.intercept, res.r_squared, strict=True
            ):
                w.writerow([f"{a:.4f}", f"{e:.4f}", f"{b:.6g}", f"{r:.6f}"])
        written.append(path)

    if results.kissinger is not None:
        path = output_dir / "kissinger.csv"
        kiss = results.kissinger
        with path.open("w", newline="", encoding="utf-8") as f:
            w = csv.writer(f)
            w.writerow(["rate_K_per_min", "Tp_K"])
            for b, T in zip(kiss.rates_K_per_min, kiss.Tp_K, strict=True):
                w.writerow([f"{b:.4f}", f"{T:.4f}"])
            w.writerow([])
            w.writerow(["Ea_kJ_per_mol", "A_per_sec_first_order", "r_squared"])
            w.writerow([
                f"{kiss.Ea_kJ_per_mol:.4f}",
                f"{kiss.pre_exponential:.4e}",
                f"{kiss.r_squared:.6f}",
            ])
        written.append(path)
    return written


# ---------- Plots ----------

def plot_ea_vs_alpha(
    results: AnalysisResults,
    output_dir: Path,
    *,
    formats: Iterable[str] = DEFAULT_FORMATS,
    dpi: int = 300,
) -> Path | None:
    """Single overlay figure of E_a(α) for every isoconversional method."""
    if not results.isoconversional:
        return None
    output_dir.mkdir(parents=True, exist_ok=True)
    with paper_style(dpi=dpi):
        fig, ax = plt.subplots(figsize=(6.0, 4.0))
        for name, res in results.isoconversional.items():
            valid = ~np.isnan(res.Ea_kJ_per_mol)
            ax.plot(
                res.alpha[valid],
                res.Ea_kJ_per_mol[valid],
                marker="o",
                markersize=4,
                label=name,
            )
        ax.set_xlabel(r"Conversion, $\alpha$")
        ax.set_ylabel(r"Activation energy, $E_{\mathrm{a}}$ (kJ/mol)")
        ax.set_xlim(0.0, 1.0)
        ax.legend(frameon=False, loc="best")
        fig.tight_layout()
        paths = save_figure(fig, output_dir / "Ea_vs_alpha", formats=formats)
        return paths[0] if paths else None


def plot_ea_per_method(
    results: AnalysisResults,
    output_dir: Path,
    *,
    formats: Iterable[str] = DEFAULT_FORMATS,
    dpi: int = 300,
) -> list[Path]:
    """One figure per isoconversional method (cleaner for journal figures)."""
    if not results.isoconversional:
        return []
    output_dir.mkdir(parents=True, exist_ok=True)
    written: list[Path] = []
    with paper_style(dpi=dpi):
        for name, res in results.isoconversional.items():
            fig, ax = plt.subplots(figsize=(5.0, 3.5))
            valid = ~np.isnan(res.Ea_kJ_per_mol)
            ax.plot(res.alpha[valid], res.Ea_kJ_per_mol[valid], marker="o", markersize=4)
            ax.set_xlabel(r"Conversion, $\alpha$")
            ax.set_ylabel(r"$E_{\mathrm{a}}$ (kJ/mol)")
            ax.set_title(name)
            ax.set_xlim(0.0, 1.0)
            fig.tight_layout()
            written.extend(save_figure(fig, output_dir / f"Ea_{name}", formats=formats))
    return written


def plot_kissinger(
    results: AnalysisResults,
    output_dir: Path,
    *,
    formats: Iterable[str] = DEFAULT_FORMATS,
    dpi: int = 300,
) -> Path | None:
    """Kissinger linearization plot ln(β/T_p²) vs 1/T_p."""
    if results.kissinger is None:
        return None
    output_dir.mkdir(parents=True, exist_ok=True)
    kiss = results.kissinger
    Tp = kiss.Tp_K
    betas_per_sec = kiss.rates_K_per_min / 60.0
    x = 1.0 / Tp
    y = np.log(betas_per_sec / Tp**2)
    slope, intercept = np.polyfit(x, y, 1)

    with paper_style(dpi=dpi):
        fig, ax = plt.subplots(figsize=(5.0, 4.0))
        ax.plot(x, y, "o", markersize=6, label="experimental")
        xs = np.linspace(x.min(), x.max(), 100)
        ax.plot(xs, slope * xs + intercept, "-", label=f"linear fit ($R^2$={kiss.r_squared:.4f})")
        ax.set_xlabel(r"$1 / T_{\mathrm{p}}$ (K$^{-1}$)")
        ax.set_ylabel(r"$\ln(\beta / T_{\mathrm{p}}^{2})$")
        ax.set_title(f"Kissinger: $E_{{\\mathrm{{a}}}}$ = {kiss.Ea_kJ_per_mol:.1f} kJ/mol")
        ax.legend(frameon=False, loc="best")
        fig.tight_layout()
        paths = save_figure(fig, output_dir / "kissinger", formats=formats)
        return paths[0] if paths else None


__all__ = [
    "plot_ea_per_method",
    "plot_ea_vs_alpha",
    "plot_kissinger",
    "write_csv",
]

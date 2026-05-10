"""Plot and CSV export helpers for analysis results."""
from __future__ import annotations

import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # safe default; CLI never opens a window
import matplotlib.pyplot as plt
import numpy as np

from .runner import AnalysisResults


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


def plot_ea_vs_alpha(results: AnalysisResults, output_dir: Path, dpi: int = 120) -> Path | None:
    """One overlay plot of E_a(α) for every isoconversional method."""
    if not results.isoconversional:
        return None
    output_dir.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(8, 5), dpi=dpi)
    for name, res in results.isoconversional.items():
        valid = ~np.isnan(res.Ea_kJ_per_mol)
        ax.plot(res.alpha[valid], res.Ea_kJ_per_mol[valid], marker="o", label=name)
    ax.set_xlabel("Conversion α")
    ax.set_ylabel("Activation energy, kJ/mol")
    ax.set_title("Isoconversional E_a(α)")
    ax.grid(True, alpha=0.3)
    ax.legend()
    path = output_dir / "Ea_vs_alpha.png"
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)
    return path


def plot_kissinger(results: AnalysisResults, output_dir: Path, dpi: int = 120) -> Path | None:
    if results.kissinger is None:
        return None
    output_dir.mkdir(parents=True, exist_ok=True)
    kiss = results.kissinger
    Tp = kiss.Tp_K
    betas_per_sec = kiss.rates_K_per_min / 60.0
    x = 1.0 / Tp
    y = np.log(betas_per_sec / Tp**2)
    slope, intercept = np.polyfit(x, y, 1)
    fig, ax = plt.subplots(figsize=(7, 5), dpi=dpi)
    ax.plot(x, y, "o", label="data")
    xs = np.linspace(x.min(), x.max(), 100)
    ax.plot(xs, slope * xs + intercept, "-", label=f"fit, R²={kiss.r_squared:.4f}")
    ax.set_xlabel("1 / T_p, 1/K")
    ax.set_ylabel("ln(β / T_p²)")
    ax.set_title(f"Kissinger: E_a = {kiss.Ea_kJ_per_mol:.2f} kJ/mol")
    ax.grid(True, alpha=0.3)
    ax.legend()
    path = output_dir / "kissinger.png"
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)
    return path


__all__ = ["write_csv", "plot_ea_vs_alpha", "plot_kissinger"]

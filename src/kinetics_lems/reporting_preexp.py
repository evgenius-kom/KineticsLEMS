"""CSV + figure exports for the pre-exponential A computation."""
from __future__ import annotations

import csv
from collections.abc import Iterable
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from .methods import PreexponentialResult
from .plotting import DEFAULT_FORMATS, paper_style, save_figure


def write_preexp_csv(result: PreexponentialResult, output_dir: Path) -> Path:
    """Per-α A values plus the summary."""
    output_dir.mkdir(parents=True, exist_ok=True)
    path = output_dir / "preexponential.csv"
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow([f"# model = {result.model_name}"])
        w.writerow([f"# log10(A) median = {result.log10_A_median:.4f}"])
        w.writerow([f"# log10(A) MAD    = {result.log10_A_mad:.4f}"])
        w.writerow([f"# A median (1/s)  = {result.A_per_sec_median:.6e}"])
        w.writerow([])
        w.writerow(["alpha", "A_per_sec_median_across_runs"])
        for a, A in zip(result.alpha, result.A_per_sec_per_alpha, strict=True):
            w.writerow([f"{a:.4f}", f"{A:.6e}"])
    return path


def plot_preexp(
    result: PreexponentialResult,
    output_dir: Path,
    *,
    formats: Iterable[str] = DEFAULT_FORMATS,
    dpi: int = 300,
) -> Path | None:
    """log10 A(α) with median band (±MAD)."""
    output_dir.mkdir(parents=True, exist_ok=True)
    A = result.A_per_sec_per_alpha
    valid = np.isfinite(A) & (A > 0)
    if not valid.any():
        return None
    log_A = np.log10(A[valid])
    a = result.alpha[valid]

    with paper_style(dpi=dpi):
        fig, ax = plt.subplots(figsize=(5.5, 4.0))
        ax.plot(a, log_A, marker="o", markersize=4, label="median across runs")
        ax.axhline(result.log10_A_median, color="0.4", linestyle="--", linewidth=0.8,
                   label=f"global median {result.log10_A_median:.2f}")
        ax.fill_between(
            a,
            result.log10_A_median - result.log10_A_mad,
            result.log10_A_median + result.log10_A_mad,
            color="0.7",
            alpha=0.25,
            label=f"± MAD {result.log10_A_mad:.2f}",
        )
        ax.set_xlabel(r"Conversion, $\alpha$")
        ax.set_ylabel(r"$\log_{10} A$ (1/s)")
        ax.set_title(f"Pre-exponential under $f(\\alpha)$ = {result.model_name}")
        ax.set_xlim(0.0, 1.0)
        ax.legend(frameon=False, loc="best", fontsize=8)
        fig.tight_layout()
        paths = save_figure(fig, output_dir / "preexponential", formats=formats)
        return paths[0] if paths else None


__all__ = ["plot_preexp", "write_preexp_csv"]

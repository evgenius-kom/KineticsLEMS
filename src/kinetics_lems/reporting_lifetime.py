"""CSV + figure exports for predictive isothermal α(t)."""
from __future__ import annotations

import csv
from collections.abc import Iterable
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from .methods import LifetimeSummary
from .plotting import DEFAULT_FORMATS, paper_style, save_figure

_SEC_PER_DAY = 86400.0
_SEC_PER_YEAR = 365.25 * _SEC_PER_DAY


def _format_seconds(t: float) -> str:
    """Render seconds in a human-friendly unit."""
    if not np.isfinite(t):
        return "inf"
    if t < 60.0:
        return f"{t:.2f} s"
    if t < 3600.0:
        return f"{t / 60.0:.2f} min"
    if t < _SEC_PER_DAY:
        return f"{t / 3600.0:.2f} h"
    if t < _SEC_PER_YEAR:
        return f"{t / _SEC_PER_DAY:.2f} d"
    return f"{t / _SEC_PER_YEAR:.2f} y"


def write_lifetime_csv(summary: LifetimeSummary, output_dir: Path) -> list[Path]:
    """One CSV with time-to-α table, one per-T α(t) curve."""
    output_dir.mkdir(parents=True, exist_ok=True)
    written: list[Path] = []

    times_path = output_dir / "lifetime_times.csv"
    with times_path.open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        header = ["T_K", "T_C"] + [
            f"t_alpha_{a:.2f}_sec" for a in summary.alpha_targets
        ] + [f"t_alpha_{a:.2f}_human" for a in summary.alpha_targets]
        w.writerow(header)
        for i, pred in enumerate(summary.predictions):
            row = [f"{pred.T_K:.2f}", f"{pred.T_K - 273.15:.2f}"]
            row.extend(f"{summary.times_at_targets[i, j]:.4e}"
                       for j in range(len(summary.alpha_targets)))
            row.extend(_format_seconds(summary.times_at_targets[i, j])
                       for j in range(len(summary.alpha_targets)))
            w.writerow(row)
    written.append(times_path)

    curves_path = output_dir / "lifetime_curves.csv"
    with curves_path.open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["T_K", "alpha", "time_sec"])
        for pred in summary.predictions:
            for a, t in zip(pred.alpha, pred.time_sec, strict=True):
                w.writerow([f"{pred.T_K:.2f}", f"{a:.4f}", f"{t:.6e}"])
    written.append(curves_path)
    return written


def plot_lifetime(
    summary: LifetimeSummary,
    output_dir: Path,
    *,
    formats: Iterable[str] = DEFAULT_FORMATS,
    dpi: int = 300,
) -> Path | None:
    """α vs log10(t) curves, one line per isothermal T."""
    output_dir.mkdir(parents=True, exist_ok=True)
    if not summary.predictions:
        return None

    with paper_style(dpi=dpi):
        fig, ax = plt.subplots(figsize=(6.0, 4.0))
        for pred in summary.predictions:
            t = pred.time_sec
            mask = t > 0
            ax.semilogx(
                t[mask],
                pred.alpha[mask],
                marker="o",
                markersize=3,
                linewidth=1.2,
                label=f"T = {pred.T_K - 273.15:.0f} °C",
            )
        ax.set_xlabel("Time (s)")
        ax.set_ylabel(r"Conversion, $\alpha$")
        ax.set_ylim(0.0, 1.0)
        ax.set_title(
            f"Isothermal prediction (f(α) = {summary.predictions[0].model_name}, "
            f"A = {summary.predictions[0].A_per_sec:.2e} /s)"
        )
        ax.legend(frameon=False, loc="best", fontsize=9)
        fig.tight_layout()
        paths = save_figure(fig, output_dir / "lifetime", formats=formats)
        return paths[0] if paths else None


__all__ = ["plot_lifetime", "write_lifetime_csv"]

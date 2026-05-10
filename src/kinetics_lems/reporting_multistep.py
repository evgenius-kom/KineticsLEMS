"""CSV + figure exports for multi-step segmentation."""
from __future__ import annotations

import csv
from collections.abc import Iterable
from pathlib import Path

import matplotlib.pyplot as plt

from .methods import IsoconversionalResult, MultiStepResult
from .plotting import DEFAULT_FORMATS, paper_style, save_figure


def write_multistep_csv(result: MultiStepResult, output_dir: Path) -> Path:
    output_dir.mkdir(parents=True, exist_ok=True)
    path = output_dir / "multistep.csv"
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow([f"# n_steps = {result.n_steps}"])
        w.writerow([f"# flatness_score = {result.flatness_score:.4f}"])
        w.writerow([])
        w.writerow([
            "step", "alpha_lo", "alpha_hi", "Ea_kJ_per_mol_median",
            "Ea_kJ_per_mol_mad", "contribution",
        ])
        for s in result.steps:
            w.writerow([
                s.index,
                f"{s.alpha_lo:.4f}",
                f"{s.alpha_hi:.4f}",
                f"{s.Ea_kJ_per_mol_median:.4f}",
                f"{s.Ea_kJ_per_mol_mad:.4f}",
                f"{s.contribution:.4f}",
            ])
    return path


def plot_multistep(
    iso: IsoconversionalResult,
    result: MultiStepResult,
    output_dir: Path,
    *,
    formats: Iterable[str] = DEFAULT_FORMATS,
    dpi: int = 300,
) -> Path | None:
    """E(α) line with shaded segment bands and per-segment medians."""
    output_dir.mkdir(parents=True, exist_ok=True)
    with paper_style(dpi=dpi):
        fig, ax = plt.subplots(figsize=(6.0, 4.0))
        ax.plot(iso.alpha, iso.Ea_kJ_per_mol, "o-", color="black",
                markersize=3, linewidth=1.0, label="E(α)")
        for k, s in enumerate(result.steps):
            ax.axvspan(s.alpha_lo, s.alpha_hi, alpha=0.10 + 0.05 * (k % 4))
            ax.hlines(
                s.Ea_kJ_per_mol_median,
                s.alpha_lo,
                s.alpha_hi,
                linestyles="--",
                linewidth=1.4,
                label=f"step {s.index}: {s.Ea_kJ_per_mol_median:.1f} kJ/mol",
            )
        ax.set_xlabel(r"Conversion, $\alpha$")
        ax.set_ylabel(r"$E_{\mathrm{a}}$ (kJ/mol)")
        ax.set_title(
            f"Multi-step segmentation — {result.n_steps} step(s), "
            f"flatness {result.flatness_score*100:.1f}%"
        )
        ax.set_xlim(0.0, 1.0)
        ax.legend(frameon=False, loc="best", fontsize=8)
        fig.tight_layout()
        paths = save_figure(fig, output_dir / "multistep", formats=formats)
        return paths[0] if paths else None


__all__ = ["plot_multistep", "write_multistep_csv"]

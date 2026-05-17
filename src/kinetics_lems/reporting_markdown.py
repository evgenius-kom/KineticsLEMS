"""Single-file Markdown report assembling every analysis output.

The report is self-contained: it cross-references the per-method CSVs and
PNG plots already written by :mod:`reporting_*`, so dropping the entire
output directory into a workshare folder gives a reviewer everything
they need.

Sections:
    1. Run summary (material, technique, rates, config snapshot).
    2. Diagnostics (endpoint warnings + ICTAC consistency).
    3. Isoconversional methods — table of mean E(α).
    4. Model identification — master_plot, Coats-Redfern, empirical.
    5. Kinetic triplet — A, multistep, compensation.
    6. Lifetime predictions.
    7. References to written artefacts.
"""
from __future__ import annotations

from datetime import UTC, datetime
from pathlib import Path

import numpy as np

from .models import CaseData
from .runner import AnalysisResults


def write_markdown_report(
    results: AnalysisResults,
    case: CaseData,
    output_dir: Path,
    *,
    filename: str = "report.md",
) -> Path:
    """Write a Markdown report into ``output_dir / filename``."""
    output_dir.mkdir(parents=True, exist_ok=True)
    path = output_dir / filename

    lines: list[str] = []
    lines.append(f"# KineticsLEMS analysis report — {case.params.material}")
    now = datetime.now(UTC).strftime("%Y-%m-%d %H:%M UTC")
    lines.append(f"_Generated {now}_")
    lines.append("")

    _section_run_summary(lines, case)
    _section_diagnostics(lines, results)
    _section_isoconversional(lines, results)
    _section_model_identification(lines, results)
    _section_kinetic_triplet(lines, results)
    _section_lifetime(lines, results)
    _section_artefacts(lines, output_dir)

    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


# --- Sections ---------------------------------------------------------------


def _section_run_summary(lines: list[str], case: CaseData) -> None:
    lines.append("## 1. Run summary")
    lines.append("")
    lines.append(f"- **Material:** {case.params.material}")
    lines.append(f"- **Technique:** {case.params.method.value}")
    lines.append(f"- **Experiment:** {case.params.experiment_type.value}")
    rates = sorted(case.params.file_to_condition.values())
    lines.append(f"- **Heating rates (K/min):** {', '.join(f'{r:g}' for r in rates)}")
    lines.append("")


def _section_diagnostics(lines: list[str], results: AnalysisResults) -> None:
    lines.append("## 2. Diagnostics")
    lines.append("")
    if results.endpoint_reliability:
        any_warnings = False
        for _name, rel in results.endpoint_reliability.items():
            if rel.warnings:
                if not any_warnings:
                    lines.append("**Endpoint / single-step warnings:**")
                    lines.append("")
                    any_warnings = True
                for msg in rel.warnings:
                    lines.append(f"- {msg}")
        if not any_warnings:
            lines.append("All isoconversional methods passed endpoint and "
                         "single-step checks.")
        lines.append("")
    if results.consistency is not None:
        c = results.consistency
        lines.append(
            f"**ICTAC consistency check** (threshold "
            f"{c.threshold * 100:.0f}%):"
        )
        lines.append("")
        lines.append("| Method A | Method B | Max ΔE/E | α of max | Exceeds? |")
        lines.append("|---|---|---:|---:|:---:|")
        for p in c.pairs:
            lines.append(
                f"| {p.method_a} | {p.method_b} | "
                f"{p.max_relative_difference * 100:.1f}% | "
                f"{p.alpha_of_max:.2f} | "
                f"{'⚠' if p.max_relative_difference > c.threshold else '—'} |"
            )
        for msg in c.warnings:
            lines.append(f"> ⚠ {msg}")
        lines.append("")


def _section_isoconversional(lines: list[str], results: AnalysisResults) -> None:
    lines.append("## 3. Isoconversional E(α)")
    lines.append("")
    if not results.isoconversional:
        lines.append("_No isoconversional methods were run._")
        lines.append("")
        return
    lines.append("| Method | Mean E_a (kJ/mol) | E_a flatness (core) |")
    lines.append("|---|---:|---:|")
    for name, res in results.isoconversional.items():
        finite = np.isfinite(res.Ea_kJ_per_mol)
        mean = float(np.mean(res.Ea_kJ_per_mol[finite])) if finite.any() else float("nan")
        rel = results.endpoint_reliability.get(name)
        flatness = (
            f"{rel.flatness_in_core * 100:.1f}%"
            if rel is not None and np.isfinite(rel.flatness_in_core)
            else "—"
        )
        lines.append(f"| {name} | {mean:.2f} | {flatness} |")
    lines.append("")
    if results.kissinger is not None:
        k = results.kissinger
        lines.append(
            f"**Kissinger:** E_a = {k.Ea_kJ_per_mol:.2f} kJ/mol "
            f"(R² = {k.r_squared:.4f})."
        )
        lines.append("")
    lines.append("![E_a vs α overlay](Ea_vs_alpha.png)")
    lines.append("")


def _section_model_identification(lines: list[str], results: AnalysisResults) -> None:
    lines.append("## 4. Model identification")
    lines.append("")
    if results.model_ranking is not None:
        ranked = results.model_ranking.ranked()
        lines.append("**Master plot Z(α) ranking** (lower RMS = better):")
        lines.append("")
        lines.append("| Model | RMS distance |")
        lines.append("|---|---:|")
        for name, rms in ranked[:6]:
            lines.append(f"| {name} | {rms:.4f} |")
        lines.append("")
        lines.append(f"Best fit: **{results.model_ranking.best_model}**.")
        lines.append("")
        lines.append("![Master plot](master_plot_z.png)")
        lines.append("")
    if results.coats_redfern is not None:
        lines.append("**Coats-Redfern ranking** (mean across heating rates):")
        lines.append("")
        lines.append("| Model | Mean R² | Mean AIC | Mean E_a (kJ/mol) | n_runs |")
        lines.append("|---|---:|---:|---:|---:|")
        for s in results.coats_redfern.summaries[:6]:
            lines.append(
                f"| {s.model} | {s.r_squared_mean:.4f} | {s.aic_mean:.2f} | "
                f"{s.Ea_kJ_per_mol_mean:.2f} | {s.n_runs} |"
            )
        lines.append("")
    if results.empirical_fits:
        lines.append("**Empirical (Prout-Tompkins / Sestak-Berggren) fits:**")
        lines.append("")
        lines.append("| Model | Parameters | R² | RMS |")
        lines.append("|---|---|---:|---:|")
        for emp in results.empirical_fits.values():
            params = ", ".join(f"{k}={v:.3f}" for k, v in emp.parameters.items())
            lines.append(
                f"| {emp.name} | {params} | {emp.r_squared:.4f} | {emp.rms:.4f} |"
            )
        lines.append("")


def _section_kinetic_triplet(lines: list[str], results: AnalysisResults) -> None:
    lines.append("## 5. Kinetic triplet")
    lines.append("")
    if results.preexponential is not None:
        pre = results.preexponential
        lines.append(
            f"- Pre-exponential under f(α) = **{pre.model_name}**: "
            f"log₁₀ A = {pre.log10_A_median:.2f} ± {pre.log10_A_mad:.2f}  "
            f"(A ≈ {pre.A_per_sec_median:.3e} 1/s)."
        )
    if results.multistep is not None:
        ms = results.multistep
        lines.append(
            f"- Multistep segmentation: **{ms.n_steps}** step(s), "
            f"flatness {ms.flatness_score * 100:.1f}%. "
            f"AIC = {ms.aic_piecewise_constant:.2f}, "
            f"BIC = {ms.bic_piecewise_constant:.2f}."
        )
        for s in ms.steps:
            lines.append(
                f"    - step {s.index}: α ∈ [{s.alpha_lo:.2f}, {s.alpha_hi:.2f}], "
                f"E_a = {s.Ea_kJ_per_mol_median:.1f} ± {s.Ea_kJ_per_mol_mad:.1f} kJ/mol "
                f"({s.contribution * 100:.0f}% of total reaction)."
            )
    if results.compensation is not None:
        c = results.compensation
        lines.append(
            f"- Compensation effect ({c.source}): "
            f"ln A = {c.slope:.3f} · E + {c.intercept:.2f}, "
            f"R² = {c.r_squared:.3f}."
        )
    if results.uncertainty is not None:
        u = results.uncertainty
        finite = np.isfinite(u.Ea_kJ_per_mol_se)
        if finite.any():
            se = float(np.mean(u.Ea_kJ_per_mol_se[finite]))
            lines.append(
                f"- Uncertainty ({u.method} jackknife): mean SE = "
                f"{se:.2f} kJ/mol (n = {u.n_runs} runs)."
            )
    lines.append("")


def _section_lifetime(lines: list[str], results: AnalysisResults) -> None:
    lines.append("## 6. Lifetime predictions")
    lines.append("")
    if results.lifetime is None:
        lines.append("_No lifetime predictions were requested._")
        lines.append("")
        return
    lt = results.lifetime
    targets = lt.alpha_targets
    header = "| T (°C) | " + " | ".join(f"t(α={a:.2f})" for a in targets) + " |"
    sep = "|---:|" + "---:|" * len(targets)
    lines.append(header)
    lines.append(sep)
    for i, pred in enumerate(lt.predictions):
        row = [f"{pred.T_K - 273.15:.1f}"]
        for j in range(len(targets)):
            t = lt.times_at_targets[i, j]
            row.append(_format_duration(t))
        lines.append("| " + " | ".join(row) + " |")
    lines.append("")
    lines.append(f"f(α) used: **{lt.predictions[0].model_name}**.")
    lines.append("")


def _section_artefacts(lines: list[str], output_dir: Path) -> None:
    lines.append("## 7. Files in this output directory")
    lines.append("")
    csv_files = sorted(p.name for p in output_dir.glob("*.csv"))
    if csv_files:
        lines.append("**CSV:**")
        for name in csv_files:
            lines.append(f"- `{name}`")
        lines.append("")
    plot_files = sorted(p.name for p in output_dir.glob("*.png"))
    if plot_files:
        lines.append("**Plots (PNG):**")
        for name in plot_files:
            lines.append(f"- `{name}`")
        lines.append("")


def _format_duration(t: float) -> str:
    if not np.isfinite(t):
        return "—"
    if t < 60.0:
        return f"{t:.1f} s"
    if t < 3600.0:
        return f"{t / 60:.1f} min"
    if t < 86400.0:
        return f"{t / 3600:.1f} h"
    if t < 365.25 * 86400.0:
        return f"{t / 86400:.1f} d"
    return f"{t / (365.25 * 86400):.1f} y"


__all__ = ["write_markdown_report"]

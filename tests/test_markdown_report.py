"""Markdown report generation: smoke + section content."""
from __future__ import annotations

from pathlib import Path

import pytest

from kinetics_lems.config import load_config
from kinetics_lems.io import load_case
from kinetics_lems.reporting_markdown import write_markdown_report
from kinetics_lems.runner import run_analysis

EXAMPLE_CASE = (
    Path(__file__).resolve().parent.parent / "examples" / "synthetic" / "F1_120kJ"
)


@pytest.fixture(scope="module")
def f1_results(tmp_path_factory):
    if not EXAMPLE_CASE.is_dir():
        pytest.skip(f"synthetic case missing at {EXAMPLE_CASE}")
    case = load_case(EXAMPLE_CASE)
    cfg = load_config()
    results = run_analysis(case, cfg)
    return case, results


def test_markdown_report_writes_file(tmp_path: Path, f1_results) -> None:
    case, results = f1_results
    path = write_markdown_report(results, case, tmp_path)
    assert path.is_file()
    text = path.read_text(encoding="utf-8")
    assert "# KineticsLEMS analysis report — synthetic-F1" in text
    # All major section headers must appear.
    for header in (
        "## 1. Run summary",
        "## 2. Diagnostics",
        "## 3. Isoconversional E",
        "## 4. Model identification",
        "## 5. Kinetic triplet",
        "## 6. Lifetime predictions",
        "## 7. Files in this output directory",
    ):
        assert header in text, f"missing section: {header}"


def test_markdown_report_contains_heating_rates(tmp_path: Path, f1_results) -> None:
    case, results = f1_results
    path = write_markdown_report(results, case, tmp_path)
    text = path.read_text(encoding="utf-8")
    for rate in case.params.file_to_condition.values():
        assert f"{rate:g}" in text


def test_markdown_report_reports_kissinger(tmp_path: Path, f1_results) -> None:
    case, results = f1_results
    if results.kissinger is None:
        pytest.skip("Kissinger not in results")
    path = write_markdown_report(results, case, tmp_path)
    assert "Kissinger" in path.read_text(encoding="utf-8")

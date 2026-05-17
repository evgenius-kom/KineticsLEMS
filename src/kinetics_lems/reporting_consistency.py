"""CSV export for the ICTAC consistency check."""
from __future__ import annotations

import csv
from pathlib import Path

from .methods import ConsistencyResult


def write_consistency_csv(result: ConsistencyResult, output_dir: Path) -> Path:
    output_dir.mkdir(parents=True, exist_ok=True)
    path = output_dir / "consistency.csv"
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(
            [
                "method_a",
                "method_b",
                "max_relative_difference",
                "alpha_of_max",
                "exceeds_threshold",
            ]
        )
        for p in result.pairs:
            w.writerow(
                [
                    p.method_a,
                    p.method_b,
                    f"{p.max_relative_difference:.6f}",
                    f"{p.alpha_of_max:.4f}",
                    "yes" if p.max_relative_difference > result.threshold else "no",
                ]
            )
        w.writerow([])
        w.writerow(["threshold", f"{result.threshold:.4f}"])
        w.writerow(["warning_count", str(len(result.warnings))])
        for msg in result.warnings:
            w.writerow(["warning", msg])
    return path


__all__ = ["write_consistency_csv"]

"""One-time extraction of ABPOT / DAPOT reference α(T) tables from
``theory/Вязовкин.xlsx`` into ``tests/fixtures/`` so the regression test
no longer needs the xlsx (which can then be deleted).

Layout produced::

    tests/fixtures/
        abpot/
            rate_0.2.tsv     # columns: alpha, T_K
            rate_1.0.tsv
            rate_2.5.tsv
            rate_5.0.tsv
            reference_Ea.tsv # columns: alpha, Ea_kJ_per_mol
            SOURCE.md
        dapot/
            rate_2.5.tsv
            rate_5.0.tsv
            rate_10.0.tsv
            reference_Ea.tsv
            SOURCE.md

Run from repo root::

    .venv/bin/python scripts/extract_reference_fixtures.py
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
from openpyxl import load_workbook

REPO = Path(__file__).resolve().parents[1]
XLSX = REPO / "theory" / "Вязовкин.xlsx"
OUT = REPO / "tests" / "fixtures"


def _read_pairs(sheet, c_alpha: int, c_T: int) -> tuple[np.ndarray, np.ndarray]:
    alpha: list[float] = []
    T: list[float] = []
    for r in range(3, sheet.max_row + 1):
        a = sheet.cell(r, c_alpha).value
        t = sheet.cell(r, c_T).value
        if a is None or t is None:
            continue
        try:
            a_f, t_f = float(a), float(t)
        except (TypeError, ValueError):
            continue
        alpha.append(a_f)
        T.append(t_f)
    a_arr = np.asarray(alpha)
    t_arr = np.asarray(T)
    order = np.argsort(a_arr)
    a_arr = np.clip(a_arr[order], 0.0, 1.0)
    t_arr = t_arr[order]
    return a_arr, t_arr


def _read_reference(sheet, col: int) -> tuple[np.ndarray, np.ndarray]:
    alpha: list[float] = []
    Ea: list[float] = []
    for r in range(3, sheet.max_row + 1):
        a = sheet.cell(r, 1).value  # column 1 is the canonical α grid
        e = sheet.cell(r, col).value
        if a is None or e is None:
            continue
        try:
            a_f, e_f = float(a), float(e)
        except (TypeError, ValueError):
            continue
        if not np.isfinite(e_f):
            continue
        alpha.append(a_f)
        Ea.append(e_f)
    return np.asarray(alpha), np.asarray(Ea)


def _save_tsv(path: Path, header: tuple[str, str], data: tuple[np.ndarray, np.ndarray]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    a, b = data
    with path.open("w", encoding="utf-8") as f:
        f.write(f"{header[0]}\t{header[1]}\n")
        for x, y in zip(a, b, strict=True):
            f.write(f"{x:.6f}\t{y:.6f}\n")


def _save_source(folder: Path, material: str, rates: list[float]) -> None:
    folder.mkdir(parents=True, exist_ok=True)
    (folder / "SOURCE.md").write_text(
        f"""# {material} reference fixtures

Source: ``theory/Вязовкин.xlsx`` sheet ``{material}``.

Each ``rate_<β>.tsv`` file contains the (α, T_K) pairs from a single
heating-rate column of the source sheet, sorted by α and clipped into [0, 1].

``reference_Ea.tsv`` contains the per-α activation energy from the
"Ea, kJ" column of the source sheet — produced by an independent
Excel-based Vyazovkin computation. Used by
``tests/test_reference_abpot_dapot.py`` to detect drift in our integral
methods (KAS / OFW / Vyazovkin) against this independent reference.

Heating rates: {rates} K/min.
""",
        encoding="utf-8",
    )


def main() -> None:
    wb = load_workbook(XLSX, data_only=True)

    # ---- ABPOT ----
    abpot = wb["АБПОТ"]
    abpot_dir = OUT / "abpot"
    abpot_rates = {0.2: (1, 2), 1.0: (4, 5), 2.5: (7, 8), 5.0: (10, 11)}
    for beta, (c_a, c_t) in abpot_rates.items():
        a, t = _read_pairs(abpot, c_a, c_t)
        _save_tsv(abpot_dir / f"rate_{beta}.tsv", ("alpha", "T_K"), (a, t))
    a_ref, ea_ref = _read_reference(abpot, col=13)
    _save_tsv(abpot_dir / "reference_Ea.tsv", ("alpha", "Ea_kJ_per_mol"), (a_ref, ea_ref))
    _save_source(abpot_dir, "АБПОТ", list(abpot_rates))

    # ---- DAPOT ----
    dapot = wb["ДАПОТ"]
    dapot_dir = OUT / "dapot"
    dapot_rates = {2.5: (4, 5), 5.0: (7, 8), 10.0: (10, 11)}
    for beta, (c_a, c_t) in dapot_rates.items():
        a, t = _read_pairs(dapot, c_a, c_t)
        _save_tsv(dapot_dir / f"rate_{beta}.tsv", ("alpha", "T_K"), (a, t))
    a_ref, ea_ref = _read_reference(dapot, col=13)
    _save_tsv(dapot_dir / "reference_Ea.tsv", ("alpha", "Ea_kJ_per_mol"), (a_ref, ea_ref))
    _save_source(dapot_dir, "ДАПОТ", list(dapot_rates))

    print(f"Fixtures written to {OUT}")


if __name__ == "__main__":
    main()

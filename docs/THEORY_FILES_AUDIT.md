# Audit: can `theory/` files be deleted?

Verdict per file (relative to current code + [docs/ALGORITHMS.md](ALGORITHMS.md)):

## DOCX — safe to delete

| File                       | What's in it                                  | Captured? |
|----------------------------|-----------------------------------------------|-----------|
| `Алгоритмы.docx`           | Step-by-step recipe for DSC pre-processing, Kissinger, Friedman, OFW, KAS, Vyazovkin (with the explicit Φ formula and Senum–Yang p(x))           | ✅ All steps and formulas reproduced verbatim in `docs/ALGORITHMS.md` §A–§F |
| `Изоконверс методы.docx`   | Master's thesis covering classification of isoconversional methods (Friedman / Flynn-Ozawa-Wall / Vyazovkin), ASTM E698, comparison & history    | ✅ Algorithmic content captured. Historical/comparative narrative — kept implicitly via references in `docs/ALGORITHMS.md` (Friedman 1964, Ozawa 1965, Vyazovkin 1996/1997/2001, ICTAC 2011) |

These are safe to delete from `theory/`.

## XLSX / XLSM — safe to delete (data already migrated)

| File                                        | What's in it                                                                    | Captured? |
|---------------------------------------------|---------------------------------------------------------------------------------|-----------|
| `Вязовкин.xlsx`                             | ABPOT and DAPOT α(T) tables + per-α reference E_a column (independent Vyazovkin computation) | ✅ Migrated to [`tests/fixtures/abpot/`](../tests/fixtures/abpot/) and [`tests/fixtures/dapot/`](../tests/fixtures/dapot/); regression-tested by [`tests/test_reference_abpot_dapot.py`](../tests/test_reference_abpot_dapot.py) |
| `для расчета по Вязовкину.xlsm`             | Macro-driven workbook for the same kind of Vyazovkin computation                | ❌ Not migrated — same kind of content as `Вязовкин.xlsx`, considered superseded by it. |
| `Сводная таблица.xlsx`                      | Per-experiment DSC summary (mass, β, peak T, integrals, polymerization peak) for ABPOT / DAPOT                  | ❌ Not migrated — experimental log, not numerical reference. Re-extractable from raw lab notebooks if ever needed. |
| `Zalfa plots Pr-Tm program.xlsm`            | Macro-driven α-vs-T plotter for Pr–Tm (separate material)                       | ❌ Not migrated — only useful if you actively work on Pr–Tm; out of current scope. |

**Recommendation.** Safe to delete. The numerically-relevant content from
`Вязовкин.xlsx` (ABPOT + DAPOT reference Eₐ) is now under
[`tests/fixtures/`](../tests/fixtures/) and is exercised by the regression
test on every `pytest` run. The other three workbooks contain
experimental-log data or material-specific analyses that are not used by
the codebase; keep them only if you personally need them as an archive.

## Other files in `theory/` (PDFs and pptx)

Not analyzed (per your instruction "do not waste time reading PDFs"). They
are bibliographic references; the ones cited in `docs/ALGORITHMS.md` are
publicly available, so deleting the local copies has no information cost.
The pptx is a slide deck on isoconversional methods — likely educational.
Both can be deleted at your discretion.

## Bottom line

- **Safe to delete:** all `.docx` and all `.xlsx` / `.xlsm` files.
  The reference Vyazovkin numbers from ABPOT / DAPOT are preserved under
  `tests/fixtures/` and exercised by `pytest`.
- **Optional / your call:** PDFs, pptx.

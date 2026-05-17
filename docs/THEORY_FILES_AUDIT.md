# Audit: can `theory/` files be deleted?

Verdict per file (relative to current code + [docs/ALGORITHMS.md](ALGORITHMS.md)):

## DOCX — safe to delete

| File                       | What's in it                                  | Captured? |
|----------------------------|-----------------------------------------------|-----------|
| `Алгоритмы.docx`           | Step-by-step recipe for DSC pre-processing, Kissinger, Friedman, OFW, KAS, Vyazovkin (with the explicit Φ formula and Senum–Yang p(x))           | ✅ All steps and formulas reproduced verbatim in `docs/ALGORITHMS.md` §A–§F |
| `Изоконверс методы.docx`   | Master's thesis covering classification of isoconversional methods (Friedman / Flynn-Ozawa-Wall / Vyazovkin), ASTM E698, comparison & history    | ✅ Algorithmic content captured. Historical/comparative narrative — kept implicitly via references in `docs/ALGORITHMS.md` (Friedman 1964, Ozawa 1965, Vyazovkin 1996/1997/2001, ICTAC 2011) |

These are safe to delete from `theory/`.

## XLSX / XLSM — safe to delete

All four workbooks were re-extracted to TSV fixtures under
[`data/reference_workbooks/`](../data/reference_workbooks/). Each material
gets its own folder with a per-folder README describing the data and
applicable tests; xlsx-internal computations land in
[`data/reference_workbooks/xlsx_crosscheck/`](../data/reference_workbooks/xlsx_crosscheck/).

| File                                        | What's in it                                                                    | Captured? |
|---------------------------------------------|---------------------------------------------------------------------------------|-----------|
| `Вязовкин.xlsx`                             | ABPOT, DAPOT and Epoxy PV15-0,5 α(T) tables + per-α reference E_a columns (independent Vyazovkin computations) | ✅ Migrated to `data/reference_workbooks/{abpot,dapot,epoxy_pv15}/`; regression-tested by `tests/test_reference_{abpot_dapot,epoxy}.py` |
| `для расчета по Вязовкину.xlsm`             | Macro-driven workbook with multi-rate Vyazovkin computations: 3-rate / 4-rate ABPOT and DAPOT, plus prepreg-СКМ (6 rates) | ✅ The 4th DAPOT rate (β = 1 K/min) and SKM 6-rate dataset extracted; SKM stress test in `tests/test_reference_skm.py` |
| `Сводная таблица.xlsx`                      | DSC experiment log + Z(α) master-plot computations (sheet "Z-plots") + per-α A calculation (sheet "A calculated") + reaction-order check (sheet "Проверка порядка") | ✅ All three xlsx-formula outputs extracted to `data/reference_workbooks/xlsx_crosscheck/`; bit-for-bit verified by `tests/test_xlsx_crosscheck_*.py` |
| `Zalfa plots Pr-Tm program.xlsm`            | Pr–Tm material's single-rate Z(α) computation                                   | ⚠️ Not extracted — Pr–Tm single-rate Z(α) is superseded by the cross-check Z-plot data extracted from `Сводная таблица.xlsx`. Re-extract only if Pr–Tm is needed as a separate test material (see [TODO_FEATURES](TODO_FEATURES.md) item 11). |

The extraction is reproducible via
[`scripts/extract_reference_fixtures.py`](../scripts/extract_reference_fixtures.py)
if the xlsx files are restored under `theory/extra/`.

## Other files in `theory/` (PDFs and pptx)

Not analyzed (per your instruction "do not waste time reading PDFs"). They
are bibliographic references; the ones cited in `docs/ALGORITHMS.md` are
publicly available, so deleting the local copies has no information cost.
The pptx is a slide deck on isoconversional methods — likely educational.
Both can be deleted at your discretion.

## Bottom line

- **Safe to delete:** all `.docx` and all `.xlsx` / `.xlsm` files
  (including `theory/extra/` if present).  The reference Vyazovkin
  numbers, xlsx Z(α)/A/n-order tables, and four reference materials
  (ABPOT, DAPOT, Epoxy, SKM) are preserved under
  `data/reference_workbooks/` and exercised by `pytest`.
- **Optional / your call:** PDFs, pptx.

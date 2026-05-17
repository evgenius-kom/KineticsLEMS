# DAPOT polymer

Source: `theory/extra/Вязовкин.xlsx` sheet `ДАПОТ`
(plus β = 1 K/min from `для расчета по Вязовкину.xlsm` sheet `3 скорости ДАПОТ`).

## Data files

| File | Columns | Notes |
|------|---------|-------|
| `rate_1.0.tsv`, `rate_2.5.tsv`, `rate_5.0.tsv`, `rate_10.0.tsv` | `alpha`, `T_K` | (α, T) pairs per heating rate |
| `reference_Ea.tsv` | `alpha`, `Ea_kJ_per_mol` | per-α reference E_a from independent Excel Vyazovkin |

Heating rates: 1.0, 2.5, 5.0, 10.0 K/min. α grid step ≈ 0.01.

## Applicable tests

* **Integral-method validation** — KAS / OFW / Vyazovkin against reference E_a
  ([`tests/test_reference_abpot_dapot.py`](../../../tests/test_reference_abpot_dapot.py)).
* **Z(α) cross-check at 5 K/min** ([`../xlsx_crosscheck/zplot_dapot_5kpm.tsv`](../xlsx_crosscheck/)).
* **Reaction-order cross-check at 5 K/min** ([`../xlsx_crosscheck/reaction_order_dapot.tsv`](../xlsx_crosscheck/)).

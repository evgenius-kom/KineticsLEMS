# ABPOT polymer

Source: `theory/extra/Вязовкин.xlsx` sheet `АБПОТ`.

## Data files

| File | Columns | Notes |
|------|---------|-------|
| `rate_0.2.tsv`, `rate_1.0.tsv`, `rate_2.5.tsv`, `rate_5.0.tsv` | `alpha`, `T_K` | (α, T) pairs per heating rate, sorted by α, clipped to [0, 1] |
| `reference_Ea.tsv` | `alpha`, `Ea_kJ_per_mol` | per-α reference E_a from an independent Excel-based Vyazovkin computation |

Heating rates: 0.2, 1.0, 2.5, 5.0 K/min. α grid step ≈ 0.01.

## Applicable tests

* **Integral-method validation** — KAS / OFW / Vyazovkin reproduce
  the per-α reference E_a (see [`tests/test_reference_abpot_dapot.py`](../../../tests/test_reference_abpot_dapot.py)).
* **Multi-rate robustness** — 4-rate dataset; useful for sanity-checking
  Vyazovkin and Vyazovkin-AIC.
* **Z(α) cross-check at 5 K/min** — paired with
  [`../xlsx_crosscheck/zplot_abpot_5kpm.tsv`](../xlsx_crosscheck/).
* **A(α) cross-check at rates 1, 2.5, 5** — paired with
  [`../xlsx_crosscheck/A_calculated_abpot.tsv`](../xlsx_crosscheck/).
* **Reaction-order cross-check at 5 K/min** — paired with
  [`../xlsx_crosscheck/reaction_order_abpot.tsv`](../xlsx_crosscheck/).

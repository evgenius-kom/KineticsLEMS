# xlsx cross-check fixtures

Source: `theory/extra/Сводная таблица.xlsx` sheets `Z-plots`,
`A calculated`, `Проверка порядка`.

These tables contain values computed by the original Excel formulas
that motivated each algorithm. Used as **bit-exactness cross-checks**
against our re-implementations: our code reads (α, T, dα/dt) from the
same rows and must reproduce the precomputed Z / A / y(n) numbers.

## Data files

| File | Used by |
|------|---------|
| `zplot_abpot_5kpm.tsv`, `zplot_dapot_5kpm.tsv` | [`tests/test_xlsx_crosscheck_zplot.py`](../../../tests/test_xlsx_crosscheck_zplot.py) — Z(α) and normalized Z(α)/Z(0.5) at 5 K/min vs our `master_plot.py` |
| `A_calculated_abpot.tsv` | [`tests/test_xlsx_crosscheck_preexp.py`](../../../tests/test_xlsx_crosscheck_preexp.py) — A(α, run) under f(α) = 1 − α (F1) vs our `preexponential.py` |
| `reaction_order_abpot.tsv`, `reaction_order_dapot.tsv` | [`tests/test_xlsx_crosscheck_reaction_order.py`](../../../tests/test_xlsx_crosscheck_reaction_order.py) — y(n; α, T) = ln(dα/dT) − n · ln(1 − α) at n ∈ {1, 2, 3, 0.7} vs our `reaction_order.py` |

The cross-check tests use *the xlsx's own (α, T, dα/dt)* as input so
that any disagreement is attributable to a difference in computation,
not in extraction.

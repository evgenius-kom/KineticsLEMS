# SKM preimpregnated polymer

Source: `theory/extra/для расчета по Вязовкину.xlsm` sheet `прпепрег-СКМ`.

Six heating rates: 1, 2.5, 5, 10, 20, 40 K/min — span of 40× between
slowest and fastest. The widest β-range in this corpus, useful for
checking multi-rate robustness of isoconversional methods.

⚠️ **No reference E_a in the source sheet.** Treat as a
*self-consistency* fixture — different methods must agree with each
other across the same α window, but there is no independent truth.

## Data files

| File | Columns | Notes |
|------|---------|-------|
| `rate_<β>.tsv` × 6 | `alpha`, `T_K` | β ∈ {1, 2.5, 5, 10, 20, 40} K/min |

α grid: 0.005 → 0.995, step 0.01.

## Applicable tests

* **Multi-rate stress test** — agreement between KAS, OFW, Vyazovkin
  on a 6-rate dataset within an interior α window.
* **Method consistency** — Vyazovkin vs Vyazovkin-AIC at α ∈ [0.1, 0.9]
  should match to ~5 kJ/mol on this benign curve.

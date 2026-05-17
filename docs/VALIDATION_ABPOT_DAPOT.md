# Cross-validation against ABPOT / DAPOT reference data

The Vyazovkin xlsx (`theory/extra/Вязовкин.xlsx`) contains α(T) tables for
two real materials at multiple heating rates, plus a column of E_a(α) values
computed by an independent (Excel-based) Vyazovkin implementation. We use
this as a regression check: feed the xlsx α(T) tables into our methods and
compare against the reference E_a column.

Reproduce:

- **Regression test (CI-runnable, no xlsx needed):** `pytest tests/test_reference_abpot_dapot.py`
- **Full per-α dump (exploratory, needs the xlsx):** `.venv/bin/python scripts/validate_abpot_dapot.py`

The regression test uses fixtures under
[`data/reference_workbooks/abpot/`](../data/reference_workbooks/abpot/) and
[`data/reference_workbooks/dapot/`](../data/reference_workbooks/dapot/), which
were extracted from the xlsx by
[`scripts/extract_reference_fixtures.py`](../scripts/extract_reference_fixtures.py).
The xlsx itself can be deleted after extraction — see
[THEORY_FILES_AUDIT.md](THEORY_FILES_AUDIT.md).

## Summary

| Method            | ABPOT mean abs Δ | ABPOT max abs Δ | DAPOT mean abs Δ | DAPOT max abs Δ |
|-------------------|-----------------:|----------------:|-----------------:|----------------:|
| KAS               |         **0.69** |          11.39  |         **1.58** |          20.38  |
| OFW               |             1.90 |          12.28  |             2.23 |          18.90  |
| **Vyazovkin**     |         **0.57** |          11.56  |         **1.55** |          20.49  |
| Friedman          |            34.20 |         143.42  |            23.03 |         216.06  |
| Vyazovkin-AIC     |            29.63 |          92.27  |            18.02 |          80.04  |

(All values in kJ/mol. ABPOT compared at 95 α-points, DAPOT at 96 α-points.
Rates: ABPOT 0.2/1/2.5/5 K/min, DAPOT 1/2.5/5/10 K/min.)

## Interpretation

**Integral methods match the independent reference to within ~1 kJ/mol mean
error.** The reference E_a column itself is a Vyazovkin computation, so our
classical Vyazovkin (mean Δ = 0.57 / 1.55 kJ/mol) is essentially a numerical
agreement check on the implementation of the temperature integral
(Senum–Yang) and the Φ(E) minimization. Both come out clean.

The max error of ~11–20 kJ/mol always falls at α ≤ 0.05 or α ≥ 0.95 —
boundary-interpolation effects on the very coarse α-grid present in the
xlsx. Inside [0.05, 0.95] all integral methods agree to a few kJ/mol.

OFW lands ~1–2 kJ/mol higher than KAS / Vyazovkin, exactly as expected from
Doyle's approximation bias documented in `docs/ALGORITHMS.md` §D.

**Friedman and Vyazovkin-AIC look bad here, and that's a property of the
input format, not a bug.** Both rely on numerical differentiation /
fine-grained integration of α(t):

- *Friedman* needs `dα/dt`, which we get via `np.gradient(α, T)·β`. On the
  coarse α=0.01-step xlsx grid this is dominated by quantization noise.
- *Vyazovkin-AIC* integrates exp(−E/RT(t)) over a window of α-half-width
  Δα = 0.02. With α-grid spacing ~0.01 the window has only 3–4 sample
  points, and the trapezoidal estimate is noisy.

Run the same methods on raw DSC heat-flow data (the path our pipeline is
designed for) and they smooth out — synthetic-data tests in
`tests/test_synthetic_recovery.py` recover E_a to within 1 kJ/mol on
identical reaction parameters. The xlsx feeds in already-processed α(T)
tables, which throws away the underlying smoothness that the differential
methods rely on.

## What this validates

- The **temperature-integral approximation** (Senum–Yang p(x)) matches the
  reference workbook.
- The **Vyazovkin Φ(E) minimization** matches the reference workbook.
- Our **KAS and OFW linear regressions** match the reference workbook.
- These four methods on pre-processed α(T) tables track the reference
  within ~1–2 kJ/mol mean error.

## What this does *not* validate

- Friedman / AIC on raw DSC data — those tests live in
  `tests/test_synthetic_recovery.py`, not here.
- The full DSC pre-processing pipeline (baseline subtraction, integration,
  normalization) — covered separately by the synthetic round-trip test.
- The unlabeled β column in the DAPOT sheet (probably 1 K/min); we used
  the three labeled rates only.

## Status

Done — fixtures live under `data/reference_workbooks/{abpot,dapot}/` and the
regression check runs as part of the normal `pytest` suite. The xlsx is
no longer required for CI / development.

## Related fixtures

The same reference-workbook directory hosts additional test materials and
xlsx cross-check data extracted in the same pass:

- [`data/reference_workbooks/epoxy_pv15/`](../data/reference_workbooks/epoxy_pv15/) —
  Epoxy PV15-0,5 (3 rates + ref E_a, non-flat E(α));
  see [`tests/test_reference_epoxy.py`](../tests/test_reference_epoxy.py).
- [`data/reference_workbooks/skm_prepreg/`](../data/reference_workbooks/skm_prepreg/) —
  SKM 6 rates (self-consistency stress test);
  see [`tests/test_reference_skm.py`](../tests/test_reference_skm.py).
- [`data/reference_workbooks/xlsx_crosscheck/`](../data/reference_workbooks/xlsx_crosscheck/) —
  precomputed Z(α), A(α), reaction-order y(n) tables used by
  `tests/test_xlsx_crosscheck_*.py` to verify our formula implementations
  bit-for-bit against the original xlsx.

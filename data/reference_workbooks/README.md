# Reference workbooks

Source: independent Excel workbooks (`theory/extra/*.xlsx`, gitignored).
This directory contains the *extracted* numerical data — no .xlsx kept
in version control.

Material data is split per folder; cross-check tables that originate
from xlsx-internal computations live in `xlsx_crosscheck/`.

| Folder | Material | Heating rates (K/min) | Has reference E_a? |
|--------|----------|-----------------------|--------------------|
| [`abpot/`](abpot/) | ABPOT polymer | 0.2, 1.0, 2.5, 5.0 | ✓ |
| [`dapot/`](dapot/) | DAPOT polymer | 1.0, 2.5, 5.0, 10.0 | ✓ |
| [`epoxy_pv15/`](epoxy_pv15/) | Epoxy PV15-0,5 | 2.5, 5.0, 10.0 | ✓ |
| [`skm_prepreg/`](skm_prepreg/) | SKM preimpregnated | 1.0, 2.5, 5.0, 10.0, 20.0, 40.0 | ✗ |
| [`xlsx_crosscheck/`](xlsx_crosscheck/) | precomputed Z, A, y(n) | n/a | n/a |

## Re-extraction

If the source xlsx files are restored under `theory/extra/`, regenerate
everything with::

    .venv/bin/python scripts/extract_reference_fixtures.py

The script is idempotent — re-running overwrites existing TSV files
with identical content (bit-exact under the same xlsx).

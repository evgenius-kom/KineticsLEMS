# KineticsLEMS

Isoconversional kinetic analysis of thermal-analysis data (DSC, TGA, FSC, POM).
Given a set of experiments at different heating rates, compute the activation
energy E_a as a function of conversion α using the canonical model-free
methods.

Implemented methods:

| Method            | Type          | What it gives             |
|-------------------|---------------|---------------------------|
| Friedman          | differential  | E_a(α)                    |
| Kissinger–Akahira–Sunose (KAS) | integral, linearized   | E_a(α) |
| Ozawa–Flynn–Wall (OFW)         | integral, Doyle approx | E_a(α) |
| Kissinger (classical)          | peak-temperature       | single E_a |
| Vyazovkin (1996/97)            | nonlinear integral, Senum–Yang | E_a(α) |
| Vyazovkin AIC (2001)           | nonlinear integral, numerical  | E_a(α) (assumes linear heating; see [docs/ALGORITHMS.md](docs/ALGORITHMS.md) §G) |

Math, derivations, and implementation choices are documented in
[docs/ALGORITHMS.md](docs/ALGORITHMS.md).

## Install

```bash
python3.12 -m venv .venv
source .venv/bin/activate
pip install -e .            # add ".[dev]" for tests
```

Requires Python ≥ 3.11 (uses `tomllib` from the standard library).

## Quick start

Generate a synthetic case with known kinetics, then analyze it:

```bash
kinetics-lems generate examples/synthetic/F1_120kJ \
    --rates 2.5 5 10 20 --ea 120 --A 1e10 --model F1

kinetics-lems analyze examples/synthetic/F1_120kJ --out out/F1_120kJ
```

Sample output:

```
Material: synthetic-F1
Method:   DSC
Type:     heating
Rates:    [2.5, 5.0, 10.0, 20.0] K/min
  friedman       mean E_a =  120.00 kJ/mol
  kas            mean E_a =  119.74 kJ/mol
  ofw            mean E_a =  121.84 kJ/mol
  vyazovkin      mean E_a =  120.00 kJ/mol
  vyazovkin_aic  mean E_a =  120.00 kJ/mol
  kissinger      E_a =  119.73 kJ/mol (R²=1.0000)
```

The output directory contains one CSV per method (E_a vs α with R²) plus
PNG plots: `Ea_vs_alpha.png` (overlay of all isoconversional methods) and
`kissinger.png` (Kissinger linear fit).

## Input format

A "case" is either a directory or a `.zip` of one. Layout:

```
case/
  settings.json
  rate_2.5.txt        # whatever filenames you reference in settings.json
  rate_5.txt
  rate_10.txt
  ...
```

`settings.json`:

```json
{
  "Settings": [
    {"name": "Material",       "value": "ABPOT"},
    {"name": "ExperimentType", "value": "heating"},
    {"name": "Method",         "value": "DSC"},
    {"name": "Conditions",     "value": {
        "rate_2.5.txt": "2.5",
        "rate_5.txt":   "5",
        "rate_10.txt":  "10",
        "rate_20.txt":  "20"
    }}
  ]
}
```

Notes:

- `Conditions` maps wave file → heating rate β (K/min) for `heating`/`cooling`,
  or → temperature (K) for `isothermal`.
- `Method` is one of `DSC`, `TGA`, `FSC`, `POM`. The current analysis pipeline
  treats them all as scalar y(T); peak interpretation differs but the math
  does not.
- Each wave file is plain text, two columns separated by spaces / tabs:
  `T (K)<TAB>y` per line. Comma decimals (`316,5`) are tolerated. Lines that
  don't parse are silently skipped.

A complete worked example lives in [`examples/synthetic/F1_120kJ/`](examples/synthetic/F1_120kJ/).

## Configuration

Algorithm parameters live in a single TOML file
([`configs/default.toml`](configs/default.toml)) — the conversion grid,
which methods to run, the E search bracket, and the AIC window width:

```toml
[conversion]
step = 0.05
min  = 0.05
max  = 0.95

[methods]
enabled = ["friedman", "kas", "ofw", "kissinger", "vyazovkin", "vyazovkin_aic"]

[methods.vyazovkin]
ea_bracket_kJ = [1.0, 600.0]

[methods.vyazovkin_aic]
delta_alpha   = 0.02
ea_bracket_kJ = [1.0, 600.0]

[output]
directory  = "out"
save_plots = true
save_csv   = true
plot_dpi   = 120
```

Override at runtime: `kinetics-lems analyze case/ --config my-config.toml`.
The case description (material, rates, file map) stays in `settings.json` —
TOML governs *how* to analyze, JSON governs *what* to analyze.

## Output

Per analysis run, in the output directory:

- `friedman.csv`, `kas.csv`, `ofw.csv`, `vyazovkin.csv`, `vyazovkin_aic.csv` —
  one row per α with `Ea_kJ_per_mol`, `intercept`, `r_squared`.
- `kissinger.csv` — peak temperatures, fitted E_a, A.
- `Ea_vs_alpha.png` — overlay of E_a(α) for every isoconversional method.
- `kissinger.png` — linear fit on Kissinger coordinates.

## Generating synthetic data

`kinetics-lems generate` integrates the Arrhenius rate equation
`dα/dt = A·f(α)·exp(-E/RT)` for given β and writes a case identical in shape
to a real experiment:

```bash
kinetics-lems generate OUT_DIR \
    --rates 5 10 20         # heating rates K/min
    --ea 150                # E_a in kJ/mol
    --A 1e12                # pre-exponential 1/s
    --model F1              # F1, F2, A2, R2, R3
    --noise 0.02            # optional relative noise on dα/dT
```

T window auto-fits the peak when `--T-start` / `--T-stop` are omitted.

## Development

```bash
pip install -e ".[dev]"
pytest                      # unit + end-to-end synthetic recovery tests
```

Repo layout:

```
src/kinetics_lems/
  __init__.py
  constants.py              # R, Doyle factor, unit conversions
  models.py                 # Wave, CaseData, CaseParams, enums
  config.py                 # TOML loading and validation
  conversion.py             # baseline subtraction, α(T), T(α), dα/dt(α)
  runner.py                 # orchestration: case + config -> AnalysisResults
  reporting.py              # CSV + PNG output
  cli.py                    # `kinetics-lems` entry point
  io/                       # case + wave file readers (folder or .zip)
  methods/                  # one file per algorithm
  synthetic/                # ground-truth data generator
configs/default.toml
docs/                       # ALGORITHMS, validation, roadmap, future-feature specs
tests/
  fixtures/                 # ABPOT/DAPOT reference α(T) tables for regression tests
scripts/                    # one-off extraction + validation utilities
examples/synthetic/         # generated by `kinetics-lems generate` (gitignored)
```

## License

MIT.

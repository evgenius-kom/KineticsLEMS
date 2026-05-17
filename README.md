# KineticsLEMS

Isoconversional kinetic analysis of thermal-analysis data (DSC, TGA, FSC, POM).
Given a set of experiments at different heating rates, compute the activation
energy E_a as a function of conversion α using the canonical model-free
methods.

Implemented methods (13):

| Method            | Type          | What it gives             |
|-------------------|---------------|---------------------------|
| Friedman          | differential  | E_a(α)                    |
| Kissinger–Akahira–Sunose (KAS) | integral, linearized   | E_a(α) |
| Ozawa–Flynn–Wall (OFW)         | integral, Doyle approx | E_a(α) |
| Kissinger (classical)          | peak-temperature       | single E_a |
| Vyazovkin (1996/97)            | nonlinear integral, Senum–Yang | E_a(α) |
| Vyazovkin AIC (2001)           | nonlinear integral, numerical  | E_a(α) (linear heating; see [docs/ALGORITHMS.md](docs/ALGORITHMS.md) §G) |
| Z(α) master plot  | model discrimination | best-fit reaction model + RMS ranking |
| Pre-exponential A | kinetic triplet | A(α) median/MAD, auto-picks f(α) from Z-plot |
| Multi-step        | E(α) segmentation | step boundaries, per-step E_a, flatness score |
| Reaction order n  | F_n linearization | best n + E_a, R²(n) sweep |
| Coats–Redfern     | per-rate model-fit | E_a + A for every g(α), independent ranking |
| Jackknife uncertainty | leave-one-run-out | mean E(α) ± 95% CI |
| Predictive lifetime | isothermal integration | α(t) and time-to-α at user temperatures |

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
Rates:    [2.5, 5.0, 10.0, 20.0] K/min
  friedman       mean E_a =  120.00 kJ/mol
  kas            mean E_a =  119.74 kJ/mol
  ofw            mean E_a =  121.84 kJ/mol
  vyazovkin      mean E_a =  120.00 kJ/mol
  vyazovkin_aic  mean E_a =  120.00 kJ/mol
  kissinger      E_a =  119.73 kJ/mol (R²=1.0000)
  master_plot    best fit: A2  [top-3: A2 (RMS=0.001), F1 (RMS=0.001), A4 (RMS=0.001)]
  preexponential under f(α)=A2: log10 A = 9.78 ± 0.15  (A ≈ 6.00e+09 1/s)
  multistep      n=1 step(s), E_a flatness 0.0%
  reaction_order best n = 0.99  (E_a = 119.5 kJ/mol, R² = 1.0000)
  coats_redfern  best model: F1  [top-3: F1 (R²=1.000), A2 (R²=1.000), A3 (R²=1.000)]
  uncertainty    vyazovkin jackknife: mean SE = 0.00 kJ/mol  (n_runs = 4)
  lifetime       isothermal time-to-α (f(α) = A2):
       T (°C)   t(α=0.05)   t(α=0.10)   t(α=0.50)   t(α=0.90)
         25.0         0.0s     559.9y    3392.4y    7233.8y
         60.0         0.0s       3.5y      21.0y      44.7y
```

The output directory contains a CSV per method plus PNG + PDF plots
(300 DPI, Type-42 fonts) for each visual artifact. See the **Output**
section below for the complete file list.

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
# Order matters: preexponential depends on vyazovkin + master_plot;
# lifetime depends on preexponential.
enabled = [
    "friedman", "kas", "ofw", "kissinger",
    "vyazovkin", "vyazovkin_aic",
    "master_plot", "preexponential", "multistep",
    "reaction_order", "coats_redfern", "uncertainty", "lifetime",
]

[methods.preexponential]
model = ""               # blank → auto-pick from master_plot

[methods.multistep]
jump_threshold   = 0.10
min_segment_size = 3

[methods.reaction_order]
alpha_min = 0.1
alpha_max = 0.9
n_min     = 0.1
n_max     = 4.0
n_steps   = 80

[methods.coats_redfern]
alpha_min = 0.1
alpha_max = 0.9
# models = ["F1", "A2", "R3"]   # optional subset; default: all 12

[methods.uncertainty]
method = "vyazovkin"     # estimator to jackknife (needs ≥ 3 runs)

[methods.lifetime]
temperatures_C = [25.0, 40.0, 60.0]
alpha_targets  = [0.05, 0.10, 0.50, 0.90]
model          = ""      # blank → use master_plot best

[methods.vyazovkin]
ea_bracket_kJ = [1.0, 600.0]

[methods.vyazovkin_aic]
delta_alpha   = 0.02
ea_bracket_kJ = [1.0, 600.0]

[output]
directory         = "out"
save_plots        = true
save_csv          = true
plot_dpi          = 300
plot_formats      = ["png", "pdf"]
per_method_panels = false
```

Override at runtime: `kinetics-lems analyze case/ --config my-config.toml`.
The case description (material, rates, file map) stays in `settings.json` —
TOML governs *how* to analyze, JSON governs *what* to analyze.

## Output

Per analysis run, in the output directory:

**CSV (one row per α unless noted):**

| File | Columns |
|------|---------|
| `friedman.csv`, `kas.csv`, `ofw.csv`, `vyazovkin.csv`, `vyazovkin_aic.csv` | `alpha`, `Ea_kJ_per_mol`, `intercept`, `r_squared` |
| `kissinger.csv` | peak temperatures + fitted E_a, A |
| `master_plot_z.csv`, `master_plot_ranking.csv` | wide Z(α) table + RMS ranking |
| `preexponential.csv` | `alpha`, `A_median_per_sec` (+ summary header) |
| `multistep.csv` | one row per detected step |
| `reaction_order.csv` | `n`, `r_squared`, `Ea_kJ_per_mol` over the swept grid |
| `coats_redfern_fits.csv`, `coats_redfern_ranking.csv` | per-(model, run) fit + mean R² ranking |
| `uncertainty.csv` | `alpha`, `Ea_mean`, `Ea_se`, `Ea_ci95_low`, `Ea_ci95_high` |
| `lifetime_times.csv`, `lifetime_curves.csv` | time-to-α table + full α(t) sweep per T |

**Plots (PNG + PDF at 300 DPI, embedded Type-42 fonts):**

- `Ea_vs_alpha` — overlay of E_a(α) for all isoconversional methods.
- `kissinger` — linear fit on Kissinger coordinates.
- `master_plot_z` — experimental Z(α) vs master curves.
- `preexponential` — A(α) profile with MAD error band.
- `multistep` — E(α) from Vyazovkin with detected step boundaries.
- `reaction_order` — R²(n) and E_a(n) over the n sweep.
- `coats_redfern` — bar chart of mean R² per model, ranked.
- `uncertainty` — E(α) with shaded 95% jackknife CI.
- `lifetime` — semilog α(t) curves, one per isothermal temperature.

## Generating synthetic data

`kinetics-lems generate` integrates the Arrhenius rate equation
`dα/dt = A·f(α)·exp(-E/RT)` for given β and writes a case identical in shape
to a real experiment:

```bash
kinetics-lems generate OUT_DIR \
    --rates 5 10 20         # heating rates K/min
    --ea 150                # E_a in kJ/mol
    --A 1e12                # pre-exponential 1/s
    --model F1              # F1, F2, F3, A2, A3, A4, R2, R3
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
  plotting.py               # paper_style() context manager, save_figure()
  runner.py                 # orchestration: case + config -> AnalysisResults
  reporting.py              # Ea vs α CSV + PNG/PDF output
  reporting_master_plot.py  # Z(α) CSV + plot
  reporting_preexp.py       # A(α) CSV + plot
  reporting_multistep.py    # multi-step CSV + plot
  reporting_reaction_order.py
  reporting_coats_redfern.py
  reporting_uncertainty.py
  reporting_lifetime.py
  cli.py                    # `kinetics-lems` entry point
  io/                       # case + wave file readers (folder or .zip)
  methods/
    friedman.py             # differential isoconversional
    kas.py                  # KAS integral isoconversional
    ofw.py                  # OFW (Doyle approx)
    kissinger.py            # peak-temperature method
    vyazovkin.py            # classical + AIC nonlinear integral
    master_plot.py          # Criado–Málek Z(α), 12-model ranking
    preexponential.py       # A from kinetic triplet, log-space median/MAD
    multistep.py            # greedy E(α) segmentation, flatness score
    reaction_order.py       # F_n linearization sweep over n
    coats_redfern.py        # per-rate model-fit baseline
    uncertainty.py          # jackknife-by-run E(α) confidence intervals
    lifetime.py             # predictive isothermal α(t) and time-to-α
  synthetic/                # ground-truth data generator
configs/default.toml
docs/                       # ALGORITHMS, validation, roadmap, future-feature specs
data/reference_workbooks/   # real-material α(T) fixtures (ABPOT, DAPOT, Epoxy, SKM)
                            # plus xlsx-formula cross-check tables
tests/                      # unit + integration + cross-validation tests
scripts/                    # one-off extraction + validation utilities
examples/synthetic/         # generated by `kinetics-lems generate` (gitignored)
```

## License

MIT.

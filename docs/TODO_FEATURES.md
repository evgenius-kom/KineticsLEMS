# Roadmap — future kinetics features

Forward-looking list of features beyond the current core. Items shipped
in code are documented in [docs/ALGORITHMS.md](ALGORITHMS.md); this file
tracks **what's next** and **why each deferred item is non-trivial**.

Status legend:
- ★ HIGH PRIORITY — strong user value, recommended next
- ◐ MEDIUM PRIORITY — useful but specialized
- ○ LOW PRIORITY — nice-to-have / niche

---

## ★ HIGH PRIORITY

### 1. Implement the model-based global fitter (`fitting/`)

**Why.** ICTAC 2020 §2 frames isoconversional analysis as *screening*,
not as a final model. The full kinetic triplet recovery for multi-step
processes (parallel / consecutive / competitive) needs a global fit
across all heating rates. The scaffold is in place
([fitting/](../src/kinetics_lems/fitting/)) — only the SINGLE topology
RHS, an `Objective.evaluate` implementation, and a thin `differential_evolution +
least_squares` driver are missing.

**Scope.**

1. SINGLE-topology ODE in `fitting/topology.py` (uses
   `MASTER_MODELS` for f(α), simple ``scipy.integrate.solve_ivp`` driver).
2. `Objective.evaluate` in `fitting/objective.py` (formula already in
   the module docstring).
3. Optimizer adapter using `scipy.optimize.differential_evolution → least_squares`.
4. PARALLEL after the 2-parallel synthetic dataset
   (already available: [synthetic/multistep.py](../src/kinetics_lems/synthetic/multistep.py)).

**Effort.** ≈ 600–800 LOC + 6–10 new tests.

**Validation.** Recover (E, A) of the F1 synthetic and the 2-parallel
synthetic to within 2–5%.

### 2. DAEM (Gaussian + discrete-grid)

**Why.** Distributed reactivity is the right model for biomass /
heterogeneous polymers / kerogens. Gaussian-DAEM synthetic data is
already shipped in [synthetic/multistep.py](../src/kinetics_lems/synthetic/multistep.py);
the fit side is missing.

**Scope.**

* Discrete-energy DAEM (free non-negative weights on a grid).
* Gaussian-parametric DAEM with shared A.
* Smoothness penalty option.

**Effort.** ≈ 250 LOC + 4 tests, ``scipy.optimize.nnls`` for discrete,
``scipy.optimize.least_squares`` for Gaussian.

**Validation.** Recover ``(Ē, σ)`` for `generate_daem_gaussian_case`.

### 3. Hold-out / leave-one-rate-out validation

**Why.** Predicting the *unfitted* β is the only rigorous out-of-sample
check for kinetics; ICTAC 2020 §6 makes it mandatory in serious reports.

**Scope.** Wrap the future fitting engine: refit on all but one rate,
predict the held-out α(T), compute RMSE; emit a separate table in the
Markdown report.

**Effort.** ≈ 80 LOC after the fitter exists.

---

## ◐ MEDIUM PRIORITY

### 4. MDA / KDA peak deconvolution

Fit `dα/dT` as a sum of asymmetric peak functions (Fraser-Suzuki,
Weibull) when E(α) shows shoulders or multiple maxima. ≈ 300 LOC; needs
only `scipy.optimize`. Warn loudly that MDA components are *not*
mechanistic.

### 5. Promote `schemas/` out of stub status

Migrate I/O boundaries (vendor adapters → `CaseData` → schema) to the
pydantic types in [schemas/canonical.py](../src/kinetics_lems/schemas/canonical.py).
This unlocks JSON-stable export of every artefact (fit results, model
specifications) and aligns with §5 of
`kinetics_research_implementation_notes.md`.

### 6. Implement at least one real vendor adapter

Pick whichever lab format actually arrives most often (NETZSCH Proteus
ASCII is the strongest candidate). The scaffold in
[io/vendors/](../src/kinetics_lems/io/vendors/) makes this drop-in:
only `load()` needs to be filled.

### 7. Coats–Redfern with explicit per-rate A reporting

Already 95% there — the per-run A is stored on `CoatsRedfernRunFit`,
just expose it in CSV / plot more prominently.

### 8. Cantera YAML export adapter

Emit the fitted kinetic triplet as a Cantera-compatible YAML reaction
mechanism so users can plug our results into reactor-network simulations.
≈ 100 LOC + dependency on `cantera` (optional).

---

## ○ LOW PRIORITY / NICE-TO-HAVE

### 9. LaTeX/BibTeX export

Generate `report.tex` + `references.bib` from `AnalysisResults`
containing every result table and method-specific citation. Targets
academic users writing papers.

### 10. ABPOT 2016 finer-grid replacement fixture

`Сводная таблица.xlsx` sheet `Новая кинетика` contains newer ABPOT
data with a denser α grid than the currently-extracted fixture. Worth
re-extracting only if the current ABPOT regression test starts showing
precision-limited failures.

### 11. Pr–Tm single-material Z(α) fixture

`Zalfa plots Pr-Tm program.xlsm` contains single-material Z(α)
calculations. Superseded by the cross-check fixtures already in
[`data/reference_workbooks/xlsx_crosscheck/`](../data/reference_workbooks/xlsx_crosscheck/);
only worth migrating if Pr-Tm is needed as a separate test material.

### 12. Bayesian / MCMC posterior over (E, A)

`pymc` or `emcee` on the SINGLE-topology fitter would give credible
intervals instead of jackknife SE. ICTAC 2020 §6 recommends *not* doing
this before deterministic identifiability is solid (Burnham & Anderson
2002, §6.3) — defer until the fitter is mature.

---

## §X. External libraries / tools — when each becomes worth the dependency cost

These are author notes on §4 of `kinetics_research_implementation_notes.md`.
Decision lens: **what marginal user problem does adding this dependency
solve?** Everything below is *opt-in* — none should join the core
``dependencies`` list.

### Already in the project

| Library | Why it's in core | Notes |
|---|---|---|
| `numpy`, `scipy` | Numerical primitives, ODE solver, optimisers, signal processing | Unavoidable. |
| `matplotlib` | Plot generation for reports | Headless-friendly; PDF + PNG output. |
| `pydantic` ≥ 2.5 | Canonical schemas at I/O boundaries | Added 2026-05; minimal runtime overhead. |

### Worth adding when we cross specific thresholds

| Library | Trigger to add | Caveats |
|---|---|---|
| `pint` | Once any vendor adapter writes non-SI units into a schema field. Currently every adapter normalises to K/s/J inside `case_loader`, so `pint` is overhead with no win. | Beware: serialising `pint.Quantity` through pydantic v2 needs a custom validator. |
| `cvxpy` | DAEM with smoothness regularisation on the energy grid. `scipy.optimize.nnls` covers the unconstrained discrete case. | Pulls in solver binaries; might fail on minimal CI. |
| `lmfit` | Never as core — keep as an *optional* adapter for users who want named-parameter UI on top of `scipy.optimize`. | Hides identifiability problems behind nice reports. |
| `numba` | Profiling shows the F1 inner loop dominates a real run. Currently the existing implementations finish in seconds even on ABPOT-sized data. | JIT compile cost on first call; CI cache matters. |
| `pandas` | If reports start needing pivot tables we cannot do with built-in `csv`. | Today pure-`csv` works for every output. |
| `pyarrow` | If we start writing Parquet for downstream analytics. | Not needed today. |

### Worth adding when the project scope expands

| Library | Triggered by which scope change | Notes |
|---|---|---|
| `pyPESTO` + `AMICI` + `PEtab` | Promoting `fitting/` from SINGLE-topology to general ODE mechanisms with profile likelihoods. | Heavy; build-from-source on macOS sometimes problematic. |
| `cantera` | If we need gas-phase combustion / reactor-network mechanisms (notes §3.3). User's A=solid-state-only decision says **no for now**. | Add only as an export-target adapter. |
| `rdkit` | If we ever take a SMILES / structure as input to predict kinetic priors. | Out of scope while A=solid-state. |
| `jax` + `diffrax` + `equinox` | When SciPy optimisation becomes the bottleneck for global multi-step fits. | Pin Python and CUDA carefully. |

### Will NOT be added under current scope (A = solid-state)

| Library / tool | Why excluded today |
|---|---|
| Cantera as a *runtime* engine | We don't need gas-phase mechanism integration; YAML export-only is the right level. |
| RMG-Py / Arkane / AutoTST | Mechanism generation belongs to combustion / atmospheric chemistry, not our solid-state thermal-analysis niche. |
| COPASI / SBML | Biochemical-network territory; A=solid-state excluded. |
| CatPred / UniKP / DLKcat / DLTKcat / CataPro | Enzyme kinetic predictors; A=solid-state excluded. |
| CatMAP / PyMKM / OpenMKM / PyCatKin | Heterogeneous-catalysis microkinetics; out of scope. |
| kmcos / kmos | Lattice kinetic Monte Carlo; out of scope. |

### Vendor binary parsers

Per `kinetics_research_implementation_notes.md` §3.2/5.1, **raw
proprietary binary formats** (e.g. NETZSCH `.ngb`, TA `.tri`) are
intentionally **out of scope**. Every vendor's ecosystem exports
ASCII/CSV/JSON — we cover the export side via
[io/vendors/](../src/kinetics_lems/io/vendors/).

### ML / surrogate models

Per [ml/](../src/kinetics_lems/ml/) docstring: ML must never become a
hard dependency. The contract for any future enzyme-predictor /
baseline-assistant / surrogate-model plugin is already in
[ml/plugin.py](../src/kinetics_lems/ml/plugin.py) — implementations join
as optional subpackages, each with their own footprint.

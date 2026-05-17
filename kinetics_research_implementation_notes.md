# Kinetics Analysis Software: Current Methods, Tools, Formats, and Implementation Plan

Last updated: 2026-05-17

Scope: computational analysis of physical/chemical kinetics from thermal analysis, calorimetry, reactor/time-series data, combustion mechanisms, catalysis/microkinetics, and enzyme kinetics. Emphasis: what to use now, what to implement internally, what to avoid, and what Claude Code should build first.

---

## 0. Bottom line

For TGA/DSC/STA kinetic analysis, the current defensible workflow is not “AI replaces kinetics”. The robust workflow is:

1. Import raw/processed data with metadata.
2. Convert signal to conversion `alpha(T,t)` and rate `dalpha/dt` with explicit baseline handling.
3. Run model-free/isoconversional screening, especially `E_alpha(alpha)`.
4. Decide whether the process is plausibly single-step, multi-step, distributed-reactivity, diffusion-controlled, or deconvolvable.
5. Fit constrained kinetic models globally across all experiments, not one curve at a time.
6. Validate by prediction on temperature programs not used in fitting.
7. Report uncertainty, identifiability, parameter correlations, and extrapolation limits.

For enzyme kinetics, open AI models are already useful as priors/screening tools. For thermal kinetics, ML should be treated as surrogate modeling, anomaly detection, baseline/peak assistance, or hypothesis generation — not as a replacement for `E`, `A`, `f(alpha)`, and validation.

---

## 1. Core kinetic analysis principles

### 1.1 Thermal-analysis kinetics: current reference workflow

The ICTAC multi-step recommendations remain the key baseline for TGA/DSC thermal kinetics. They identify four main approaches for multi-step kinetic analysis:

- Multi-step model-fitting.
- Distributed reactivity / distributed activation energy models.
- Isoconversional analysis.
- Deconvolution analysis.

The target is the full kinetic triplet, not only activation energy:

```text
kinetic triplet = E, A, f(alpha)
```

where:

- `E` = activation energy.
- `A` = pre-exponential factor.
- `f(alpha)` = reaction model / conversion dependence.

Key operating rules:

- Do not infer single-step behavior from one smooth DSC/DTG peak. Multi-step kinetics can look like a single peak.
- Use multiple heating rates or multiple temperature programs. For non-isothermal data, the slowest and fastest heating/cooling rates should preferably differ by roughly one order of magnitude.
- For screening, compute `E_alpha(alpha)` over `alpha = 0.1–0.9`. Variability larger than the typical uncertainty band, roughly >10–20% of average `E_alpha`, is a strong warning that the process is not single-step.
- Prefer Friedman or flexible/advanced Vyazovkin isoconversional methods over rigid integral methods when `E_alpha` varies strongly.
- Multi-step fitting must be physically constrained: fewest meaningful steps, realistic parameter bounds, and validation outside the fitting conditions.

Primary source: ICTAC 2020 multi-step recommendations, Thermochimica Acta, DOI `10.1016/j.tca.2020.178597`.

### 1.2 Practical decision tree

```text
Input: multiple TGA/DSC/STA/MS/FTIR experiments

1. Normalize and preprocess
   -> alpha, dalpha/dt, T(t), beta(t), metadata

2. Visual diagnostics
   -> peaks, shoulders, mass-gain/loss, exo/endo overlap, baseline quality

3. Isoconversional diagnostics
   -> E_alpha(alpha), uncertainty, alpha range reliability

4. Classify
   a. E_alpha almost constant, one peak, no metadata conflict
      -> candidate single-step model
   b. E_alpha varies strongly; shoulders/peaks or rate-dependent total effect
      -> multi-step model-fitting
   c. Heterogeneous material: biomass, coal, char, solid fuel, broad site distribution
      -> DAEM / distributed reactivity
   d. Strongly overlapping peaks with interpretable components
      -> MDA/KDA deconvolution, then kinetic fitting
   e. Curing/crosslinking with late slowdown/vitrification
      -> Kamal-Sourour + diffusion/vitrification correction
   f. Crystallization/nucleation
      -> Avrami/JMAK, Hoffman-Lauritzen, Turnbull-Fisher variants

5. Fit globally across all experiments

6. Validate
   -> hold-out heating rate, isothermal run, arbitrary T(t), residual diagnostics

7. Report
   -> parameter estimates, confidence/profile intervals, correlations, prediction limits
```

---

## 2. Methods to support

### 2.1 Model-free / isoconversional methods

Use for initial diagnosis and activation-energy profiles.

| Method | Data | Output | Priority | Notes |
|---|---|---:|---:|---|
| Friedman differential method | `alpha`, `dalpha/dt`, `T` at several programs | `E_alpha(alpha)` | P0 | Sensitive to derivative noise; strong diagnostic power. |
| Vyazovkin nonlinear integral | non-isothermal/isothermal data | `E_alpha(alpha)` | P0 | Better for arbitrary temperature programs; implement with robust numerical integration. |
| Advanced Vyazovkin | variable `E_alpha` | `E_alpha(alpha)` | P1 | Best practice when `E_alpha` changes with conversion. |
| KAS | constant heating rates | approximate `E_alpha(alpha)` | P1 | Useful benchmark; not enough alone. |
| OFW/Flynn-Wall-Ozawa | constant heating rates | approximate `E_alpha(alpha)` | P2 | Include for compatibility/literature comparison. |
| Kissinger peak method | peak temperatures at several `beta` | one effective `E` | P2 | Use only as a quick diagnostic; not a full kinetic model. |

Implementation notes:

- `alpha` grid: default `0.05–0.95`, report primary range `0.1–0.9`.
- Always expose uncertainty via bootstrap/resampling.
- Warn on unreliable endpoints: low signal-to-noise near `alpha < 0.1` and `alpha > 0.9`.
- Treat `E_alpha` as a diagnostic plus possible constraint for later model fitting, not as the final model.

### 2.2 Single-step reaction models

Support as building blocks, not as blind model choices.

| Model family | Typical use | Formula sketch |
|---|---|---|
| nth-order `Fn` | decomposition, simple apparent kinetics | `f(alpha) = (1-alpha)^n` |
| Contracting geometry `Rn` | phase-boundary processes | geometry-dependent |
| Diffusion models `Dn` | solid-state diffusion limitation | model-specific |
| Avrami/JMAK `An` | nucleation and growth, crystallization | `alpha = 1 - exp(-(kt)^n)` in isothermal form |
| Prout-Tompkins / autocatalytic | sigmoid/autocatalytic processes | `f(alpha)=alpha^m(1-alpha)^n` |
| Sestak-Berggren | empirical flexible model | `f(alpha)=alpha^m(1-alpha)^n[-ln(1-alpha)]^p` |
| Kamal-Sourour | thermoset curing/crosslinking | `dalpha/dt=(k1+k2 alpha^m)(1-alpha)^n` |

Rules:

- If several models fit similarly, prefer the model with fewer parameters and stronger physical support.
- Do not interpret Sestak-Berggren parameters mechanistically unless independent evidence supports it.
- Include physical bounds: `E > 0` for most chemical reactions, `A > 0`, `0 <= alpha <= 1`, weights sum to 1, typical reaction orders not arbitrarily huge.

### 2.3 Multi-step model-fitting

Use when `E_alpha(alpha)` varies significantly, the signal shows multiple peaks/shoulders, or chemistry implies parallel/consecutive/competitive steps.

Model classes to implement:

```text
Independent parallel:
A_i -> P_i
rate_total = sum_i w_i rate_i

Competitive:
A -> B
A -> C
shared reactant pool

Consecutive:
A -> B -> C
intermediate species balance required

Mixed:
A -> B
A -> C -> D
parallel + consecutive branches
```

Fit targets:

- `alpha(T,t)`.
- `dalpha/dt(T,t)`.
- Optionally simultaneous mass + heat-flow + evolved gas channels.

Objective function should support:

```text
loss = w_alpha * normalized_SSE(alpha)
     + w_rate  * normalized_SSE(dalpha/dt)
     + w_peak  * peak-position/shape penalties, optional
     + regularization/constraints
```

Minimum implementation requirements:

- Simultaneous global fit across all experiments.
- Multi-start optimization.
- Bounds and parameter transformations: fit `logA`, not `A` directly.
- Constraints: `w_i >= 0`, `sum(w_i)=1`, realistic order/exponent bounds.
- Residual plots per experiment and per signal type.
- Fit reproducibility: fixed seeds, saved initial guesses, saved optimizer settings.

### 2.4 Distributed activation energy / DAEM

Use for heterogeneous materials where a continuum of reactivities is more physically honest than many pseudo-discrete steps.

High-value use cases:

- Biomass pyrolysis.
- Coal/kerogen/oil shale pyrolysis.
- Polymer char devolatilization.
- Desorption from heterogeneous sites.
- Solid fuels and complex organic residues.

Models to support:

| Distribution | Priority | Notes |
|---|---:|---|
| Discrete activation energy distribution | P1 | Most flexible and practical; optimize weights over energy grid. |
| Gaussian DAEM | P1 | Common baseline; symmetric limitation. |
| Weibull / Gamma | P2 | Better asymmetry; harder optimization. |
| Logistic / asymmetric logistic | P3 | Useful research option; less critical. |

Implementation details:

- Energy grid spacing: start with ~2–5 kJ/mol; expose as setting.
- Default reaction order: first-order channels; allow nth-order channels later.
- Common `A` default; optional compensation relation `ln(A_i)=a+bE_i`.
- Fit both conversion and rate when data quality allows.
- Add sparsity/regularization to avoid noisy artificial distributions.

### 2.5 Deconvolution: MDA and KDA

Use when overlapping peaks can plausibly be separated into individual kinetic events.

MDA = mathematical deconvolution analysis:

- Fit rate peak as a sum of asymmetric peak functions.
- Candidate functions: Fraser-Suzuki, Weibull, asymmetric Gaussian/logistic.
- Use only as preprocessing or hypothesis generation.

KDA = kinetic deconvolution analysis:

- Fit total rate directly as a sum of kinetic step rates.
- Prefer KDA over purely mathematical deconvolution when physical models are available.

Warnings:

- Mathematical peak components are not automatically reaction steps.
- Opposite-sign contributions, such as exothermic/endothermic overlap or mass loss/gain overlap, require special handling.
- Deconvolution without independent evidence is high-risk.

### 2.6 Mechanistic ODE fitting

Use outside thermal-analysis-only problems, or when concentration/speciation data exist.

Data:

- Concentration time series.
- Reactor outlet composition.
- Spectroscopy-derived species trajectories.
- Pressure/temperature/flow histories.
- Calorimetry with known enthalpy model.

Recommended formalism:

```text
dC/dt = S * r(C, T, theta)
r_j = k_j(T) * product_i C_i^nu_ij
k_j(T) = A_j * T^b_j * exp(-E_j / RT)
```

Use tools:

- `pyPESTO + AMICI + PEtab` for serious parameter estimation and uncertainty.
- `Cantera` for gas-phase/combustion/reactor networks.
- `COPASI` for biochemical networks and SBML workflows.

---

## 3. Tools and libraries to know/use

### 3.1 Open-source thermal kinetics tools

| Tool | Language | License/status | Best use | Link |
|---|---|---|---|---|
| pICNIK | Python | Open-source; GitHub | Isoconversional computations: OFW, KAS, Friedman, Vyazovkin, advanced Vyazovkin | https://github.com/ErickErock/pICNIK |
| pkynetics | Python | Open-source; GitHub | DSC/TGA/dilatometry preprocessing, Friedman/KAS/OFW, Kissinger, JMAK | https://github.com/PPeitsch/pkynetics |
| KinOpt | Python | BSD-3-Clause | Isoconversional analysis + model optimization | https://github.com/alan-tabore/KinOpt |
| THINKS | Mixed/app suite | Free/open GitHub project; check redistribution license before embedding | Thermal kinetic extraction; isoconversional and model-fitting workflows | https://github.com/muravyevV/thinks |
| Kinetic Calculation | Python? | MIT in SoftwareX repo | TGA kinetic triplets, Vyazovkin + OFW | https://github.com/ElsevierSoftwareX/SOFTX_2019_190 |
| Thermokinetics-Toolset | MATLAB | Open GitHub scripts | Semi-automated TGA/DSC kinetic triplets | https://github.com/Xenozite/Thermokinetics-Toolset |
| tga_data_analysis | Python | MIT | Automated TGA preprocessing/comparability, not full kinetics | https://github.com/mpecchi/tga_data_analysis |

Recommendation: do not depend on one of these as the core engine. Use them as references, benchmarks, and test-case sources. Build your own core with explicit data models, unit handling, multi-experiment fitting, and validation.

### 3.2 Commercial thermal kinetics software used in labs

| Software | Vendor | Current role | What to support on our side |
|---|---|---|---|
| NETZSCH Kinetics Neo | NETZSCH | Industrial/research TGA/DSC/DEA/ARC/MS kinetic analysis, model-free/model-based, predictions | Import/export ASCII/TXT/CSV; support their common 3-column data style; support arbitrary external temperature profile `T,t` CSV. |
| NETZSCH Proteus | NETZSCH | Instrument control/evaluation for NETZSCH thermal instruments | Support exported ASCII/CSV rather than raw proprietary files. |
| NETZSCH Termica Neo | NETZSCH | Newer thermal kinetics workflow using ICTAC-type methods; compatible with Kinetics Neo | Same as above; treat as external comparison source. |
| AKTS Thermokinetics TK | AKTS | Advanced thermokinetics, DSC/TG/TG-MS/TG-FTIR/DSF, safety/storage/lifetime prediction | Support ASCII/CSV import/export and reproduce core model-free/model-fitting outputs. |
| TA Instruments TRIOS | TA/Waters | Instrument control + analysis for DSC/TGA/SDT/rheology/microcalorimetry | Support TRIOS JSON export and CSV/text exports. Do not rely on raw proprietary files. |
| Mettler Toledo STARe + Model Free Kinetics | Mettler Toledo | Instrument ecosystem; TGA/DSC prediction under different conditions | Support ASCII exports from STARe; expect user-defined export protocols. |
| PerkinElmer Pyris / Pyris Player | PerkinElmer | DSC/TGA instrument ecosystem | Support CSV/TXT exports. |
| Shimadzu TA software | Shimadzu | DSC/TGA instrument ecosystem | Support CSV/TXT exports. |
| Setaram Calisto + AKTS integration | Setaram/Kep | Thermal analysis + AKTS workflow | Support ASCII/CSV exports. |

Most practical compatibility target: **generic ASCII/CSV + TRIOS JSON + well-documented canonical JSON/YAML**. Raw proprietary binary formats should be explicitly out of scope unless a stable public parser exists.

### 3.3 General chemical kinetics and mechanism tools

| Tool | Type | Use | Link |
|---|---|---|---|
| Cantera | Open-source kinetics/thermo/transport engine | Gas-phase/liquid kinetics, reactor networks, flames, sensitivity, YAML/CHEMKIN-style mechanisms | https://cantera.org/ and https://github.com/Cantera/cantera |
| RMG-Py | Automatic reaction mechanism generation | Pyrolysis, combustion, atmospheric chemistry; elementary-step mechanisms | https://github.com/ReactionMechanismGenerator/RMG-Py |
| Arkane | RMG subpackage | Thermochemistry, high-pressure-limit rates, pressure-dependent rates from quantum chemistry | https://reactionmechanismgenerator.github.io/RMG-Py/users/arkane/introduction.html |
| AutoTST | Automated transition-state theory | Automated TST for selected reaction families | https://github.com/ReactionMechanismGenerator/AutoTST |
| T3 | RMG/ARC orchestration | Iterative, sensitivity-driven kinetic model development | https://github.com/ReactionMechanismGenerator/T3 |
| COPASI | Biochemical network simulator | SBML, ODE/SDE/Gillespie, parameter estimation | https://copasi.org/ |
| pyPESTO | Parameter estimation framework | Multi-start/global optimization, uncertainty, AMICI/PEtab integration | https://github.com/icb-dcm/pypesto |
| AMICI | ODE sensitivity/simulation | Efficient model simulation/sensitivities, SBML | https://github.com/AMICI-dev/AMICI |
| PEtab | Model/data parameter-estimation format | Reproducible parameter estimation problems | https://petab.readthedocs.io/ |
| PySINDy | Sparse model discovery | Candidate ODE discovery from time-series data | https://github.com/dynamicslab/pysindy |

### 3.4 Catalysis, microkinetics, and kinetic Monte Carlo

| Tool | Type | Use | Link |
|---|---|---|---|
| CatMAP | Descriptor-based microkinetics | Catalyst screening, descriptor maps, reaction-condition studies | https://github.com/SUNCAT-Center/catmap |
| PyMKM | Microkinetic modeling | Surface coverages, rates, selectivity, apparent activation energy, reaction orders | https://github.com/LopezGroup-ICIQ/pymkm |
| OpenMKM | C++/Cantera-based microkinetics | Multiscale heterogeneous catalytic reaction modeling | https://github.com/VlachosGroup/openmkm |
| PyCatKin | Microkinetics + energy span | Heterogeneous catalysis modeling | https://github.com/aab64/PyCatKin |
| kmcos | Lattice kMC | Kinetic Monte Carlo for lattice chemical systems | https://github.com/kmcos/kmcos |
| kmos | Lattice kMC | Surface science / heterogeneous catalysis kMC | https://github.com/mhoffman/kmos |

---

## 4. AI / ML models with practical relevance

### 4.1 Enzyme kinetic parameter models

These are the most directly useful open AI models for kinetic parameters.

| Model/tool | Input | Output | Status | Use in our pipeline |
|---|---|---|---|---|
| CatPred | Enzyme sequence/structure + substrate; supports batch workflows | `kcat`, `Km`, `Ki`, `kcat/Km` depending mode | Open GitHub + web app | Use as external predictor/prior for enzyme module; store uncertainty and training-domain flags. |
| UniKP / EF-UniKP | Protein sequence + substrate structure; EF-UniKP includes pH/T factors | `kcat`, `Km`, `kcat/Km` | GitHub research code | Use as enzyme screening/prior; compare against CatPred. |
| DLKcat | Protein sequence + substrate SMILES | `kcat` | GitHub, widely cited | Use with caution; include similarity/domain checks due published critiques. |
| DLTKcat | Sequence + substrate + temperature | temperature-dependent `kcat` | GPL-3 GitHub | Useful for bioprocess temperature dependence; license matters if embedding. |
| CataPro | Protein/small molecule language models + fingerprints | kinetic parameters | GitHub research code | Evaluate as newer baseline for enzyme discovery. |
| OpenKinetics Predictor | Web predictor | `kcat`, `Km`, `kcat/Km` | Web access | Useful for manual benchmarking, not core dependency. |

Important caveat: enzyme models predict **enzyme-catalyzed kinetic constants**, not TGA/DSC thermal decomposition parameters.

### 4.2 AI for reaction rate constants and mechanism generation

Useful directions:

- RMG/Arkane/ARC/AutoTST ecosystem for automated mechanism and rate estimation.
- ML rate-rule estimators such as structurally interpolating decision trees in RMG-type workflows.
- Graph neural networks / transformer models for reaction property prediction are active research, but not yet a universal production-grade replacement for mechanism generation plus validation.

Practical use in our software:

- Use AI/ML predictions as priors, initial guesses, or candidate mechanism suggestions.
- Always tag AI-derived parameters as `source=predicted`, with model name/version, input hash, domain/similarity diagnostics, and uncertainty if available.
- Never silently mix predicted and experimentally fitted parameters.

### 4.3 AI for thermal-analysis curves

Current role:

- Baseline suggestion.
- Peak/shoulder detection.
- Curve classification.
- Fast surrogate prediction once a validated dataset exists.
- Initial guess generation for model fitting.

Not acceptable as standalone:

- “ANN predicted TGA curve” as a mechanistic kinetic model.
- Inferring `E, A, f(alpha)` without multi-rate validation.
- Storage/safety/lifetime prediction outside the training domain.

---

## 5. File formats and output compatibility

### 5.1 Formats likely encountered from instruments/software

| Source | Likely accessible export | What to implement |
|---|---|---|
| NETZSCH Kinetics Neo | ASCII/TXT/CSV; supported data types include DSC, TGA, DIL, DEA, ARC, MS; external temperature profile ASCII/CSV | Generic delimited importer with column mapping; specialized Kinetics Neo template. |
| NETZSCH Proteus | ASCII/CSV export from proprietary measurement files | Generic importer; preserve instrument metadata when present. |
| TA Instruments TRIOS | JSON export; CSV/text exports | Native TRIOS JSON parser; CSV fallback. |
| Mettler Toledo STARe | ASCII export via export protocol | Generic delimited importer; metadata sidecar support. |
| AKTS TK | Vendor-neutral thermoanalytical data; ASCII/CSV workflows | Import/export ASCII/CSV; allow comparison against AKTS reports. |
| Older/other thermal software | TXT, CSV, DAT, XLSX | Robust table sniffing + manual column assignment. |
| Custom labs | Excel/CSV with arbitrary headers | Interactive mapping + saved parser profiles. |

### 5.2 Canonical internal data model

Use a strict schema. Do not let vendor-specific quirks leak into numerical code.

```yaml
Experiment:
  id: string
  source_file: string
  source_vendor: optional string
  instrument: optional string
  technique: enum[TGA, DSC, DTA, STA, DEA, ARC, MS, FTIR, DIL, TMA, CUSTOM]
  sample:
    name: string
    mass_initial_mg: optional float
    mass_final_mg: optional float
    geometry: optional string
    particle_size: optional string
  atmosphere:
    gas: optional string
    flow_ml_min: optional float
    pressure: optional float
  program:
    mode: enum[dynamic, isothermal, step, arbitrary]
    heating_rate_K_min: optional float
    time_s: array[float]
    temperature_K: array[float]
  signals:
    - name: string
      unit: string
      values: array[float]
      role: enum[mass, mass_loss, heat_flow, dsc, dtg, conversion, rate, ms_channel, ftir_band, custom]
  preprocessing:
    baseline: optional object
    smoothing: optional object
    conversion_definition: optional object
```

### 5.3 Canonical outputs to support

#### `fit_result.json`

```json
{
  "model_id": "multi_step_parallel_v1",
  "software_version": "0.1.0",
  "data_hashes": ["..."],
  "fit_scope": {
    "experiments_used": ["exp_001", "exp_002", "exp_003"],
    "holdout_experiments": ["exp_004"]
  },
  "parameters": [
    {"name": "E_1", "value": 125.4, "unit": "kJ/mol", "stderr": 4.2, "bounds": [0, 400]},
    {"name": "logA_1", "value": 18.2, "unit": "ln(1/s)", "stderr": 0.8},
    {"name": "w_1", "value": 0.62, "unit": "dimensionless"}
  ],
  "parameter_correlation": [[1.0, 0.92], [0.92, 1.0]],
  "objective": {
    "loss_total": 0.0142,
    "loss_alpha": 0.0081,
    "loss_rate": 0.0061
  },
  "diagnostics": {
    "aic": 123.4,
    "bic": 140.2,
    "rss": 0.012,
    "warnings": ["High E-logA correlation in step 1"]
  }
}
```

#### `model.yaml`

```yaml
model_type: multi_step
topology: parallel
steps:
  - id: step1
    reaction_model: Fn
    parameters:
      E:
        value: 125.4
        unit: kJ/mol
      logA:
        value: 18.2
        unit: ln(1/s)
      n:
        value: 1.0
      weight:
        value: 0.62
  - id: step2
    reaction_model: Fn
    parameters:
      E:
        value: 178.0
        unit: kJ/mol
      logA:
        value: 25.7
        unit: ln(1/s)
      n:
        value: 1.2
      weight:
        value: 0.38
validity:
  temperature_range_K: [420, 760]
  heating_rates_K_min: [2, 20]
  atmosphere: nitrogen
```

#### `prediction.csv`

```text
time_s,temperature_K,alpha,dalpha_dt,signal_predicted,experiment_id,model_id
0,300.0,0.0000,0.0,0.0,sim_001,multi_step_parallel_v1
...
```

---

## 6. Recommended implementation plan for Claude Code

### P0 — Build the data foundation first

Deliverables:

1. Python package structure:

```text
kinetics_lems/
  io/
  models/
  preprocessing/
  isoconversional/
  fitting/
  prediction/
  diagnostics/
  reports/
  schemas/
  tests/
```

2. Strict schemas with `pydantic`:

- `Experiment`.
- `Signal`.
- `TemperatureProgram`.
- `ConversionCurve`.
- `KineticModel`.
- `FitResult`.
- `PredictionResult`.

3. Importers:

- Generic CSV/TXT/ASCII with delimiter/decimal/header detection.
- TA TRIOS JSON parser.
- NETZSCH/Kinetics Neo-style 2–3 column ASCII parser.
- Mettler STARe ASCII parser profile.
- Manual column mapping saved as YAML.

4. Unit handling with `pint`:

- Temperature: C/K.
- Time: s/min/h.
- Heating rate: K/min, K/s.
- Energy: J/mol, kJ/mol.
- Mass: mg/g/%.
- Heat flow: W/mW/mW mg^-1.

Acceptance criteria:

- Can import at least 5 exported example files with different separators and decimal formats.
- Converts all temperatures to Kelvin internally.
- Stores original columns and preprocessing provenance.
- Fails loudly on ambiguous columns.

Recommended libraries:

```text
pandas, numpy, pydantic, pint, pyarrow, openpyxl, scipy
```

### P1 — Preprocessing and conversion engine

Deliverables:

1. Baseline correction:

- Manual endpoints.
- Linear baseline.
- Piecewise baseline.
- Optional asymmetric least-squares baseline for signals where appropriate.

2. Smoothing:

- Savitzky-Golay with full provenance.
- No smoothing by default for fitting unless user explicitly enables it.

3. Conversion calculation:

For TGA:

```text
alpha = (m0 - m_t) / (m0 - m_inf)
```

For DSC:

```text
alpha = partial_area / total_area
```

4. Rate calculation:

- Numerical derivative from `alpha(t)`.
- Direct DTG/DSC-derived rate when available.
- Rate uncertainty estimate if replicate data exist.

Acceptance criteria:

- Every conversion has a documented definition and selected baseline range.
- Endpoint instability warning if `alpha` normalization is sensitive to baseline choice.
- User can regenerate results from saved preprocessing config.

Recommended libraries:

```text
scipy.signal, numpy, pandas, pydantic
```

### P2 — Isoconversional module

Deliverables:

1. Implement:

- Friedman.
- KAS.
- OFW.
- Vyazovkin nonlinear integral.
- Advanced/incremental Vyazovkin.

2. Diagnostics:

- `E_alpha(alpha)` plot.
- Confidence intervals via bootstrap.
- Endpoint reliability flags.
- Single-step/multi-step warning rule.

3. Export:

- `e_alpha.csv`.
- `isoconversional_result.json`.
- plots as SVG/PNG.

Acceptance criteria:

- Benchmarks against pICNIK on known examples.
- Unit tests for synthetic first-order data recover true `E` within tolerance.
- Handles arbitrary non-isothermal temperature profile for Vyazovkin methods.

Recommended libraries:

```text
numpy, scipy.integrate, scipy.optimize, joblib, matplotlib or plotly
```

### P3 — Single-step and multi-step fitting engine

Deliverables:

1. Reaction model registry:

```python
reaction_models = {
    "Fn": f_fn,
    "An": f_avrami_like,
    "PT": f_prout_tompkins,
    "SB": f_sestak_berggren,
    "Kamal": f_kamal_sourour,
    "DiffusionLimited": f_diffusion_limited,
}
```

2. ODE solver for:

- Single-step.
- Parallel independent steps.
- Consecutive steps.
- Competitive steps.
- Mixed topologies.

3. Fit API:

```python
fit_model(
    experiments=[...],
    model_spec=..., 
    objective={"alpha": 1.0, "rate": 1.0},
    optimizer="differential_evolution_then_least_squares",
    bounds=...,
    constraints=...,
)
```

4. Optimizers:

- `scipy.optimize.differential_evolution` for global search.
- `scipy.optimize.least_squares` for local refinement.
- Multi-start local optimization.
- Optional `pyPESTO` adapter for advanced optimization/UQ.

5. Fit diagnostics:

- residual plots.
- parameter covariance/correlation.
- AIC/BIC.
- profile likelihood or bootstrap intervals.
- hold-out prediction report.

Acceptance criteria:

- Synthetic datasets: recover parameters for 1-step, 2-parallel, 2-consecutive cases.
- Reject/flag overparameterized models with non-identifiable parameters.
- All fits save optimizer seed, bounds, initial guesses, and loss components.

Recommended libraries:

```text
scipy.optimize, scipy.integrate, numpy, pandas, joblib, numba, pydantic
```

Do not start with full Bayesian inference. Add it later after deterministic fitting is stable.

### P4 — DAEM / distributed reactivity

Deliverables:

1. Discrete DAEM first:

```text
alpha_total = sum_i w_i alpha_i(E_i, A, f_i)
```

2. Supported distributions:

- Discrete grid with nonnegative weights.
- Gaussian.
- Weibull/Gamma later.

3. Regularization:

- Nonnegative weights.
- Sum weights = 1.
- Smoothness penalty or sparsity penalty configurable.

4. Outputs:

- `D(E)` table.
- predicted rates/conversion.
- uncertainty bands if bootstrapped.

Acceptance criteria:

- Fits benchmark biomass/coal-like synthetic data.
- Does not create unstable spiky distributions by default.
- Supports fitting conversion and rate simultaneously.

Recommended libraries:

```text
numpy, scipy.optimize, scipy.sparse, cvxpy optional, numba optional
```

### P5 — Deconvolution module

Deliverables:

1. MDA peak fitting:

- Fraser-Suzuki.
- Weibull peak.
- asymmetric Gaussian/logistic optional.

2. KDA step decomposition:

- Use kinetic model registry.
- Support positive and negative contributions.

3. Warnings:

- Component count sensitivity.
- Nonphysical negative/huge parameters.
- Lack of independent mechanistic support.

Acceptance criteria:

- Can reproduce synthetic overlapped peaks.
- Clearly labels MDA components as mathematical unless KDA/mechanistic validation exists.

Recommended libraries:

```text
scipy.optimize, lmfit optional, numpy
```

### P6 — Prediction/simulation engine

Deliverables:

1. Predict under:

- Constant heating rate.
- Isothermal hold.
- Step programs.
- Arbitrary external `T(t)` from CSV.

2. Outputs:

- `alpha(t,T)`.
- `dalpha/dt`.
- predicted mass/heat-flow if conversion-to-signal mapping known.
- time-to-conversion and temperature-to-conversion.

3. Validity warnings:

- Extrapolation outside fitted temperature range.
- Extrapolation outside fitted heating-rate range.
- Phase transition/vitrification flags if model includes them.

Acceptance criteria:

- Can run Kinetics-Neo-like external temperature profile prediction from ASCII/CSV.
- Report clearly distinguishes interpolation vs extrapolation.

### P7 — Reporting and interoperability

Deliverables:

1. Markdown/HTML/PDF report generator.
2. Export artifacts:

- `model.yaml`.
- `fit_result.json`.
- `e_alpha.csv`.
- `prediction.csv`.
- `report.md`.

3. Optional import/export:

- PEtab for mechanistic ODE models.
- Cantera YAML adapter for gas-phase mechanisms.
- SBML adapter for biochemical networks.

Acceptance criteria:

- A complete run is reproducible from a project folder and one config file.
- Reports include all preprocessing and fitting settings.

### P8 — AI/ML layer

Implement only after P0–P7 are stable.

Features:

1. Baseline/peak assistant:

- Suggest baseline regions.
- Suggest number of peaks/shoulders.
- Never silently apply changes.

2. Initial-guess generator:

- Use isoconversional outputs and peak features to propose model candidates.

3. Surrogate model:

- Train only on validated fitted models or high-quality replicated datasets.
- Output prediction uncertainty and domain distance.

4. Enzyme predictor integration:

- Wrap CatPred/UniKP/DLKcat/DLTKcat as optional plugins.
- Store model version, input hash, sequence identity/domain similarity when possible.

Do not make ML a dependency for core thermal kinetics.

---

## 7. Specific libraries recommended for our own implementation

### 7.1 Core stack

| Need | Recommended library | Why |
|---|---|---|
| Arrays/numerics | `numpy` | Standard, stable. |
| Data tables | `pandas`, `pyarrow` | CSV/Excel/Parquet workflows. |
| Schemas | `pydantic` | Strict validation and JSON export. |
| Units | `pint` | Avoid C/K, min/s, J/kJ mistakes. |
| ODE/integration | `scipy.integrate.solve_ivp` | Good first choice. |
| Optimization | `scipy.optimize` | Global + local optimizers. |
| Speed | `numba` | Useful for repeated model evaluation. |
| Parallelism | `joblib`, `concurrent.futures` | Multi-start and bootstrap. |
| Plots | `matplotlib`, `plotly` | Static report + interactive UI. |
| Testing | `pytest`, `hypothesis` | Synthetic data and property tests. |
| CLI | `typer` | Clean command-line interface. |
| Config | `ruamel.yaml` or `pyyaml` | YAML project configs. |

### 7.2 Optional advanced stack

| Need | Library | When to use |
|---|---|---|
| Advanced parameter estimation | `pyPESTO` | Multi-start/global fitting, uncertainty, black-box objectives. |
| ODE sensitivities | `AMICI` | Large mechanistic ODEs, SBML, sensitivity-based fitting. |
| Reproducible estimation specs | `PEtab` | If supporting systems biology / mechanistic ODE exchange. |
| Differentiable/JIT modeling | `jax`, `diffrax`, `equinox` | Later, if SciPy/Numba becomes bottleneck. |
| Convex subproblems | `cvxpy` | DAEM nonnegative/smooth weight fitting. |
| Reaction mechanisms | `cantera` | Combustion/gas-phase/reactor networks. |
| Cheminformatics | `rdkit` | Molecule parsing, substrate features, reaction metadata. |
| Sparse ODE discovery | `pysindy` | Hypothesis generation from time series. |

### 7.3 Libraries/tools to avoid as core dependencies

| Avoid as core | Why |
|---|---|
| Raw vendor binary parsers | Fragile, often undocumented/licensed; use exported ASCII/CSV/JSON. |
| Excel as internal format | Good import/export only; bad for reproducibility and units. |
| `lmfit` as the only optimizer | Convenient but can hide constraints/identifiability issues; okay as optional wrapper. |
| Neural networks for first release | Not needed for scientifically defensible core; high validation burden. |
| Single-method kinetic packages as hard dependency | Too narrow; use for benchmarking only. |
| Bayesian MCMC before deterministic stability | Expensive and misleading if model/identifiability is bad. |

---

## 8. What to avoid scientifically

| Avoid | Why | Better alternative |
|---|---|---|
| Single heating-rate kinetic parameters | Non-identifiable, often meaningless. | Multiple heating rates/programs + isoconversional screening. |
| Treating Kissinger `E` as complete kinetics | One effective peak parameter, not full mechanism. | Use `E_alpha(alpha)` and model fitting. |
| Blind Coats-Redfern model selection | Can produce false precision and wrong mechanisms. | ICTAC-style model-free + constrained model-based workflow. |
| Overfitting many pseudo-steps | Low RSS does not imply physical meaning. | Fewest physically justified steps; F-test/AIC/BIC + validation. |
| Reporting only `R^2` | High `R^2` is easy for smooth curves. | Residuals, hold-out prediction, parameter uncertainty. |
| Ignoring `E-logA` compensation/correlation | Parameters may be non-identifiable. | Correlation matrix, profile likelihood, constrained fitting. |
| Heavy smoothing before derivatives | Derivative kinetics can be distorted. | Store raw + processed, test sensitivity to smoothing. |
| Extrapolating high-rate DSC/TGA to low-temperature storage blindly | Mechanism/phase/diffusion can change. | Add low-rate/isothermal validation; flag extrapolation. |
| Interpreting mathematical peak deconvolution as mechanism | Peak functions are descriptive. | Use KDA + independent MS/FTIR/XRD/Raman/microscopy evidence. |
| Treating ML thermal predictions as mechanistic kinetics | ML surrogate may fail outside domain. | ML only with domain checks and uncertainty. |

---

## 9. Minimal product specification

### 9.1 CLI commands

```bash
kinetics import data/*.csv --profile generic_tga --out project/
kinetics preprocess project/config.yaml
kinetics iso project/config.yaml --methods friedman,vyazovkin,kas
kinetics fit project/config.yaml --model models/two_parallel_fn.yaml
kinetics predict project/config.yaml --temperature-profile profiles/storage_profile.csv
kinetics report project/config.yaml --out report.md
```

### 9.2 Python API sketch

```python
from kinetics_lems.io import load_experiments
from kinetics_lems.preprocessing import compute_conversion
from kinetics_lems.isoconversional import friedman, vyazovkin
from kinetics_lems.models import MultiStepModel, Fn
from kinetics_lems.fitting import fit_model
from kinetics_lems.prediction import simulate_temperature_program

experiments = load_experiments("project/raw", profile="generic_tga")
curves = [compute_conversion(exp, method="tga_mass_loss") for exp in experiments]
iso = vyazovkin(curves, alpha_grid="0.05:0.95:0.01")

model = MultiStepModel.parallel([
    Fn(id="step1"),
    Fn(id="step2"),
])

result = fit_model(
    curves,
    model,
    objective={"alpha": 1.0, "rate": 1.0},
    optimizer="de_then_ls",
    multistart=100,
)

prediction = simulate_temperature_program(result.model, "profiles/profile.csv")
```

### 9.3 UI priorities

1. Data import wizard with column mapping and unit confirmation.
2. Baseline/conversion editor with side-by-side raw and processed signal.
3. Isoconversional dashboard: `E_alpha(alpha)`, warnings, method comparison.
4. Model builder: single/parallel/consecutive/competitive topology.
5. Fit dashboard: residuals, parameters, correlations, hold-out prediction.
6. Prediction tab: arbitrary `T(t)` file, isothermal, constant-rate.
7. Export/report tab.

---

## 10. Suggested repository issues for Claude Code

### Issue 1 — Create canonical schemas and import pipeline

**Goal:** create validated data objects and generic ASCII/CSV import.

**Tasks:**

- Add `Experiment`, `Signal`, `TemperatureProgram` pydantic models.
- Add CSV/TXT sniffer.
- Add unit normalization via `pint`.
- Add parser profiles in YAML.
- Add tests with synthetic CSVs.

**Acceptance:** `pytest` passes; `kinetics import` produces canonical `experiment.json` and normalized Parquet/CSV.

### Issue 2 — Implement conversion and preprocessing

**Goal:** convert TGA/DSC signals to `alpha` and `dalpha/dt`.

**Tasks:**

- TGA mass-loss conversion.
- DSC integrated-area conversion.
- Baseline selection config.
- Savitzky-Golay optional smoothing.
- Store preprocessing provenance.

**Acceptance:** synthetic first-order mass-loss data converts to monotonic `alpha` in `[0,1]` and stable derivative.

### Issue 3 — Implement Friedman/KAS/OFW

**Goal:** first model-free diagnostics.

**Tasks:**

- Alpha-grid interpolation across experiments.
- Friedman linear regression per alpha.
- KAS/OFW regressions for constant-rate data.
- CI via bootstrap.
- Plots and CSV/JSON export.

**Acceptance:** synthetic Arrhenius first-order data recovers `E` within 5–10% depending noise level.

### Issue 4 — Implement Vyazovkin nonlinear isoconversional method

**Goal:** support high-quality isoconversional analysis.

**Tasks:**

- Temperature integral numerical implementation.
- Arbitrary `T(t)` support.
- Incremental/advanced variant.
- Benchmark against pICNIK where possible.

**Acceptance:** first-order synthetic and variable-`E_alpha` synthetic tests pass.

### Issue 5 — Implement kinetic model registry

**Goal:** reusable reaction models.

**Tasks:**

- Implement `Fn`, `An`, Prout-Tompkins, Sestak-Berggren, Kamal-Sourour.
- Define parameter bounds and defaults.
- Add model YAML serialization.

**Acceptance:** unit tests verify finite rates and expected qualitative peak shapes.

### Issue 6 — Implement global fitting engine

**Goal:** fit single-step and multi-step models across experiments.

**Tasks:**

- ODE system generator for single, parallel, consecutive, competitive models.
- Objective on `alpha` and `rate`.
- Differential evolution + least-squares refinement.
- Multi-start.
- Parameter correlation and AIC/BIC.

**Acceptance:** recover known parameters from synthetic two-step data; warnings for non-identifiable overparameterized model.

### Issue 7 — Implement prediction engine

**Goal:** simulate fitted models under new temperature programs.

**Tasks:**

- Constant heating rate.
- Isothermal.
- Arbitrary `T(t)` CSV.
- Extrapolation warnings.

**Acceptance:** can reproduce fitted experiments and generate predictions for external profile.

### Issue 8 — Implement report generation

**Goal:** one-command reproducible report.

**Tasks:**

- Markdown report template.
- Embed plots.
- Include preprocessing, methods, parameters, warnings, validation.
- Export model/result files.

**Acceptance:** `kinetics report` creates a self-contained project report.

### Issue 9 — Implement DAEM module

**Goal:** support distributed reactivity.

**Tasks:**

- Discrete energy grid.
- Nonnegative weights and sum-to-one constraint.
- Common `A`; optional compensation relation later.
- Smoothness regularization.

**Acceptance:** synthetic distributed-energy data recovered without spurious oscillatory weights.

### Issue 10 — Add optional enzyme AI plugin

**Goal:** integrate open enzyme kinetic predictors as optional tools.

**Tasks:**

- Define predictor interface.
- Add CatPred wrapper if installable/reproducible.
- Add UniKP/DLKcat wrappers as optional plugins.
- Store prediction provenance and domain warnings.

**Acceptance:** plugin returns predictions with model version and input hash; core thermal package does not depend on ML packages.

---

## 11. Validation datasets to create internally

Create synthetic datasets before relying on real lab exports.

| Dataset | Purpose |
|---|---|
| `synthetic_1step_fn_dynamic` | Validate Friedman/KAS/OFW/Vyazovkin and single-step fitting. |
| `synthetic_2parallel_fn_dynamic` | Validate multi-step fitting and `E_alpha(alpha)` variability. |
| `synthetic_2consecutive_fn_dynamic` | Validate species-balance ODE topology. |
| `synthetic_kamal_curing` | Validate autocatalytic curing and DSC conversion. |
| `synthetic_daem_gaussian` | Validate DAEM. |
| `synthetic_overlapped_peaks` | Validate MDA/KDA. |
| `synthetic_arbitrary_temperature_profile` | Validate prediction under non-linear `T(t)`. |

Then add real/exported examples from:

- NETZSCH/Proteus ASCII/CSV.
- TA TRIOS JSON.
- Mettler STARe ASCII.
- Public examples from pICNIK/pkynetics/KinOpt where licenses permit.

---

## 12. Reference links

### Standards/recommendations and reviews

- ICTAC multi-step kinetics recommendations, Thermochimica Acta, 2020, DOI: `10.1016/j.tca.2020.178597`.
- ICTAC thermal decomposition kinetics recommendations, 2022: https://www.osti.gov/biblio/1964012
- Pyrolysis kinetics quality review, 2022: https://www.mdpi.com/2673-7264/2/4/29
- TGA + ANN biomass pyrolysis review, 2025: https://pubs.acs.org/doi/10.1021/acsomega.5c05098

### Thermal kinetics tools

- pICNIK: https://github.com/ErickErock/pICNIK
- pkynetics: https://github.com/PPeitsch/pkynetics
- KinOpt: https://github.com/alan-tabore/KinOpt
- THINKS: https://github.com/muravyevV/thinks and https://chemphys.space/thinks/
- Kinetic Calculation: https://github.com/ElsevierSoftwareX/SOFTX_2019_190
- tga_data_analysis: https://github.com/mpecchi/tga_data_analysis
- NETZSCH Kinetics Neo import docs: https://kinetics.netzsch.com/en/learn/import-of-user-defined-data
- NETZSCH Kinetics Neo technical datasheet: https://kinetics.netzsch.com/_Resources/Persistent/f/2/f/0/f2f0ca1432943a23a00415c7c177a7eb16f874e7/Technical_Datasheet__Kinetics_Neo.pdf
- AKTS Thermokinetics: https://www.akts.com/tk/thermokinetics-software-thermal-analysis-isoconversional-model-fitting-dsc-tg-detailed-description/
- TA TRIOS JSON export: https://www.tainstruments.com/applications-notes/trios-software-and-json-export-overview-of-json-export-from-trios-software-and-import-into-python-tb105/
- Mettler Toledo Model Free Kinetics: https://www.mt.com/gb/en/home/products/Laboratory_Analytics_Browse/TA_Family_Browse/TA_software_browse/STARe_Software_Option_Model_Free_Kinetics_1.html

### General kinetics/mechanisms

- Cantera: https://cantera.org/ and https://github.com/Cantera/cantera
- RMG-Py: https://github.com/ReactionMechanismGenerator/RMG-Py
- RMG website: https://rmg.mit.edu/
- Arkane docs: https://reactionmechanismgenerator.github.io/RMG-Py/users/arkane/introduction.html
- AutoTST: https://github.com/ReactionMechanismGenerator/AutoTST
- T3: https://github.com/ReactionMechanismGenerator/T3
- COPASI: https://copasi.org/
- pyPESTO: https://github.com/icb-dcm/pypesto
- AMICI: https://github.com/AMICI-dev/AMICI
- PEtab: https://petab.readthedocs.io/
- PySINDy: https://github.com/dynamicslab/pysindy

### Catalysis/microkinetics/kMC

- CatMAP: https://github.com/SUNCAT-Center/catmap
- PyMKM: https://github.com/LopezGroup-ICIQ/pymkm
- OpenMKM: https://github.com/VlachosGroup/openmkm
- PyCatKin: https://github.com/aab64/PyCatKin
- kmcos: https://github.com/kmcos/kmcos
- kmos: https://github.com/mhoffman/kmos

### Enzyme AI kinetic predictors

- CatPred paper: https://www.nature.com/articles/s41467-025-57215-9
- CatPred GitHub: https://github.com/maranasgroup/CatPred
- CatPred database/code: https://github.com/maranasgroup/CatPred-DB
- UniKP: https://github.com/Luo-SynBioLab/UniKP
- UniKP paper: https://pmc.ncbi.nlm.nih.gov/articles/PMC10713628/
- DLKcat paper: https://www.nature.com/articles/s41929-022-00798-z
- DLKcat GitHub: https://github.com/SysBioChalmers/DLKcat
- DLKcat critique: https://academic.oup.com/biomethods/article-abstract/9/1/bpae061/7740676
- DLTKcat GitHub: https://github.com/SizheQiu/DLTKcat
- CataPro: https://github.com/zchwang/CataPro
- OpenKinetics Predictor: https://predictor.openkinetics.org/

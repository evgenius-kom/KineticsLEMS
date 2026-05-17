# Roadmap — future kinetics features

Forward-looking list of features beyond the current 13-method core.
Items already shipped are documented in [docs/ALGORITHMS.md](ALGORITHMS.md)
(sections H–N) and not repeated here. Priorities reflect value-to-user
weighed against implementation cost.

Status legend:
- ★ HIGH — strong user value, recommended next
- ◐ MEDIUM — useful but specialized
- ○ LOW — nice-to-have / niche

---

## ★ HIGH PRIORITY

### 1. Nonlinear T(t) support (modulated DSC / temperature jumps)

**Why.** AIC already handles non-linear heating mathematically, but the
input pipeline assumes constant β. Modern thermal-analysis data
(modulated DSC, T-jump, fast-cycling FSC) needs raw `(t, T)` samples.

**Scope.**

1. Extend `settings.json` to allow a wave file with three columns
   `t, T, y` instead of two-column `T, y`.
2. Update [`io/wave_reader.py`](../src/kinetics_lems/io/) to autodetect
   2- vs 3-column input.
3. Update [`conversion.py`](../src/kinetics_lems/conversion.py): when t
   is provided, compute dα/dt by numerical differentiation of α(t)
   instead of `(y/total)·β`.
4. Keep backwards compatibility: 2-column files still use the
   linear-heating shortcut.

**Effort.** ≈ 200 LOC + 3–5 new tests. Schema change in
`settings.json` (additive, non-breaking).

### 2. Distributed Activation Energy Model (DAEM)

**Why.** Standard for biomass pyrolysis, fossil fuels, complex polymers
where activation energy is *continuously distributed* rather than from
N discrete steps. Often gives much better fits than multi-step
deconvolution in those domains.

**Scope.**

* Gaussian distribution of E around E̅ with stddev σ; fit (E̅, σ, A)
  to the experimental dα/dt curves across all rates.
* Nonlinear least squares (`scipy.optimize.least_squares`) on pooled
  (T, dα/dt, β) data.
* Report (E̅, σ, A, R²) and the implied f(α).

**Effort.** ≈ 250 LOC. Well-defined math; reference implementations
exist in academic Python.

**Reference.** Cai et al., *Renewable Sustainable Energy Rev.* 36 (2014) 236.

### 3. ICTAC consistency check

**Why.** ICTAC 2011/2020 recommend running ≥ 2 isoconversional methods
and reporting agreement. The CLI prints them side-by-side but does not
quantify divergence.

**Scope.**

* Compute pairwise relative differences between E_a(α) curves of every
  pair of isoconversional methods.
* If max-pair-diff > threshold (default 10 %), emit a structured warning
  with the worst α and the methods involved.
* Add `consistency.csv` and the warning to CLI stdout.

**Effort.** ≈ 80 LOC + 2 tests. Pure post-processing.

---

## ◐ MEDIUM PRIORITY

### 4. Sestak–Berggren empirical model

**Why.** Three-parameter f(α) = α^m · (1−α)^n · [−ln(1−α)]^p that
interpolates between F_n, A_m, R_n, D_n. Useful when none of the 12
canonical models matches Z(α) well.

**Scope.** Add to `MASTER_MODELS` as a parametric fit; expose (m, n, p)
fitting via `scipy.optimize.curve_fit`. Allow toggle in
`[methods.master_plot] use_sestak_berggren = true`.

**Reference.** Šesták & Berggren, *Thermochim. Acta* 3 (1971) 1.

### 5. Autocatalytic models (Prout–Tompkins, Kamal–Sourour)

**Why.** Industry-standard for epoxy resin cure kinetics and other
self-catalysed reactions where dα/dt ≠ 0 at α = 0.

**Scope.** Add two new models:
* Prout–Tompkins: f(α) = α^m · (1−α)^n.
* Kamal–Sourour: dα/dt = (k₁ + k₂·α^m) · (1−α)^n  (two Arrhenius pairs).

The second is a model-fitting target (not just an f(α)); needs its own
solver since it has two competing reaction channels.

**Effort.** ≈ 200 LOC + 4 tests.

### 6. Coats–Redfern with explicit per-rate A reporting

**Why.** Current Coats–Redfern aggregates A across runs. Some users want
the per-rate A_i to assess Arrhenius compensation directly.

**Scope.** Already 90 % there — `CoatsRedfernRunFit` carries per-run A;
expose it through the CSV/plot more prominently.

**Effort.** ≈ 20 LOC + 1 test.

### 7. Savitzky–Golay smoothing for Friedman

**Why.** Differentiation amplifies noise; SavGol smooths α(T) and dα/dt
before linear regression, dramatically improving Friedman robustness on
noisy DSC traces.

**Scope.** Optional `[methods.friedman] smooth_window`, `smooth_poly`
parameters. Pre-process `dalpha_dt` with `scipy.signal.savgol_filter`
when enabled.

**Effort.** ≈ 30 LOC + 2 tests.

---

## ○ LOW PRIORITY / NICE-TO-HAVE

### 8. Compensation effect — ln A vs E plot

Plotting ln A vs E_a across α (or across Coats–Redfern models) often
reveals a linear "compensation" relationship. Diagnostic of whether a
kinetic triplet is physically meaningful or just a fitting artefact.
No new computation — just a new plot consumer of existing results.

### 9. Model-free isothermal prediction (no f(α) assumption)

Vyazovkin (2000) showed that `t(α; T_iso)` can be computed from E(α)
without choosing f(α) explicitly, by using the ratio of temperature
integrals. More robust than picking the wrong model in [`lifetime.py`](../src/kinetics_lems/methods/lifetime.py).
Add as an alternate path that takes only E(α).

### 10. LaTeX/BibTeX export

Generate `report.tex` + `references.bib` from `AnalysisResults`
containing every result table and method-specific citation. Targets
academic users writing papers.

### 11. ABPOT 2016 finer-grid replacement fixture

`Сводная таблица.xlsx` sheet `Новая кинетика` contains newer ABPOT
data with a denser α grid than the currently-extracted fixture.
Worth re-extracting only if the current ABPOT regression test starts
showing precision-limited failures — at present the existing
[`data/reference_workbooks/abpot/`](../data/reference_workbooks/abpot/) data is adequate.

### 12. Pr–Tm single-material Z(α) fixture

`Zalfa plots Pr-Tm program.xlsm` contains single-material Z(α)
calculations. Superseded by the cross-check fixtures already in
[`data/reference_workbooks/xlsx_crosscheck/`](../data/reference_workbooks/xlsx_crosscheck/);
only worth migrating if Pr-Tm is needed as a separate test material.

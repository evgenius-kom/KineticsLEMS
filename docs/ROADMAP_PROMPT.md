# Prompt for a research-focused chat — KineticsLEMS roadmap review

> **Note (2026-05).** Most of the original questions in this prompt have
> been answered and shipped — see [docs/ALGORITHMS.md](ALGORITHMS.md) for
> the implemented set and [docs/TODO_FEATURES.md](TODO_FEATURES.md) for
> what's still pending. The latest research-chat output is preserved
> verbatim in [`kinetics_research_implementation_notes.md`](../kinetics_research_implementation_notes.md)
> (repo root). Use the template below when you want a *new* review pass
> against the current state.

Copy the block below into a fresh chat with a research-strong model
(Claude Opus / GPT-5 / etc.) to get a status check and roadmap.

---

I maintain a small Python project, **KineticsLEMS**, for thermo-kinetic
analysis of solid-state reactions from DSC/TGA/FSC/POM data.

**Currently implemented:** see [docs/ALGORITHMS.md](ALGORITHMS.md) — it
is the canonical reference for every method, with formulas and source
citations. The README's first table also lists the methods at a glance.
**Infrastructure scaffolds in progress** (stubs only — see the package
`__init__.py` docstrings for status): `fitting/`, `schemas/`,
`io/vendors/`, `ml/`.

**Inputs:** zip or folder with `settings.json` + per-rate `.txt` waves —
either 2-column `T y` or 3-column `t T y` (arbitrary T(t) programs).
**Outputs:** CSV per method + matplotlib plots + a single
self-contained `report.md`.

I want a focused, no-fluff review on four questions. Be specific, include
references where you can, and rank suggestions by **value / implementation
cost**.

### 1. State of the art — what's still missing?

What modern isoconversional methods or refinements are commonly used
today that I do NOT yet implement?
- Cai–Chen, Starink, Tang and other improved integral approximations vs. KAS / Doyle.
- Nonparametric kinetics (Sempere & Nomen, NPK).
- ASTM E698 (peak-shift method — practitioner standard).
- Modulated DSC / temperature-modulated kinetics — anything actionable beyond passing
  the 3-column data through AIC?
- Anything from the ICTAC 2020 recommendations beyond ICTAC 2011 that I
  should add?

### 2. Adjacent physical-chemistry analysis I could plug into the same pipeline

Given that the input is "y(T) at multiple β" (or "y(t,T)" for arbitrary
programs) and the output is Eₐ(α), what *other* thermal-analysis or
related techniques can be processed with the same isoconversional
machinery, possibly with small adapters?
- TMA (thermomechanical analysis), DMA (dynamic mechanical analysis)?
- Conductivity / dielectric relaxation under T-ramp?
- IR / Raman peak-area kinetics (chemometric maps over T)?
- Mass spectrometry (TG-MS) coupled curves?
- Pyrolysis Py-GC/MS — single rate vs multiple rates?
- Curing kinetics / glass-transition kinetics — what specializations apply?

For each: what changes in the math vs my current pipeline, what changes in
data ingestion, and is it worth doing?

### 3. Model-based fitting — feasibility and complexity

The scaffold ([`src/kinetics_lems/fitting/`](../src/kinetics_lems/fitting/))
declares SINGLE / PARALLEL / CONSECUTIVE / COMPETITIVE / MIXED topologies
but currently `build_ode_system` raises `NotImplementedError` for all of
them. I want to know:
- Recommended SciPy-only path: `differential_evolution → least_squares`
  on a SINGLE-topology RHS first, then PARALLEL. Pitfalls?
- For multi-step reactions: is multi-stage fitting feasible in pure
  NumPy/SciPy or does it need a proper library (pyPESTO + AMICI,
  Burnham's KIN-SOL)?
- Distributed reactivity models (DAEM): how complex is a baseline
  Gaussian DAEM fit on top of my existing isoconversional pipeline?
- Profile-likelihood vs jackknife vs MCMC for uncertainty — when does
  each pay off?

### 4. Architecture & data-quality side issues

- Baseline subtraction is currently linear. Is sigmoidal / spline
  baseline a meaningful improvement for DSC, or does it just trade
  apparent smoothness for tunable bias?
- I don't yet handle truly partial peaks (where reaction is incomplete
  in the recorded T window). What's the standard practitioner remedy?
- For the 3-column (t, T, y) path: the Vyazovkin classical / KAS / OFW
  still depend on the nominal β from `settings.json`. Is there a clean
  way to make every isoconversional method *fully* T(t)-aware without
  duplicating implementations?
- Vendor data — do you know a canonical (open) sample dataset for
  NETZSCH Proteus, TA TRIOS JSON, or Mettler STARe that we could use to
  validate the adapters before they touch real lab files?

Format your answer as four numbered sections matching the questions
above, each ≤ 8 bullets, each bullet ≤ 2 sentences. Skip preambles.

# Prompt for a research-focused chat — KineticsLEMS roadmap review

Copy the block below into a fresh chat with a research-strong model
(Claude Opus / GPT-5 / etc.) to get a status check and roadmap.

---

I maintain a small Python project, **KineticsLEMS**, for thermo-kinetic
analysis of solid-state reactions from DSC/TGA/FSC/POM data. Current scope:

**Implemented (model-free, isoconversional):**
- Friedman (differential)
- KAS (Kissinger–Akahira–Sunose, integral, linearized)
- OFW (Ozawa–Flynn–Wall, integral, Doyle approximation)
- Kissinger (classical peak-temperature, single Eₐ)
- Vyazovkin (1996/1997, classical nonlinear integral with Senum–Yang p(x))
- Vyazovkin AIC (advanced isoconversional, 2001) — currently with linear-heating reconstruction of t

**Inputs:** zip or folder with `settings.json` + per-rate `.txt` waves
(T vs y). **Outputs:** Eₐ(α) per method as CSV + matplotlib plots.

I want a focused, no-fluff review on four questions. Be specific, include
references where you can, and rank suggestions by **value / implementation
cost**.

### 1. State of the art — what's missing?

What modern isoconversional methods or refinements are commonly used today
that I do NOT yet implement? In particular:
- Cai–Chen, Starink, Tang and other improved integral approximations vs. KAS / Doyle.
- Nonparametric kinetics (Sempere & Nomen, NPK).
- ASTM E698 (peak-shift method — practitioner standard).
- Modulated DSC / temperature-modulated kinetics — anything actionable?
- Master-plot / model-discrimination workflows (Criado, Málek).
- Anything from the ICTAC 2020 recommendations beyond ICTAC 2011 that I
  should add?

### 2. Adjacent physical-chemistry analysis I could plug into the same pipeline

Given that the input is "y(T) at multiple β" and the output is Eₐ(α), what
*other* thermal-analysis or related techniques can be processed with the
same isoconversional machinery, possibly with small adapters?
- TMA (thermomechanical analysis), DMA (dynamic mechanical analysis)?
- Conductivity / dielectric relaxation under T-ramp?
- IR / Raman peak-area kinetics (chemometric maps over T)?
- Mass spectrometry (TG-MS) coupled curves?
- Pyrolysis Py-GC/MS — single rate vs multiple rates?
- Curing kinetics / glass-transition kinetics — what specializations apply?

For each: what changes in the math vs my current pipeline, what changes in
data ingestion, and is it worth doing?

### 3. Model-based fitting — feasibility and complexity

Right now I'm purely model-free. I want to know:
- What's the simplest, modern model-based workflow that complements
  isoconversional analysis? (e.g. fit f(α) and A after Eₐ(α) is known, then
  compare predicted DSC to data.)
- For multi-step reactions: is multi-stage fitting (Vyazovkin's "predicting
  from kinetic triplet") feasible in pure NumPy/SciPy or does it need a
  proper library (e.g., kissinger, Burnham's KIN-SOL, etc.)?
- Distributed reactivity models (DAEM): how complex is a baseline DAEM fit?
- Master-curve ranking (Criado y(α), z(α) plots) for f(α) identification —
  how much code? Worth it?
- Provide pseudo-code or library pointers for each.

### 4. Architecture & data-quality side issues

- Baseline subtraction is currently linear (between first/last points).
  Is sigmoidal / spline baseline a meaningful improvement for DSC?
- I don't yet handle truly partial peaks (where reaction is incomplete in
  the recorded T window). What's the standard practitioner remedy?
- For Vyazovkin AIC under truly non-linear heating (modulated DSC,
  isothermal jumps), I would need to pass measured t. Is there a clean
  data-format convention for that, or should I invent one?

Format your answer as four numbered sections matching the questions above,
each ≤ 8 bullets, each bullet ≤ 2 sentences. Skip preambles.

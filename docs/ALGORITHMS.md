# Kinetic algorithms — reference

This file is a self-contained crib sheet for every algorithm implemented in
KineticsLEMS. It is written so that the local `theory/` folder (papers, xlsm,
docx — gitignored) can be deleted without losing know-how; original sources
are listed at the bottom.

Notation:

- α — degree of conversion (extent of reaction), 0 ≤ α ≤ 1.
- T — absolute temperature, K.
- t — time, s.
- β = dT/dt — heating rate. **Always converted from K/min to K/s before use.**
- E_a, A — Arrhenius activation energy (J/mol) and pre-exponential factor (1/s).
- R = 8.314 462 618 J·mol⁻¹·K⁻¹.
- f(α), g(α) — differential and integral reaction-model functions.

The general non-isothermal rate equation:

    dα/dt = A · f(α) · exp(-E_a / (R·T))                           (1)

Under linear heating T(t) = T₀ + β·t this becomes

    dα/dT = (A/β) · f(α) · exp(-E_a / (R·T)).                      (2)

The integral form:

    g(α) = ∫₀^α dα'/f(α')  =  (A/β) · ∫_{T₀}^{T} exp(-E_a/(R·T')) dT'
                            ≈ (A·E_a / (β·R)) · p(x),  x = E_a/(R·T)        (3)

Senum–Yang rational approximation of p(x) (accurate for x > 20):

    p(x) = (exp(-x) / x) · (x² + 10x + 18) / (x³ + 12x² + 36x + 24).        (4)

---

## A. Pre-processing of a DSC / TGA wave

For each heating rate β, raw "y vs T" curves (DSC heat-flow or TGA mass-loss)
are turned into α(T):

1. **Baseline subtraction.** Linear baseline drawn between the first and last
   points of the peak; subtracted pointwise. Real DSC software offers nonlinear
   baselines (sigmoid, polynomial); a linear baseline is sufficient for clean
   data and matches the protocol described in the local algorithms note.
2. **Cumulative integration.** I(T) = ∫_{T_start}^{T} y(T') dT' (trapezoid).
3. **Normalization.** α(T) = I(T) / I(T_end).
4. **dα/dt.** Under linear heating, dα/dt = (y / I(T_end)) · β (with β in K/s).
5. **Sampling at fixed α.** For a user-supplied conversion grid {α_k}, look up
   T_α and (dα/dt)_α by linear interpolation on the monotone branch of α(T).

Implementation: [`conversion.py`](../src/kinetics_lems/conversion.py).

**Why the formula in step 4 is unit-agnostic.** The same
``dα/dt = (y / total) · β`` works whether the input ``y`` is

* raw DSC heat flow Q̇ (mW): ``∫y dT = β · ΔH``, so ``y/total · β = y/ΔH``;
* dα/dT (already normalized per K): ``∫y dT = 1``, so ``y/total · β = y · β``;
* dα/dt directly: ``∫y dT = β``, so ``y/total · β = y``.

In every case the proportionality constant cancels in the ratio
``y / total``, leaving the correct rate. So you don't need to know the
absolute calibration of your DSC heat flow, only that ``y`` is proportional
to the reaction rate.

### A.1 Arbitrary T(t) — 3-column input

When a wave file carries a recorded time column (3-column `t  T  y`
input — modulated DSC, T-jump, fast-cycling FSC), the same α(T)
construction is done against the recorded ``t`` instead of against T
under a linear-heating assumption:

* ``I(t) = ∫₀^t y(t') dt'`` (trapezoid in t),
* ``α(t) = I(t) / I(t_end)``,
* ``dα/dt = y / I(t_end)`` (in 1/s),
* α(T) by interpolation along the recorded T(t).

The 2-column code path is preserved bit-for-bit when no ``t`` column is
present. Detection happens in [`io/wave_reader.py`](../src/kinetics_lems/io/wave_reader.py)
on the first non-comment line. Caveat: classical Vyazovkin (§F), KAS, and
OFW still use the *nominal* β from `settings.json` for their β-dependent
math — for strongly non-linear programs prefer Vyazovkin-AIC (§G), whose
J(E) integral is over the recorded T(t) directly.

---

## B. Kissinger (1957) — peak-temperature method

**Goal.** A single E_a from peak temperatures across heating rates.

**Algorithm.**

1. For each run, find T_p — the temperature of the maximum of dα/dt
   (equivalently the DSC/DTG peak).
2. Linear regression:

       y = ln(β / T_p²),   x = 1 / T_p,
       slope = -E_a / R,   intercept ≈ ln(A·R / E_a)   (assumes f(α) = 1 − α). (5)

**Caveat.** Equation (5) is exact only for first-order kinetics. For other
models the activation energy is recovered correctly to a few percent, but the
intercept-to-A conversion is biased. Use Kissinger as a sanity check, not as
the primary E_a estimate.

Implementation: [`methods/kissinger.py`](../src/kinetics_lems/methods/kissinger.py).

---

## C. Friedman (1964) — differential isoconversional

**Goal.** E_a(α) without assuming f(α).

For each fixed α, equation (1) gives

    ln(dα/dt)_α = ln{A·f(α)} − E_a(α) / (R·T_α).                   (6)

**Algorithm.**

1. For every run k and every α_i, look up T_{i,k} = T(α_i, run k) and
   (dα/dt)_{i,k}.
2. For each α_i, linear-regress y = ln(dα/dt)_{i,k} vs x = 1/T_{i,k}:

       slope = -E_a(α_i) / R,
       intercept ≈ ln{A · f(α_i)}.

**Pros.** Model-free, numerically simple, no temperature-integral assumption.
**Cons.** Differentiation amplifies noise.

**Optional Savitzky–Golay pre-smoothing.** Setting `smooth_window`
(odd integer ≥ `smooth_poly + 2`) in `[methods.friedman]` applies a
SavGol filter to dα/dt **per run** before α-resampling and regression.
Default is no smoothing — opt-in by config so prior outputs stay
bit-reproducible. Window 11–21 with `smooth_poly = 3` is a sensible
starting point on noisy DSC traces; larger windows distort the peak.

Implementation: [`methods/friedman.py`](../src/kinetics_lems/methods/friedman.py).

---

## D. Ozawa–Flynn–Wall (OFW) — integral isoconversional

For fixed α, integrating (3) and applying Doyle's approximation log p(x) ≈ −2.315 − 0.4567·x:

    ln β = const − 1.052 · E_a(α) / (R·T_α).                       (7)

**Algorithm.** For each α_i, regress ln β (β in K/s) vs 1/T_{α_i}:

    slope = -1.052 · E_a / R    ⇒    E_a(α) = -slope · R / 1.052.

**Caveat.** Doyle's approximation is biased outside the range x ∈ [20, 60].
Modern recommendations (ICTAC 2011) prefer KAS or Vyazovkin over OFW.

Implementation: [`methods/ofw.py`](../src/kinetics_lems/methods/ofw.py).

---

## E. Kissinger–Akahira–Sunose (KAS) — integral isoconversional

Improved Coats–Redfern integral form. For fixed α:

    ln(β / T_α²) = const − E_a(α) / (R·T_α).                       (8)

**Algorithm.** For each α_i, regress ln(β/T_α²) (β in K/s) vs 1/T_{α_i}:

    slope = -E_a / R    ⇒    E_a(α) = -slope · R.

**Why preferred over OFW.** KAS comes from an integral approximation that is
accurate over a wider x range and does not need Doyle's correction.

Implementation: [`methods/kas.py`](../src/kinetics_lems/methods/kas.py).

---

## F. Vyazovkin (1996, 1997) — nonlinear integral

**Idea.** Avoid linearization. Use the full temperature integral

    I(E, T_α)  =  (E / R) · p(x),   x = E / (R · T_α),

and exploit the fact that for any pair of runs (i, j),

    [I(E_a, T_{α,i}) / β_i]  =  [I(E_a, T_{α,j}) / β_j]            (9)

if E_a is the true activation energy at α. Define the residual

    Φ(E)  =  Σ_{i ≠ j}  [ I(E, T_{α,i}) · β_j ] / [ I(E, T_{α,j}) · β_i ]. (10)

For n runs this is n(n-1) terms (e.g. n=4 gives the 12-term sum).
At the true E_a, every ratio equals 1 and Φ(E) = n(n-1).

**Algorithm.**

1. For each α, build the n-vector { T_{α,k} } across runs.
2. Find E that minimizes Φ(E) on a 1-D bracket (default 1–600 kJ/mol),
   using bounded Brent / golden-section.
3. Repeat for every α in the grid.

**Approximation.** I uses Senum–Yang p(x), eq. (4), valid for E_a/(R·T_α)
roughly 20–100 — covers virtually all thermal-analysis kinetics.

**Validity.** Strictly assumes linear heating in each run. If runs use
arbitrary T(t) — modulated DSC, isothermal jumps — switch to the AIC variant.

Implementation: [`methods/vyazovkin.py::vyazovkin`](../src/kinetics_lems/methods/vyazovkin.py).

---

## G. Advanced Isoconversional (AIC, Vyazovkin 2001)

**Idea.** Replace the closed-form temperature integral with a numerical
integral over a small α window. Removes the linear-heating assumption and
reduces systematic error from finite Δα windows.

For each α_k pick a half-width Δα (default 0.02). Define

    J_k(E)  =  ∫_{t(α_k − Δα)}^{t(α_k + Δα)} exp(-E / (R · T(t))) dt.   (11)

Then minimize

    Φ(E)  =  Σ_{i ≠ j}  J_i(E) / J_j(E),                           (12)

exactly as in (10) but with a numerically integrated J in place of
I(E, T_α).

**Algorithm.**

1. For each run, reconstruct (t, T) from the recorded T(t) (under linear
   heating: t = (T − T_start)/β; otherwise pass the measured t directly).
2. For each α_k ∈ [Δα, 1 − Δα]:
   a. For each run, find t at α_k − Δα and α_k + Δα by interpolation.
   b. Sample T(t) on a dense grid of that t-window (default 64 points).
   c. Compute J_k(E) by trapezoid for any candidate E.
3. Minimize (12) over E on a bounded bracket (default 1–600 kJ/mol).

**Cost.** O(n_α · n_runs · n_window · n_iter); cheap (≪ 1 s for typical cases).

**When to prefer AIC over plain Vyazovkin.** Always, in principle.
The classical Vyazovkin (Senum–Yang p(x)) accumulates a small bias when E(α)
varies significantly across the α window — AIC has none.

Implementation: [`methods/vyazovkin.py::vyazovkin_aic`](../src/kinetics_lems/methods/vyazovkin.py).

---

## H. Criado–Málek master plot — Z(α) model discrimination

**Goal.** Identify f(α) by comparing experimental Z(α) to 12 theoretical
master curves.

**Definition.** Z(α) = (dα/dt) · T² is proportional to f(α)·g(α) up to
constants that cancel after normalization at α = 0.5:

    Z_norm(α) = Z(α) / Z(0.5) = f(α)·g(α) / [f(0.5)·g(0.5)]         (13)

Every master curve and the experimental curve equal 1 at α = 0.5, making
them directly comparable without knowing A, E_a, or β.

**12-model standard set** (Málek 1992 / ICTAC 2011 §5):

| Model | f(α)                                       | g(α)                          |
|-------|--------------------------------------------|-------------------------------|
| F1    | (1 − α)                                    | −ln(1 − α)                    |
| F2    | (1 − α)²                                   | 1/(1 − α) − 1                 |
| F3    | (1 − α)³                                   | ½[1/(1 − α)² − 1]             |
| A2    | 2(1 − α)[−ln(1 − α)]^(1/2)                 | [−ln(1 − α)]^(1/2)            |
| A3    | 3(1 − α)[−ln(1 − α)]^(2/3)                 | [−ln(1 − α)]^(1/3)            |
| A4    | 4(1 − α)[−ln(1 − α)]^(3/4)                 | [−ln(1 − α)]^(1/4)            |
| R2    | 2(1 − α)^(1/2)                             | 1 − (1 − α)^(1/2)             |
| R3    | 3(1 − α)^(2/3)                             | 1 − (1 − α)^(1/3)             |
| D1    | 1/(2α)                                     | α²                            |
| D2    | [−ln(1 − α)]⁻¹                             | (1−α)ln(1−α) + α              |
| D3    | 3(1−α)^(2/3) / (2[1−(1−α)^(1/3)])          | [1−(1−α)^(1/3)]²              |
| D4    | 3 / (2[(1−α)^(−1/3) − 1])                  | 1 − ⅔α − (1−α)^(2/3)          |

**Algorithm.** Average Z_k_norm across runs; compute the 12 master
curves normalized at α = 0.5; rank by RMS distance over α
(smaller = better fit).

**Known degeneracy.** F1 and A_m share identical Z_norm shapes
(Z = m·(1−α)·[−ln(1−α)] for every A_m, equal to Z_F1 up to a constant
prefactor that cancels after normalization). The Z-plot separates model
*families* (F_n vs R_n vs D_n vs A_m) but cannot distinguish F1 from
A_m by itself — combine with Friedman intercepts for that.

Implementation: [`methods/master_plot.py`](../src/kinetics_lems/methods/master_plot.py).

---

## I. Pre-exponential A from kinetic triplet

**Goal.** Given E_a(α) from Vyazovkin and a reaction model f(α), recover
A(α) = A · f(α) / f(α) = A pointwise using the rate equation:

    A(α, run) = (dα/dt)_{α,run} / [ f(α) · exp(−E_a(α) / (R · T_{α,run})) ] (15)

**Algorithm.**

1. For each α_k and each run, look up (dα/dt)_{α,k,run} and T_{α,k,run}.
2. Compute ln A_raw = ln(dα/dt) − ln f(α) + E_a / (RT) pointwise.
3. Take the median and MAD across all (α, run) pairs in log-space (robust
   to outliers near α → 0 or 1 where dα/dt is noisy).
4. Report: `A_median_per_sec` and `A_mad_per_sec` per α, plus the
   grand median across all α.

**f(α) selection.** By default, `model = ""` in the config makes the runner
auto-pick the top-ranked model from the Z(α) master plot; falls back to F1
if master_plot is disabled.

**Caveat.** Equation (15) is exact only under the assumption that E_a(α)
is a single-step quantity. If the reaction is multi-step (non-flat E(α)),
A(α) absorbs the composition of A values and loses physical meaning.
Use the multi-step result first.

Implementation: [`methods/preexponential.py`](../src/kinetics_lems/methods/preexponential.py).

---

## J. Multi-step detection via E(α) segmentation

**Goal.** Decide whether E(α) is flat enough to be single-step, and if
not, identify the α boundaries of each quasi-constant segment.

**Flatness score.**

    flatness = (max E − min E) / median E                              (16)

Thresholds (rule-of-thumb, ICTAC 2011 §6):
- < 0.05 → effectively single-step;
- 0.05–0.20 → mild multi-step;
- > 0.20 → clear multi-step.

**Algorithm (greedy segmentation).**

1. Discard NaN points. Compute `median_E` over all valid α.
2. Walk through α in order, extending the current segment while
   `|E[i] − median(segment)| ≤ threshold · |median(segment)|`.
   Start a new segment when the condition breaks.
3. Merge any segment with fewer than `min_segment_size` points into
   its neighbour (prevents single-outlier artefacts).
4. For each segment, report: `alpha_lo`, `alpha_hi`, `Ea_median`,
   `Ea_mad` (robust spread within the segment), `contribution`
   (= alpha_hi − alpha_lo).

**Recommended input.** Vyazovkin or Vyazovkin-AIC output (cleanest E(α)
profile); Friedman tends to be noisier near the tails.

Implementation: [`methods/multistep.py`](../src/kinetics_lems/methods/multistep.py).

---

## K. Reaction-order n via linearization sweep

**Goal.** For F_n kinetics — f(α) = (1 − α)ⁿ — directly recover n
without committing to a single model in advance.

**Linearization.** Starting from the Friedman form (6) and moving the
β-dependent term to the LHS:

    ln(dα/dT) + ln(β) − n · ln(1 − α)  =  ln A − E_a / (R·T)            (17)

The (constant) ln A intercept no longer depends on β, so multiple
heating rates can be **pooled** in the same regression. For the correct
n, the LHS is linear in 1/T with slope −E_a/R.

**Algorithm.**

1. For each candidate n on a grid, build the LHS over all (α, run)
   points in the working α window [α_min, α_max].
2. Linear regression LHS vs 1/T; record R² and slope.
3. Pick the n with maximum R²; report E_a = −slope · R, R²_best.

**Caveat.** Only valid for the F_n family. Avrami, contracting-geometry,
diffusion kinetics give biased n with R² < 1 — combine with the
master-plot Z(α) ranking to confirm the model family first.

Implementation: [`methods/reaction_order.py`](../src/kinetics_lems/methods/reaction_order.py).

---

## L. Coats–Redfern — single-rate model-fitting

**Goal.** Cross-check the master-plot ranking with a different
methodology. Coats–Redfern fits each model independently per heating
rate; the master plot uses ratios across rates.

**Linearization.** For each model g(α) and each rate β:

    ln[g(α) / T²]  ≈  ln(A · R / (β · E))  −  E / (R · T)               (18)

Linear in 1/T with slope −E/R and intercept ≈ ln(AR/(βE)).

**Algorithm.**

1. For every (model, run) pair, fit (18) over [α_min, α_max].
2. Extract E from the slope; back out A from the intercept assuming
   linear heating: A = β · E · exp(intercept) / R.
3. Per model, average E and log10 A across runs; rank models by mean R².

**Interpretation.** A good f(α) gives high R² *and* consistent E, A
across heating rates (small std). Disagreement signals multi-step
kinetics or a non-canonical model.

Implementation: [`methods/coats_redfern.py`](../src/kinetics_lems/methods/coats_redfern.py).

---

## M. Confidence intervals via leave-one-run-out jackknife

**Goal.** Quantify E(α) uncertainty — no closed-form covariance exists
for isoconversional methods that pool across runs.

**Formula** (Efron & Tibshirani 1993 §11, ICTAC 2011 §3). For n runs:

    SE_jack(α) = sqrt[ (n − 1) / n · Σ_i (E_(i)(α) − E̅(α))² ]          (19)

where E_(i)(α) is the estimator refit on the (n − 1)-run subset with
run i removed, and E̅ is their mean. Report mean ± 1.96·SE for ~95 % CI.

**Caveats.** Needs n ≥ 3 runs (with 2 every LOO subset has 1 run which
isoconversional methods reject — SE reported as NaN). The jackknife SE
underestimates true uncertainty for highly nonlinear estimators, but is
order-of-magnitude correct for E_a methods.

Implementation: [`methods/uncertainty.py`](../src/kinetics_lems/methods/uncertainty.py).

---

## N. Predictive isothermal kinetics — α(t) at fixed T

**Goal.** Given the kinetic triplet (E_a(α), A, f(α)), predict α(t) at a
user-chosen storage / operating temperature ("shelf life at 25 °C").

**Method.** Equation (1) at constant T rearranges to

    t(α)  =  ∫_{α_start}^{α}  dα' / [A · f(α') · exp(−E_a(α')/(R·T))]   (20)

A single cumulative trapezoid over the α grid gives the full α → t map;
linear interpolation yields time-to-α_target. Repeat at several T to
tabulate the lifetime curve.

**Caveats.** Prediction quality is bounded by (a) uncertainty in E_a(α)
— quantify with §M; (b) correctness of f(α) — verify §H and §L agree;
(c) extrapolation distance — predicting at 25 °C from 450 °C data
extrapolates the Arrhenius factor over ~10 orders of magnitude.

Reference: Vyazovkin (2000); ICTAC 2011 §7. Implementation:
[`methods/lifetime.py::predict_alpha_of_t`](../src/kinetics_lems/methods/lifetime.py).

---

## O. Model-free isothermal prediction (Vyazovkin 2000)

**Goal.** Predict t(α) at a fixed T_iso from E(α) **alone**, without
committing to an f(α). Useful when the master-plot / Coats-Redfern
ranking is ambiguous and the wrong f(α) would dominate the prediction
error of §N.

**Identity.** Under the isoconversional assumption, for any α the
"time spent at α" is

    t(α; T_iso)  =  (I(E_α, T_α) / β_exp) · exp(E_α / (R · T_iso))     (21)

where `I(E, T) ≈ (E/R)·p(x)` is the same Senum–Yang temperature integral
as in §F, evaluated at the temperature `T_α` at which conversion α was
reached in one of the *experimental* runs at rate `β_exp`. The unknown
`A · f(α)` cancels because the same α is matched on both sides.

**Algorithm.**

1. Pick one experimental run as the "reference" — typically the slowest
   heating rate, since it spans the widest T-range.
2. For every α on the grid, look up `T_α` on that run.
3. Compute `t(α; T_iso)` via (21); enforce monotonicity (numerical
   artefacts at α near the endpoints can produce tiny dips).

**When to prefer over §N.** Use §O when the f(α) family is *ambiguous*
(several master-plot models within RMS≈ε of each other) or when the
single-step assumption is suspect (non-flat E(α), §J). §N is sharper
when f(α) is well-identified.

**Arbitrary T(t) variant.** For a non-isothermal storage profile, the
same identity gives (21) extended:

    τ(α; program)  =  ∫₀^t exp(−E_α / R·T(t')) dt'                     (22)

Solve `τ(α; program) = I_ref(α) / β_exp` for ``t`` per α to obtain the
predicted t(α) under the program.

Reference: Vyazovkin (2000) eq. (16). Implementation:
[`methods/prediction_modelfree.py`](../src/kinetics_lems/methods/prediction_modelfree.py).

---

## P. α(t) under arbitrary T(t) — triplet-based extension

**Goal.** When the full kinetic triplet is known (E_a(α), A, f(α)), and
a recorded or designed T(t) profile is supplied (storage with diurnal
variation, multi-step cure cycle, modulated DSC), forward-integrate
dα/dt directly.

**Method.** Explicit midpoint rule on the supplied (time_s, T_K) grid:

    α_{k+1} = α_k + A · f(α_k) · exp(−E_a(α_k) / (R · T_mid)) · Δt     (23)

with `T_mid = ½(T_k + T_{k+1})`. E_a(α) is interpolated linearly off
the isoconversional grid; α is clipped to [0, 1].

**Caveats.** First-order accurate in Δt — fine for storage profiles
(hours/days) where T changes slowly. For modulated DSC at ~60 s
periods, the time step should resolve the modulation (≤ 5 s).
Identical f(α) and A caveats as §N apply.

Implementation: [`methods/lifetime.py::predict_under_program`](../src/kinetics_lems/methods/lifetime.py).

---

## Q. ICTAC cross-method consistency check

**Goal.** Quantify the agreement (or lack thereof) between
isoconversional methods. ICTAC 2020 §3.2 makes pairwise comparison
mandatory in a serious report — large disagreement is a red flag for
bad data (insufficient β range, poor baseline) or kinetic pathology
(multi-step, distributed reactivity, diffusion control).

**Metric.** For every unordered pair (M_i, M_j) of methods:

    Δ_ij(α)  =  | E_i(α) − E_j(α) | / max(|E_i(α)|, |E_j(α)|, ε_E)    (24)
    Δ_ij     =  max_α Δ_ij(α)                                          (25)

where `ε_E = 1 kJ/mol` (in J/mol) is a floor so that ratios stay bounded
at α-endpoints where E may be near zero from noise.

**Reporting.** Emit a warning for every pair with `Δ_ij > threshold`
(default 0.10 = 10 %, ICTAC §3.2 rough cutoff). Pairs are reported with
the worst-case α — that tells the user *where* the methods disagree,
which usually points to a tail / baseline issue rather than a method bug.

Implementation: [`methods/consistency.py`](../src/kinetics_lems/methods/consistency.py).

---

## R. Endpoint reliability + single-step diagnostic

**Goal.** Per-method reliability scoring derived from a single E(α) curve.
Three flags get raised:

1. **Low-α tail invalid.** Fraction of finite E values for α < α_low
   below threshold (default 0.5). Signal-to-noise near α = 0 is poor;
   integral methods underflow `g(α)`.
2. **High-α tail invalid.** Same for α > α_high (default 0.9). Baseline
   drift dominates near α = 1.
3. **E_a varies in the core.** Relative spread

       flatness  =  (max E − min E) / max(|median E|, 1 kJ/mol)        (26)

   over α ∈ [α_low, α_high]. Above `flatness_threshold` (default 0.10),
   the single-step assumption is flagged as suspect.

The α window comes from `[conversion]` in the TOML config. Warnings
are folded into the CLI summary and the Markdown report.

Implementation: [`methods/diagnostics.py`](../src/kinetics_lems/methods/diagnostics.py).

---

## S. AIC / BIC for model comparison

**Goal.** Add an information-theoretic ranking on top of mean R² for
both Coats–Redfern model fits (§L) and multi-step segmentations (§J).

**Formulae.** Under Gaussian residuals the maximised log-likelihood
reduces to `−n/2 · ln(RSS/n)`, so up to a model-independent constant
(Burnham & Anderson 2002 §2.2):

    AIC  =  n · ln(RSS / n)  +  2 · k                                  (27)
    BIC  =  n · ln(RSS / n)  +  k · ln(n)                              (28)

with `k` the number of free parameters. For Coats–Redfern, `k = 2`
(slope + intercept of the linear regression). For multi-step
piecewise-constant approximation of E(α), `k = n_steps`.

**Use.** Only AIC *differences* between competitors on the same data
are meaningful (the dropped constant is identical). Lower = better.
Combine with R² and physical plausibility — a tiny AIC win for an
exotic f(α) does not justify it.

Implementation: `_aic_bic_ols` in [`methods/coats_redfern.py`](../src/kinetics_lems/methods/coats_redfern.py),
`_aic_bic_piecewise` in [`methods/multistep.py`](../src/kinetics_lems/methods/multistep.py).

---

## T. Compensation effect — ln A vs E identifiability flag

**Goal.** Flag when the Arrhenius parameters (E, ln A) are statistically
*non-identifiable* on the available data.

**Background.** Across α (isoconversional A) and across reaction-model
assumptions (Coats–Redfern with all 12 models × N heating rates), the
pairs (E_i, ln A_i) often fall on a straight line:

    ln A  ≈  a  +  b · E                                               (29)

This is *not* a physical law — it's an artefact of the Arrhenius
linearization at a common temperature window. A tight line (R² → 1) is
therefore a **red flag**: it means many (E, A) combinations fit the data
equally well, and the recovered "kinetic triplet" is statistically
degenerate.

**Algorithm.** Collect every Coats–Redfern (E, A) per (model, run)
pair with A > 0; regress ln A against E.

**Reading the output.**
- R² ≳ 0.95 → strong compensation, parameters non-identifiable.
- R² ≈ 0.7–0.9 → moderate compensation, expected on smooth single-step
  data, but flag for review.
- R² < 0.7 → little compensation, parameters more independently
  identified.

Reference: Vyazovkin & Wight (2000), Galwey & Brown (2002).
Implementation: [`methods/compensation.py`](../src/kinetics_lems/methods/compensation.py).

---

## U. Empirical Z(α) fits — Sestak-Berggren and Prout-Tompkins

**Goal.** When none of the 12 canonical master-plot models matches the
experimental Z(α) well — typical for autocatalytic cure (epoxy, thermoset),
sigmoidal nucleation, or genuinely mixed-mechanism reactions — fit a
parametric empirical f(α) to Z(α) directly.

**Models.**

- **Prout-Tompkins:** `f(α) = α^m · (1 − α)^n`              (2 params)
- **Sestak-Berggren:** `f(α) = α^m · (1 − α)^n · [−ln(1−α)]^p`   (3 params)

Both interpolate between F_n, A_m and R_n at specific (m, n[, p])
values — Prout-Tompkins is a strict subset of Sestak-Berggren (p = 0).

**Algorithm.**

1. Compute experimental Z(α) the same way as §H (Criado-Málek
   normalisation at α = 0.5).
2. For a given (m, n[, p]), generate the analytical model curve:
   `f_model(α) → g_model(α)` by numerical cumulative integration of
   `1/f`, then `Z_model = f · g`, normalised at α = 0.5.
3. Least-squares fit (m, n[, p]) on bounded `scipy.optimize.least_squares`.
4. Report parameters, RMS distance from experimental Z (directly
   comparable to the master-plot RMS ranking), and R².

**Identifiability warning.** Sestak-Berggren is famously
**non-identifiable** — many (m, n, p) combinations produce nearly
identical Z(α). Always report R² and the compensation diagnostic (§T)
alongside. Prefer Prout-Tompkins when it gives a comparable R² with
fewer parameters.

**Opt-in by design.** Both fits are off by default; enable via
`[methods.empirical_models]` in `configs/default.toml`. They never
replace the 12-model master-plot ranking; they augment it.

Reference: Šesták & Berggren (1971); Prout & Tompkins (1944).
Implementation: [`methods/empirical_models.py`](../src/kinetics_lems/methods/empirical_models.py).

---

## Practical notes / pitfalls

- **Units of β.** All formulae are derived for β in K/s. The CLI / settings
  accept K/min for ergonomics; division by 60 happens in
  [`conversion.py`](../src/kinetics_lems/conversion.py).
- **Edge α values.** All isoconversional methods become noisy near α → 0 and
  α → 1 because T(α) is poorly determined there. Default α grid is [0.05, 0.95].
- **Synthetic ground truth.** All methods recover E_a to < 0.1 kJ/mol on
  noiseless single-step data and to a few kJ/mol with 3 % noise; see
  [`tests/test_synthetic_recovery.py`](../tests/test_synthetic_recovery.py).
- **Multi-step reactions.** A non-flat E(α) is the diagnostic: > 20 % variation
  signals a multi-step process (§J gives a quantitative flatness score).
  Single-E_a methods (Kissinger, Coats–Redfern) become misleading; report
  E(α) instead. Use A only after confirming single-step character.

---

## References

- Friedman, H. L. (1964). *J. Polymer Sci. Part C* 6, 183.
- Ozawa, T. (1965). *Bull. Chem. Soc. Jpn.* 38, 1881.
- Flynn, J. H.; Wall, L. A. (1966). *J. Res. NBS A* 70 (6), 487.
- Akahira, T.; Sunose, T. (1971). *Sci. Tech. Energ. Mater.* 22, 254.
- Kissinger, H. E. (1957). *Anal. Chem.* 29 (11), 1702.
- Vyazovkin, S. (1996). *Int. J. Chem. Kinet.* 28, 95.
- Vyazovkin, S. (1997). *J. Comput. Chem.* 18, 393.
- Vyazovkin, S. (2001). *J. Comput. Chem.* 22, 178 (advanced isoconversional).
- ICTAC Kinetics Committee (2011). *Thermochim. Acta* 520, 1.
- Senum, G. I.; Yang, R. T. (1977). *J. Therm. Anal.* 11, 445.
- Criado, J. M. (1978). *Thermochim. Acta* 24, 186.
- Criado, J. M. et al. (1989). *Thermochim. Acta* 147, 75.
- Málek, J. (1992). *Thermochim. Acta* 200, 257.
- Coats, A. W.; Redfern, J. P. (1964). *Nature* 201, 68.
- Vyazovkin, S. (2000). *Thermochim. Acta* 355, 155 — model-free prediction.
- Efron, B.; Tibshirani, R. J. (1993). *An Introduction to the Bootstrap.* Chapman & Hall.
- ICTAC Kinetics Committee (2020). *Thermochim. Acta* 689, 178597 — updated recommendations.
- Burnham, K. P.; Anderson, D. R. (2002). *Model Selection and Multimodel Inference.* Springer.
- Vyazovkin, S.; Wight, C. A. (2000). *Thermochim. Acta* 340-341, 53 — compensation effect.
- Galwey, A. K.; Brown, M. E. (2002). *Thermochim. Acta* 386, 91 — compensation effect critique.
- Šesták, J.; Berggren, G. (1971). *Thermochim. Acta* 3, 1 — empirical Sestak–Berggren model.
- Prout, E. G.; Tompkins, F. C. (1944). *Trans. Faraday Soc.* 40, 488 — autocatalytic model.
- Savitzky, A.; Golay, M. J. E. (1964). *Anal. Chem.* 36 (8), 1627 — polynomial smoothing filter.

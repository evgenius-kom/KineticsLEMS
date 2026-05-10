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

**Goal.** Identify f(α) by comparing the experimental Z(α) curve to the
theoretical master curves of 12 standard models.

**Definition.**

    Z(α) = (dα/dt) · T²                                               (13)

This is proportional to f(α)·g(α) up to constants that cancel after
normalization. After dividing by Z at α = 0.5:

    Z_norm(α) = Z(α) / Z(0.5) = f(α)·g(α) / [f(0.5)·g(0.5)]         (14)

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
| D4    | 3 / (2[(1−α)^(−1/3) − 1])                  | 1 − ⅔α − (1−α)^(2/3)         |

**Algorithm.**

1. For each run k, evaluate Z_k(α) = (dα/dt)_α · T_α² at every α in the grid.
2. Normalize: Z_k_norm(α) = Z_k(α) / Z_k(0.5).
3. Average across runs: Z_exp(α) = mean_k Z_k_norm(α).
4. Compute the 12 master curves and normalize each at α = 0.5.
5. Rank by RMS distance: RMSD_m = sqrt(mean_α (Z_exp − Z_master_m)²).
   Smaller RMSD → better fit.

**Known degeneracy.** F1 and A_m (m = 2, 3, 4) have identical Z_norm
shapes: Z = m·(1−α)·[−ln(1−α)] for A_m, which equals Z_F1 up to a
constant prefactor that cancels after normalization. The Z-plot reliably
separates model *families* (F_n vs R_n vs D_n vs A_m) but cannot
distinguish F1 from A_m by itself. Combine with Friedman intercepts for
that.

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

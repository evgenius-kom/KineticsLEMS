# TODO — features captured from Excel workbooks (not yet implemented)

These three feature ideas were extracted from `theory/Сводная таблица.xlsx`
and `theory/Zalfa plots Pr-Tm program.xlsm`. The Excel files themselves
can be deleted: the formulas, algorithm sketches, and references below are
sufficient to implement them later. None of these is implemented today —
this file is the spec to implement them when needed.

Cross-reference: [docs/ROADMAP_PROMPT.md](ROADMAP_PROMPT.md) section 3
covers model-based features more broadly. The three items here are the
specific Excel-derived ones that have a known formula and a clear
input/output already.

---

## 1. Criado–Málek master plot — Z(α)

**Goal.** Identify the reaction-model function f(α) by comparing the
experimental Z(α) curve against the master Z(α) curves of standard
reaction models (F1, F2, A2, R2, R3, D1–D4, …).

**Formula** (from `Сводная таблица.xlsx` sheet "Z-plots", verified
against Criado et al. 1989 / Málek 1992):

    Z(α) = (dα/dt) · T²

Two normalizations are common:

* Unnormalized: plot Z(α) vs α directly.
* Normalized: divide by Z at α = 0.5 → curves of all models intersect at α = 0.5.

`Сводная таблица.xlsx` row 3 formula was `=C3*B3^2/$H$3` — column C is
dα/dt, B is T, $H$3 is the normalizing constant Z(α=0.5). Also implemented
that way in `Zalfa plots Pr-Tm program.xlsm` as `=(C/$H$2)*(B/$G$2)^2`,
which is the same thing rewritten with separate dα/dt and T² normalizers.

**Master curves to compare against** (from Málek 1992):

| Model | f(α)                        | g(α)                          | Z(α) = f·g            |
|-------|-----------------------------|-------------------------------|-----------------------|
| F1    | 1 − α                       | −ln(1 − α)                    | (1 − α)·[−ln(1 − α)]  |
| F2    | (1 − α)²                    | (1 − α)⁻¹ − 1                 | (1 − α)·[1 − (1 − α)] |
| A2    | 2(1 − α)·[−ln(1 − α)]^(1/2) | [−ln(1 − α)]^(1/2)            | …                     |
| A3    | 3(1 − α)·[−ln(1 − α)]^(2/3) | [−ln(1 − α)]^(1/3)            | …                     |
| R2    | 2·(1 − α)^(1/2)             | 1 − (1 − α)^(1/2)             | …                     |
| R3    | 3·(1 − α)^(2/3)             | 1 − (1 − α)^(1/3)             | …                     |
| D1    | 1/(2α)                      | α²                            | α / 2                 |
| D2    | [−ln(1 − α)]⁻¹              | (1 − α)·ln(1 − α) + α         | …                     |
| D3    | 3(1 − α)^(2/3) / (2·(1 − (1 − α)^(1/3))) | …               | …                     |

**Implementation sketch (≈ 60 LOC, ≈ 1 hour):**

1. Add `methods/master_plot.py` with `def z_alpha(run, alphas, normalize=True) -> np.ndarray`.
2. Add `MASTER_CURVES: dict[str, callable]` returning Z(α) for each model.
3. Add a `plot_master_zalpha` reporter that overlays the experimental Z(α)
   from one run on top of the master curves.
4. Pick the model whose master curve sits closest to experimental — print
   ranking by mean abs distance.

**References.** Criado, J. M. et al. *Thermochim. Acta* 147 (1989) 75. Málek, J.
*Thermochim. Acta* 200 (1992) 257.

---

## 2. Pre-exponential A from per-α E_a (assuming f(α) = 1 − α)

**Goal.** Once Vyazovkin gives E_a(α), recover A under a chosen reaction
model, then average across α and rates.

**Formula** (from `Сводная таблица.xlsx` sheet "A calculated", row 4
formula `=C4 / ((1 − A4) · EXP(−D4·1000 / (8.314·B4)))`):

    A(α) = (dα/dt)_α  /  [ f(α) · exp(−E_a(α) / (R·T_α)) ]

with f(α) = 1 − α (first-order assumption baked into the spreadsheet).
`D4*1000` converts kJ/mol → J/mol because the Vyazovkin column is in kJ.

**Implementation sketch (≈ 30 LOC, ≈ 30 min):**

1. After running Vyazovkin to get `Ea_J_per_mol[α]`, compute
   `A_per_alpha[α, run] = (dα/dt)_α,run / (f(α) · exp(−E/(R·T_α,run)))`.
2. Default f(α) is F1 (matches the spreadsheet); allow other models from
   the master-plot module.
3. Report median A across α and rates with an MAD-based dispersion.

**Caveat.** Bias depends entirely on the chosen f(α). Combine with the
master-plot output above before quoting a number.

---

## 3. Reaction-order n via linearization

**Goal.** Find n in f(α) = (1 − α)ⁿ by trying several values and picking
the one that gives the straightest plot.

**Formula** (from `Сводная таблица.xlsx` sheet "Проверка порядка"):

    y(n; α, T) = ln(dα/dT) − n · ln(1 − α)
    plot y vs 1/T  for each candidate n ∈ {0.7, 1, 2, 3, …}

For the correct n this is linear with slope −E_a/R (exactly the Friedman
form rewritten in T-domain). The spreadsheet hard-codes n ∈ {1, 2, 3, 0.7}.

**Implementation sketch (≈ 25 LOC, ≈ 30 min):**

1. For a single rate (or all rates pooled), sweep n over a grid (e.g.
   0.5 → 3.0, step 0.1).
2. For each n, fit a line to (1/T, ln(dα/dT) − n·ln(1 − α)) and record R².
3. Return the n that maximizes R²; also report the slope (= −E_a/R) of
   the best fit.

**Caveat.** Only valid if the reaction is well-described by an n-th-order
model. For Avrami / diffusion / contracting-geometry kinetics this gives
biased n. Use master-plot ranking (item 1) first to confirm the model
family before quoting n.

---

## Other extra data (not migrated, recoverable from xlsx if anyone needs)

The following auxiliary data was visible in the workbooks but was not
turned into a fixture because the project doesn't currently need it:

* **СКМ (`прпепрег-СКМ`)** — preimpregnated polymer, 6 heating rates
  (1, 2.5, 5, 10, 20, 40 K/min), α grid 0.005 → 0.995 step 0.01. No
  reference E_a in the sheet. Useful as a 6-rate stress-test if
  multi-rate robustness ever becomes an issue.
* **Pr–Tm** (`Zalfa plots Pr-Tm program.xlsm`) — single-material Z(α)
  computation; only useful as a worked example for item 1 above.
* **`Сводная таблица.xlsx` sheet "Новая кинетика"** — newer ABPOT 2016
  data with finer α grid. Possibly higher-quality replacement for the
  current ABPOT fixture; would need re-extraction if re-used.

# Epoxy PV15-0,5

Source: `theory/extra/Вязовкин.xlsx` sheet `Эпоксидка PV15-0,5`.

Epoxy resin with low (E_a ≈ 50 kJ/mol) activation energy that
*decreases* with α — classic Avrami-style autocatalytic-looking
fingerprint that is interesting because it is *different* from the
ABPOT/DAPOT polymers which have flat E_a(α).

## Data files

| File | Columns | Notes |
|------|---------|-------|
| `rate_2.5.tsv`, `rate_5.0.tsv`, `rate_10.0.tsv` | `alpha`, `T_K` | (α, T) per heating rate |
| `reference_Ea.tsv` | `alpha`, `Ea_kJ_per_mol` | independent Excel Vyazovkin reference |

Heating rates: 2.5, 5.0, 10.0 K/min. α grid step ≈ 0.01.

## Applicable tests

* **Integral-method validation** under declining E_a(α) — KAS / OFW /
  Vyazovkin should track the reference even though the E curve is not
  flat. Catches drift on materials where multi-step character is mild.
* **Reaction-model identification** — Z(α) master-plot should rank
  Avrami / nth-order higher than diffusion models, which the reference
  E_a profile expects.

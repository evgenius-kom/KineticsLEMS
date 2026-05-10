# ДАПОТ reference fixtures

Sources:

- ``theory/Вязовкин.xlsx`` sheet ``ДАПОТ`` — provides α(T) at β = 2.5, 5, 10 K/min
  and the reference E_a column.
- ``theory/для расчета по Вязовкину.xlsm`` sheet ``4 скорости ДАПОТ`` —
  provides the additional β = 1 K/min column (``rate_1.0.tsv``).

Each ``rate_<β>.tsv`` file contains the (α, T_K) pairs from a single
heating-rate column of the source sheet, sorted by α and clipped into [0, 1].

``reference_Ea.tsv`` contains the per-α activation energy from the
"Ea, kJ" column of ``Вязовкин.xlsx`` — produced by an independent
Excel-based Vyazovkin computation **using only the 3 rates 2.5/5/10**.
Used by ``tests/test_reference_abpot_dapot.py`` to detect drift in our
integral methods (KAS / OFW / Vyazovkin).

Heating rates available: [1.0, 2.5, 5.0, 10.0] K/min.
Rates used by the regression test: [2.5, 5.0, 10.0] (so the comparison
matches the reference's exact configuration).

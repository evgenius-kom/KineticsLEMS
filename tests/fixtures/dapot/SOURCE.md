# ДАПОТ reference fixtures

Source: ``theory/Вязовкин.xlsx`` sheet ``ДАПОТ``.

Each ``rate_<β>.tsv`` file contains the (α, T_K) pairs from a single
heating-rate column of the source sheet, sorted by α and clipped into [0, 1].

``reference_Ea.tsv`` contains the per-α activation energy from the
"Ea, kJ" column of the source sheet — produced by an independent
Excel-based Vyazovkin computation. Used by
``tests/test_reference_abpot_dapot.py`` to detect drift in our integral
methods (KAS / OFW / Vyazovkin) against this independent reference.

Heating rates: [2.5, 5.0, 10.0] K/min.

"""Read a two- or three-column text file into a Wave.

* 2-column: ``T  y`` (or ``x  y`` in general) — the historical format,
  assumes linear heating; the analysis pipeline reconstructs t from
  ``T = T₀ + β·t`` using the rate from settings.json.

* 3-column: ``t  T  y`` — arbitrary T(t) program (modulated DSC, T-jump,
  fast-cycling FSC). ``t`` in seconds, ``T`` in K, ``y`` proportional to
  rate or dα/dT. The recorded ``t`` is preserved on the :class:`Wave`
  and used downstream instead of the linear-heating assumption.

Detection: first non-comment line wins — if it parses to ≥3 floats the
file is treated as 3-column for the entire scan; otherwise 2-column.
This is the same convention NETZSCH Kinetics Neo uses for its
user-defined data import.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

from ..models import Wave


def read_wave(path: Path | str) -> Wave:
    """Read a 2- or 3-column whitespace/comma-delimited text file into a Wave.

    Lines that don't parse are skipped silently.
    """
    path = Path(path)
    parsed_rows: list[list[float]] = []
    detected_cols: int | None = None

    with path.open("r", encoding="utf-8") as f:
        for raw in f:
            line = raw.strip().replace(",", ".")
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            try:
                row = [float(p) for p in parts[:3]]
            except ValueError:
                continue
            if detected_cols is None:
                detected_cols = 3 if len(row) >= 3 else 2
            parsed_rows.append(row[:detected_cols])

    if not parsed_rows:
        raise ValueError(f"No numeric data found in {path}")

    arr = np.asarray(parsed_rows, dtype=float)
    if detected_cols == 3:
        t = arr[:, 0]
        T = arr[:, 1]
        y = arr[:, 2]
        return Wave(x=T, y=y, t_seconds=t)
    return Wave(x=arr[:, 0], y=arr[:, 1])

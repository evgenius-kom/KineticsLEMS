"""Read a two-column (x, y) text file into a Wave."""
from __future__ import annotations

from pathlib import Path

import numpy as np

from ..models import Wave


def read_wave(path: Path | str) -> Wave:
    """Read whitespace- or comma-decimal text file with two columns into a Wave.

    Lines that don't parse as two floats are skipped silently.
    """
    path = Path(path)
    xs: list[float] = []
    ys: list[float] = []
    with path.open("r", encoding="utf-8") as f:
        for raw in f:
            line = raw.strip().replace(",", ".")
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            try:
                xs.append(float(parts[0]))
                ys.append(float(parts[1]))
            except ValueError:
                continue
    if not xs:
        raise ValueError(f"No numeric data found in {path}")
    return Wave(np.asarray(xs, dtype=float), np.asarray(ys, dtype=float))

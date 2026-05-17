"""AI/ML auxiliaries — placeholder, not used by the core pipeline.

This package is intentionally minimal and **must not** become a hard
dependency of the thermal-kinetics core. The whole point of the ICTAC
recommendations the core implements is *defensible, validated kinetics*;
ML enters only as:

1. **Baseline / peak-detection assistant** — suggest baseline regions and
   peak/shoulder count to the user, never apply changes silently.
2. **Initial-guess generator** — derive sensible (E, A, n) starting points
   for model-based fitting from isoconversional outputs and peak features.
3. **Surrogate model** — fast prediction-after-fit for parameter studies,
   trained only on validated fits.

What is intentionally *out of scope* (mirroring our A=solid-state decision):

* Enzyme kinetic predictors (CatPred / UniKP / DLKcat / DLTKcat / CataPro) —
  these belong to bio-kinetics, not solid-state TGA/DSC.
* Mechanism generation (RMG-Py / ARC / AutoTST).
* Property-prediction graph neural networks for gas-phase reactions.

If the project ever expands beyond solid-state thermal analysis, the
enzyme / mechanism / GNN integrations would live as separate optional
sub-packages each with their own dependency footprint, all sharing the
:class:`PredictorPlugin` contract below.

Status: **placeholder**. The classes here document the eventual plugin
contract; none of them do anything yet.
"""
from __future__ import annotations

from .plugin import (
    PredictorMetadata,
    PredictorPlugin,
    PredictorPrediction,
    PredictorRegistryError,
    register_predictor,
    registry,
)

__all__ = [
    "PredictorMetadata",
    "PredictorPlugin",
    "PredictorPrediction",
    "PredictorRegistryError",
    "register_predictor",
    "registry",
]

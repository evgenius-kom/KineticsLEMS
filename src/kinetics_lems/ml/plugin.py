"""Predictor plugin contract — see :mod:`kinetics_lems.ml` for context.

Every ML predictor that wants to feed into the pipeline (baseline
assistant, initial-guess generator, surrogate model, ...) implements
:class:`PredictorPlugin` and registers itself with :func:`register_predictor`.
The contract is intentionally narrow:

* metadata up-front (name, version, domain, license);
* a single ``predict()`` call returning typed :class:`PredictorPrediction`
  with mandatory ``source`` provenance and optional uncertainty.

The core pipeline never imports this module. The CLI / report layer can
optionally consult ``registry()`` and annotate outputs with predictions,
clearly tagged as such.
"""
from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from typing import Any


class PredictorRegistryError(RuntimeError):
    """Raised when two plugins try to register under the same name."""


@dataclass(frozen=True)
class PredictorMetadata:
    name: str
    version: str
    description: str
    domain: str  # e.g. "baseline_assist", "initial_guess", "surrogate"
    license: str | None = None
    reference: str | None = None
    implemented: bool = False


@dataclass(frozen=True)
class PredictorPrediction:
    """One predictor's output. Always tagged with provenance.

    The pipeline must never silently mix predicted values with
    experimentally fitted ones — `source` is what makes that auditable.
    """

    source: str  # plugin name + version
    input_hash: str  # SHA-256 of the input payload
    value: Any
    uncertainty: float | None = None
    domain_warnings: list[str] = field(default_factory=list)


class PredictorPlugin(ABC):
    """Base class for any ML-backed auxiliary."""

    METADATA: PredictorMetadata

    @abstractmethod
    def predict(self, payload: Any) -> PredictorPrediction:
        """Run the model on ``payload`` and return a tagged prediction.

        Implementations should:
        - hash ``payload`` deterministically into ``PredictorPrediction.input_hash``;
        - report ``domain_warnings`` whenever the input is outside the model's
          training distribution (sequence identity < N%, T outside training
          range, unusual atmosphere, etc.).
        """


# --- Registry ---------------------------------------------------------------

_REGISTRY: dict[str, type[PredictorPlugin]] = {}


def register_predictor(cls: type[PredictorPlugin]) -> type[PredictorPlugin]:
    """Class decorator: add ``cls`` to the global predictor registry."""
    name = cls.METADATA.name
    if name in _REGISTRY:
        raise PredictorRegistryError(f"Predictor already registered: {name!r}")
    _REGISTRY[name] = cls
    return cls


def registry() -> dict[str, type[PredictorPlugin]]:
    """Shallow copy of the predictor registry."""
    return dict(_REGISTRY)


__all__ = [
    "PredictorMetadata",
    "PredictorPlugin",
    "PredictorPrediction",
    "PredictorRegistryError",
    "register_predictor",
    "registry",
]

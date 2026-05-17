"""Smoke tests for the ML plugin contract."""
from __future__ import annotations

import pytest

from kinetics_lems.ml import (
    PredictorMetadata,
    PredictorPlugin,
    PredictorPrediction,
    PredictorRegistryError,
    register_predictor,
    registry,
)


def test_registry_starts_empty() -> None:
    assert registry() == {}


def test_register_and_lookup() -> None:
    @register_predictor
    class _Toy(PredictorPlugin):
        METADATA = PredictorMetadata(
            name="toy",
            version="0.1.0",
            description="constant predictor",
            domain="surrogate",
        )

        def predict(self, payload):
            return PredictorPrediction(
                source="toy@0.1.0",
                input_hash="x",
                value=42.0,
            )

    try:
        assert "toy" in registry()
        result = _Toy().predict({"any": "thing"})
        assert result.value == 42.0
        assert result.source == "toy@0.1.0"
    finally:
        # Keep the global registry clean for other tests.
        from kinetics_lems.ml.plugin import _REGISTRY

        _REGISTRY.pop("toy", None)


def test_duplicate_registration_raises() -> None:
    @register_predictor
    class _A(PredictorPlugin):
        METADATA = PredictorMetadata(
            name="dup_test", version="0.1.0", description="", domain="x"
        )

        def predict(self, payload):
            return PredictorPrediction(source="x", input_hash="", value=None)

    try:
        with pytest.raises(PredictorRegistryError):

            @register_predictor
            class _B(PredictorPlugin):
                METADATA = PredictorMetadata(
                    name="dup_test", version="0.2.0", description="", domain="x"
                )

                def predict(self, payload):
                    return PredictorPrediction(source="y", input_hash="", value=None)
    finally:
        from kinetics_lems.ml.plugin import _REGISTRY

        _REGISTRY.pop("dup_test", None)

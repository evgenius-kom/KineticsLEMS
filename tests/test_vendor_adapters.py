"""Vendor-adapter infrastructure: registry, generic dispatch, stubs error politely."""
from __future__ import annotations

from pathlib import Path

import pytest

from kinetics_lems.io.vendors import (
    GenericKineticsLEMSAdapter,
    NetzschKineticsNeoAdapter,
    NetzschProteusAdapter,
    TATriosJSONAdapter,
    VendorAdapterError,
    registry,
)

EXAMPLE_CASE = Path(__file__).resolve().parent.parent / "examples" / "synthetic" / "F1_120kJ"


def test_registry_lists_all_known_vendors() -> None:
    names = set(registry())
    # If a new vendor adapter is added, update this assertion intentionally.
    assert names == {
        "akts",
        "generic",
        "mettler_stare",
        "netzsch_kinetics_neo",
        "netzsch_proteus",
        "perkinelmer_pyris",
        "shimadzu",
        "ta_trios_json",
    }


def test_only_generic_adapter_is_implemented_today() -> None:
    implemented = {
        name for name, cls in registry().items() if cls.FORMAT.implemented
    }
    assert implemented == {"generic"}


def test_generic_adapter_loads_synthetic_case() -> None:
    if not EXAMPLE_CASE.is_dir():
        pytest.skip(f"Example case missing at {EXAMPLE_CASE}")
    case = GenericKineticsLEMSAdapter.load(EXAMPLE_CASE)
    assert case.params.material == "synthetic-F1"
    assert len(case.waves) >= 2


def test_generic_adapter_can_handle_directory_with_settings() -> None:
    if not EXAMPLE_CASE.is_dir():
        pytest.skip(f"Example case missing at {EXAMPLE_CASE}")
    assert GenericKineticsLEMSAdapter.can_handle(EXAMPLE_CASE)


def test_generic_adapter_rejects_unrelated_directory(tmp_path: Path) -> None:
    assert not GenericKineticsLEMSAdapter.can_handle(tmp_path)


@pytest.mark.parametrize(
    "adapter_cls",
    [NetzschProteusAdapter, NetzschKineticsNeoAdapter, TATriosJSONAdapter],
)
def test_stub_adapters_raise_with_helpful_message(adapter_cls, tmp_path: Path) -> None:
    fake = tmp_path / "fake_export.txt"
    fake.write_text("noop", encoding="utf-8")
    with pytest.raises(VendorAdapterError) as exc:
        adapter_cls.load(fake)
    # Error message must point users to the working fallback.
    assert "GenericKineticsLEMSAdapter" in str(exc.value)


def test_vendor_format_metadata_is_complete() -> None:
    for name, cls in registry().items():
        fmt = cls.FORMAT
        assert fmt.vendor == name
        assert fmt.description.strip(), f"empty description for {name!r}"
        assert fmt.extensions, f"empty extension list for {name!r}"

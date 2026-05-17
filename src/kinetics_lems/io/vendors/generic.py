"""Generic KineticsLEMS adapter — the format actually used today.

Delegates to the existing :func:`kinetics_lems.io.case_loader.load_case`.
This wrapper exists so the public dispatch path (vendor registry) is the
same as for any future vendor-specific importer.
"""
from __future__ import annotations

from pathlib import Path

from ...models import CaseData
from ..case_loader import load_case
from .base import VendorAdapter, VendorFormat, register


@register
class GenericKineticsLEMSAdapter(VendorAdapter):
    FORMAT = VendorFormat(
        vendor="generic",
        extensions=(".zip",),  # also accepts directories; see can_handle below
        description=(
            "Native KineticsLEMS layout: directory or .zip containing "
            "settings.json + per-rate two-column (T, y) .txt files. "
            "T in K, y proportional to dα/dT (heat-flow or mass-loss). "
            "See README.md for the schema."
        ),
        reference=None,
        implemented=True,
    )

    @classmethod
    def load(cls, path: Path | str) -> CaseData:
        return load_case(path)

    @classmethod
    def can_handle(cls, path: Path | str) -> bool:
        p = Path(path)
        if p.is_dir() and (p / "settings.json").is_file():
            return True
        return p.is_file() and p.suffix.lower() == ".zip"


__all__ = ["GenericKineticsLEMSAdapter"]

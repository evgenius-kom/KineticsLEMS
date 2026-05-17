"""PerkinElmer Pyris — ASCII/CSV export.

NOT YET IMPLEMENTED. Stub only.
"""
from __future__ import annotations

from pathlib import Path

from ...models import CaseData
from .base import VendorAdapter, VendorAdapterError, VendorFormat, register


@register
class PerkinElmerPyrisAdapter(VendorAdapter):
    FORMAT = VendorFormat(
        vendor="perkinelmer_pyris",
        extensions=(".txt", ".csv"),
        description=(
            "PerkinElmer Pyris DSC/TGA ASCII export. Header block contains "
            "sample/method metadata; data section is typically Temperature "
            "(°C or K), Time, Heat Flow or Mass."
        ),
        reference=None,
        implemented=False,
    )

    @classmethod
    def load(cls, path: Path | str) -> CaseData:
        raise VendorAdapterError(
            "PerkinElmer Pyris adapter is not yet implemented."
        )


__all__ = ["PerkinElmerPyrisAdapter"]

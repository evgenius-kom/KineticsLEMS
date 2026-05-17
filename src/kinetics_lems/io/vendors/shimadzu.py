"""Shimadzu TA software — ASCII/CSV export.

NOT YET IMPLEMENTED. Stub only.
"""
from __future__ import annotations

from pathlib import Path

from ...models import CaseData
from .base import VendorAdapter, VendorAdapterError, VendorFormat, register


@register
class ShimadzuAdapter(VendorAdapter):
    FORMAT = VendorFormat(
        vendor="shimadzu",
        extensions=(".txt", ".csv"),
        description=(
            "Shimadzu DSC/TGA ASCII export. Multiple comma-separated columns "
            "after a small header; exact column set depends on instrument model."
        ),
        reference=None,
        implemented=False,
    )

    @classmethod
    def load(cls, path: Path | str) -> CaseData:
        raise VendorAdapterError(
            "Shimadzu adapter is not yet implemented."
        )


__all__ = ["ShimadzuAdapter"]

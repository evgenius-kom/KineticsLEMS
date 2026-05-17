"""Mettler Toledo STARe — ASCII export.

NOT YET IMPLEMENTED. Stub only.
"""
from __future__ import annotations

from pathlib import Path

from ...models import CaseData
from .base import VendorAdapter, VendorAdapterError, VendorFormat, register


@register
class MettlerSTAReAdapter(VendorAdapter):
    FORMAT = VendorFormat(
        vendor="mettler_stare",
        extensions=(".txt", ".csv"),
        description=(
            "Mettler Toledo STARe TGA/DSC ASCII export. Format depends entirely "
            "on the lab-defined export protocol; typical layout is a 5-10 line "
            "header (Method, Sample, Date, ...) followed by tab-delimited "
            "Time/s, Temp/°C, Sample/mg or HeatFlow/mW columns. Multiple curves "
            "per file are possible when STARe is set to combine runs."
        ),
        reference=(
            "https://www.mt.com/gb/en/home/products/Laboratory_Analytics_Browse/TA_Family_Browse/TA_software_browse.html"
        ),
        implemented=False,
    )

    @classmethod
    def load(cls, path: Path | str) -> CaseData:
        # TODO: implement after collecting at least 2-3 export samples covering
        # different export protocols (each lab configures STARe differently).
        raise VendorAdapterError(
            "Mettler STARe adapter is not yet implemented. "
            "Pre-process to 2-column .txt and use the GenericKineticsLEMSAdapter."
        )


__all__ = ["MettlerSTAReAdapter"]

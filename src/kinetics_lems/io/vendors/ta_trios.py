"""TA Instruments TRIOS — JSON export.

NOT YET IMPLEMENTED. Stub only — see TODOs below.
"""
from __future__ import annotations

from pathlib import Path

from ...models import CaseData
from .base import VendorAdapter, VendorAdapterError, VendorFormat, register


@register
class TATriosJSONAdapter(VendorAdapter):
    FORMAT = VendorFormat(
        vendor="ta_trios_json",
        extensions=(".json",),
        description=(
            "TA Instruments TRIOS software exports DSC/TGA/SDT/rheology runs "
            "as JSON. The schema is documented by TA in their 'TB105' note: "
            "top-level Sample / Procedure / Steps blocks plus an Arrays section "
            "with the raw (Time, Temp, Heat Flow) columns. One TRIOS JSON file "
            "is typically ONE experiment; multi-rate cases need several files "
            "matched by heating rate from the Procedure block."
        ),
        reference=(
            "https://www.tainstruments.com/applications-notes/trios-software-and-json-export-overview-of-json-export-from-trios-software-and-import-into-python-tb105/"
        ),
        implemented=False,
    )

    @classmethod
    def load(cls, path: Path | str) -> CaseData:
        # TODO: implement once a real TRIOS JSON sample is in the repo.
        # Implementation outline:
        #   1. json.load → top-level dict.
        #   2. Pull metadata: Sample.Name → material, Procedure.HeatingRate → β.
        #   3. Arrays section: column names → indices; extract Temp (°C → K),
        #      Time (min → s if needed), and the primary signal (DSC or TGA).
        #   4. If multiple files are passed, group by β and produce a CaseData
        #      with one wave per heating rate.
        #   5. Preserve provenance under CaseData.params (when CaseParams gets
        #      an optional `metadata` dict — see schemas/ TODOs).
        raise VendorAdapterError(
            "TA TRIOS JSON adapter is not yet implemented. "
            "Export your TRIOS data as ASCII/CSV and use the GenericKineticsLEMSAdapter."
        )


__all__ = ["TATriosJSONAdapter"]

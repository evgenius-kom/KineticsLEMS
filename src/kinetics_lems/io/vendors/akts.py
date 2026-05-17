"""AKTS Thermokinetics — ASCII export.

NOT YET IMPLEMENTED. Stub only.
"""
from __future__ import annotations

from pathlib import Path

from ...models import CaseData
from .base import VendorAdapter, VendorAdapterError, VendorFormat, register


@register
class AKTSAdapter(VendorAdapter):
    FORMAT = VendorFormat(
        vendor="akts",
        extensions=(".txt", ".csv"),
        description=(
            "AKTS Thermokinetics (TK) accepts vendor-neutral 2-/3-column ASCII "
            "input and produces similar exports for its model-free / model-fitting "
            "results. Useful as a cross-check destination format: emit our results "
            "in an AKTS-compatible layout so users can re-import them."
        ),
        reference=(
            "https://www.akts.com/tk/thermokinetics-software-thermal-analysis-isoconversional-model-fitting-dsc-tg-detailed-description/"
        ),
        implemented=False,
    )

    @classmethod
    def load(cls, path: Path | str) -> CaseData:
        # TODO: AKTS input is already close to our generic format; the main work
        # is metadata mapping. Likely the AKTSAdapter ends up being a thin alias
        # over the generic adapter once we settle the metadata schema.
        raise VendorAdapterError(
            "AKTS adapter is not yet implemented. Use the GenericKineticsLEMSAdapter "
            "for plain ASCII; the AKTS export layout is essentially compatible."
        )


__all__ = ["AKTSAdapter"]

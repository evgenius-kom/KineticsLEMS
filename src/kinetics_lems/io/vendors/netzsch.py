"""NETZSCH instrument family — Proteus / Kinetics Neo exports.

NOT YET IMPLEMENTED. Stubs only — see TODOs below.
"""
from __future__ import annotations

from pathlib import Path

from ...models import CaseData
from .base import VendorAdapter, VendorAdapterError, VendorFormat, register


@register
class NetzschProteusAdapter(VendorAdapter):
    FORMAT = VendorFormat(
        vendor="netzsch_proteus",
        extensions=(".txt", ".csv"),
        description=(
            "Proteus ASCII export from NETZSCH STA/DSC/TGA instruments. "
            "Multi-line text header (Sample, Instrument, Atmosphere, Heating Rate, ...) "
            "followed by a tab- or semicolon-delimited table whose first columns are "
            "typically Temp/°C, Time/min, DSC/(mW/mg) or Mass/% — the exact column "
            "set depends on the export protocol the lab configured."
        ),
        reference="https://www.netzsch.com/en/products/thermal-analysis/proteus-software/",
        implemented=False,
    )

    @classmethod
    def load(cls, path: Path | str) -> CaseData:
        # TODO: implement once a real Proteus export sample is in the repo.
        # Implementation outline:
        #   1. Read the header block (lines starting with '##' or a sentinel
        #      key:value section); pull out heating rate, sample mass, gas.
        #   2. Detect delimiter (tab / semicolon / comma) on the first data line.
        #      Decimal separator may be ',' (German locale) — normalise.
        #   3. Locate the Temp column (auto-K conversion: °C + 273.15) and the
        #      signal column based on either the technique (DSC/TGA) or explicit
        #      user mapping in settings.json.
        #   4. Each Proteus file is ONE heating rate — caller supplies the rate
        #      via settings.json or it is read from the header.
        raise VendorAdapterError(
            "NETZSCH Proteus adapter is not yet implemented. "
            "Export your data to the generic 2-column (T, y) .txt format and "
            "use the GenericKineticsLEMSAdapter for now."
        )


@register
class NetzschKineticsNeoAdapter(VendorAdapter):
    FORMAT = VendorFormat(
        vendor="netzsch_kinetics_neo",
        extensions=(".txt", ".csv"),
        description=(
            "NETZSCH Kinetics Neo user-defined data import format. Either "
            "2-column (T °C, signal) or 3-column (Time, T °C, signal) ASCII; "
            "comma decimal separator common. One file per heating rate; an "
            "external 'temperature profile' is supported as a separate (t, T) CSV "
            "for arbitrary T(t) prediction."
        ),
        reference=(
            "https://kinetics.netzsch.com/en/learn/import-of-user-defined-data"
        ),
        implemented=False,
    )

    @classmethod
    def load(cls, path: Path | str) -> CaseData:
        # TODO: implement once at least one Kinetics Neo example is available.
        # Specific differences from Proteus:
        #   - Kinetics Neo files tend to be minimal-header (almost like ours).
        #   - 3-column files appear when modulated/arbitrary T(t) was used.
        #     Once the 3-column wave_reader path is in place (item #7), this
        #     adapter mostly becomes a header-stripper.
        raise VendorAdapterError(
            "NETZSCH Kinetics Neo adapter is not yet implemented. "
            "Strip the header lines manually and use the GenericKineticsLEMSAdapter."
        )


__all__ = ["NetzschKineticsNeoAdapter", "NetzschProteusAdapter"]

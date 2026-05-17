"""Vendor-specific input adapters for thermal-analysis data.

The :class:`VendorAdapter` base defines the contract every importer must
satisfy: take a path, return a :class:`CaseData` ready for the rest of
the pipeline.

Currently only the *generic* KineticsLEMS adapter (settings.json + 2-column
.txt waves) is implemented. The vendor-specific stubs raise
``NotImplementedError`` with a description of the format they will eventually
parse, so the registry already exposes the full surface area and downstream
code can dispatch on vendor name today.

The shared format conventions every adapter resolves into:

* temperature → Kelvin
* time → seconds
* heating rate → K/min in :class:`CaseParams.file_to_condition`
* signal → arbitrary (proportional to dα/dT or dα/dt; see
  :func:`kinetics_lems.conversion._build_run` for why absolute calibration
  doesn't matter for isoconversional analysis)

TODO (cross-cutting):
    - Replace dict-based metadata with the pydantic schemas in
      :mod:`kinetics_lems.schemas` once that module is promoted out of
      stub status.
    - Add ``--vendor <name>`` flag to the CLI for explicit dispatch.
    - Add ``sniff(path)`` heuristics so the CLI can auto-detect vendor.
"""
from __future__ import annotations

from .akts import AKTSAdapter
from .base import VendorAdapter, VendorAdapterError, VendorFormat, register, registry
from .generic import GenericKineticsLEMSAdapter
from .mettler_stare import MettlerSTAReAdapter
from .netzsch import NetzschKineticsNeoAdapter, NetzschProteusAdapter
from .perkinelmer import PerkinElmerPyrisAdapter
from .shimadzu import ShimadzuAdapter
from .ta_trios import TATriosJSONAdapter

__all__ = [
    "AKTSAdapter",
    "GenericKineticsLEMSAdapter",
    "MettlerSTAReAdapter",
    "NetzschKineticsNeoAdapter",
    "NetzschProteusAdapter",
    "PerkinElmerPyrisAdapter",
    "ShimadzuAdapter",
    "TATriosJSONAdapter",
    "VendorAdapter",
    "VendorAdapterError",
    "VendorFormat",
    "register",
    "registry",
]

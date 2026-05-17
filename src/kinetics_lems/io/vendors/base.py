"""Base class and registry for vendor-specific input adapters."""
from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass
from pathlib import Path

from ...models import CaseData


class VendorAdapterError(RuntimeError):
    """Raised when a vendor adapter cannot parse the supplied path."""


@dataclass(frozen=True)
class VendorFormat:
    """Static description of a vendor's data export — intentionally lightweight.

    Used both as documentation in code and as the source of the ``--help``-style
    listing that ``kinetics-lems vendors`` will eventually print.
    """

    vendor: str
    """Display name of the vendor / instrument family."""

    extensions: tuple[str, ...]
    """File extensions this adapter understands (lower-case, with leading dot)."""

    description: str
    """One-paragraph description of the export format and quirks."""

    reference: str | None = None
    """URL to vendor documentation, if any."""

    implemented: bool = False
    """``True`` if the adapter can load files today; ``False`` for stubs."""


class VendorAdapter(ABC):
    """Convert a vendor's exported file/folder into a canonical :class:`CaseData`."""

    #: subclasses MUST set ``FORMAT`` to a :class:`VendorFormat` describing them.
    FORMAT: VendorFormat

    @classmethod
    @abstractmethod
    def load(cls, path: Path | str) -> CaseData:
        """Parse ``path`` (file or folder) and return canonical :class:`CaseData`."""

    @classmethod
    def can_handle(cls, path: Path | str) -> bool:
        """Best-effort heuristic: does this adapter recognize ``path``?

        Default implementation matches by extension. Subclasses with stronger
        signals (magic bytes, sentinel header lines, manifest files) should
        override.
        """
        p = Path(path)
        if p.is_dir():
            return False
        return p.suffix.lower() in cls.FORMAT.extensions


# --- Registry ---------------------------------------------------------------

_REGISTRY: dict[str, type[VendorAdapter]] = {}


def register(adapter_cls: type[VendorAdapter]) -> type[VendorAdapter]:
    """Class decorator: add ``adapter_cls`` to the global vendor registry."""
    name = adapter_cls.FORMAT.vendor
    if name in _REGISTRY:
        raise ValueError(f"Vendor adapter already registered: {name!r}")
    _REGISTRY[name] = adapter_cls
    return adapter_cls


def registry() -> dict[str, type[VendorAdapter]]:
    """Return a shallow copy of the vendor-name → adapter-class registry."""
    return dict(_REGISTRY)


__all__ = [
    "VendorAdapter",
    "VendorAdapterError",
    "VendorFormat",
    "register",
    "registry",
]

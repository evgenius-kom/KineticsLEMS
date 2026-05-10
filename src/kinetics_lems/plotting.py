"""Paper-grade plot styling and multi-format saving.

Use this module from any reporter/plotter to get a consistent look across
figures (font sizes, palette, line widths, grid) and to save the same
figure to PNG + PDF + SVG with one call.

Example::

    from .plotting import paper_style, save_figure

    with paper_style():
        fig, ax = plt.subplots()
        ax.plot(x, y)
        save_figure(fig, output_dir / "plot", formats=("png", "pdf"))
"""
from __future__ import annotations

from collections.abc import Iterable
from contextlib import contextmanager
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from cycler import cycler

# ColorBrewer "Set1"-inspired palette: distinct, print-friendly, color-blind
# safe up to ~5 categories. Falls back to matplotlib's tab10 beyond.
_PAPER_PALETTE: tuple[str, ...] = (
    "#1f78b4",  # blue
    "#e31a1c",  # red
    "#33a02c",  # green
    "#ff7f00",  # orange
    "#6a3d9a",  # purple
    "#b15928",  # brown
    "#a6cee3",  # light blue
    "#fb9a99",  # light red
)

DEFAULT_FORMATS: tuple[str, ...] = ("png", "pdf")


def palette() -> tuple[str, ...]:
    """Return the project's default color palette."""
    return _PAPER_PALETTE


@contextmanager
def paper_style(dpi: int = 300, font_size: int = 11):
    """Temporarily apply publication-grade matplotlib rcParams.

    Defaults are tuned for single-column journal figures (~3.5 in wide):

    * sans-serif font, 11 pt body / 10 pt ticks / 9 pt legend;
    * 1.2 pt lines, 5 pt markers, 0.8 pt axis spines;
    * subtle grid (gray, dashed, alpha=0.3);
    * tight tick paddings, no top/right spines.
    """
    rc_overrides: dict[str, object] = {
        "figure.dpi": dpi,
        "savefig.dpi": dpi,
        "savefig.bbox": "tight",
        "font.family": "sans-serif",
        "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
        "font.size": font_size,
        "axes.labelsize": font_size,
        "axes.titlesize": font_size,
        "xtick.labelsize": font_size - 1,
        "ytick.labelsize": font_size - 1,
        "legend.fontsize": font_size - 2,
        "axes.linewidth": 0.8,
        "lines.linewidth": 1.2,
        "lines.markersize": 5,
        "axes.grid": True,
        "grid.linestyle": ":",
        "grid.alpha": 0.4,
        "grid.color": "#888888",
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.prop_cycle": cycler(color=list(_PAPER_PALETTE)),
        "axes.formatter.use_mathtext": True,
        "pdf.fonttype": 42,  # editable text in vector PDFs (Type 42 / TrueType)
        "ps.fonttype": 42,
        "svg.fonttype": "none",
    }
    with plt.rc_context(rc_overrides):
        yield


def save_figure(
    fig,
    base_path: Path,
    formats: Iterable[str] = DEFAULT_FORMATS,
) -> list[Path]:
    """Save ``fig`` to multiple formats sharing the same stem.

    ``base_path`` is treated as a stem (with or without extension); the
    extension is replaced for each requested format.
    Returns the list of paths actually written.
    """
    base_path = Path(base_path)
    base_path.parent.mkdir(parents=True, exist_ok=True)
    written: list[Path] = []
    for fmt in formats:
        out = base_path.with_suffix(f".{fmt.lstrip('.')}")
        fig.savefig(out)
        written.append(out)
    plt.close(fig)
    return written


__all__ = ["DEFAULT_FORMATS", "paper_style", "palette", "save_figure"]

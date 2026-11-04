# post_processing/graphical_visualisers/figure_style.py
"""One figure style for the examples and validation plots, and the marks that say what each value is.

Every value in a 1d-fem figure is one of a few kinds of finite element data, and each kind
has its own mark, so a reader can tell computed values from interpolated ones and from
references:

========================  ==========================================================
mark                      what it is
========================  ==========================================================
thin solid line           an element field drawn through the element's own shape
                          functions (``field``), broken at element boundaries
solid triangle            an explicit nodal value (``nodal``)
small solid circle        a Gauss point value (``gauss``)
hollow circle             a nodal value projected from Gauss points (``projected``)
filled marker, joined     a sequence of solves: load steps, speeds, lengths
by nothing                (``solves``)
black dashed line         an analytic reference (``analytic``)
hollow black diamond      a computed reference: shell, solid, finite strip (``reference``)
========================  ==========================================================

``use()`` sets the rcParams; ``save()`` writes a figure at 200 dpi. ``MODEL_COLOURS``
gives each beam theory the same colour in every example.

Public API
----------
use, save, MODEL_COLOURS, MATERIAL_COLOURS, colour, field, nodal, gauss, projected,
solves, analytic, reference, marks_legend, panel_title.
"""

from __future__ import annotations

from typing import Iterable, Optional, Sequence

import numpy as np
import numpy.typing as npt

__all__ = [
    "use", "save", "MODEL_COLOURS", "MATERIAL_COLOURS", "colour", "field", "nodal", "gauss",
    "projected", "solves", "analytic", "reference", "marks_legend", "panel_title",
]

DPI = 200
INK = "#1b1f24"
MUTED = "#6b7480"
GRID = "#d9dde2"

#: one colour per model, the same in every example (colour-blind safe)
MODEL_COLOURS = {
    "Euler-Bernoulli": "#8a8f98",
    "Timoshenko": "#0b6fb0",
    "Benscoter": "#109b72",
    "GBT": "#d4570b",
    "GBT, free root warping": "#e9a23b",
    "2D sections": "#a2427e",
    "corotational": "#0b6fb0",
    "linear": "#8a8f98",
}
_CYCLE = ["#0b6fb0", "#d4570b", "#109b72", "#a2427e", "#e9a23b", "#5a5fc4", "#8a8f98"]

#: laminate materials as they are drawn in sections
MATERIAL_COLOURS = {
    "CFRP": "#28323d",
    "carbon": "#28323d",
    "GFRP": "#5f9fca",
    "glass": "#5f9fca",
    "foam": "#e8c766",
    "gelcoat": "#c9ced4",
    "steel": "#7d8792",
    "aluminium": "#9aa5b1",
}


def colour(name: str, i: int = 0) -> str:
    """The model's colour, or the ``i``-th colour of the cycle for any other series."""
    return MODEL_COLOURS.get(name, _CYCLE[i % len(_CYCLE)])


def use() -> None:
    """Set the rcParams of the style (Agg backend, 200 dpi, quiet axes)."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from cycler import cycler

    plt.rcParams.update({
        "figure.dpi": 100,
        "savefig.dpi": DPI,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.06,
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "font.size": 9,
        "axes.titlesize": 9.5,
        "axes.titleweight": "bold",
        "axes.titlelocation": "left",
        "axes.labelsize": 9,
        "axes.labelcolor": INK,
        "axes.edgecolor": MUTED,
        "axes.linewidth": 0.7,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": True,
        "axes.axisbelow": True,
        "axes.prop_cycle": cycler(color=_CYCLE),
        "grid.color": GRID,
        "grid.linewidth": 0.5,
        "xtick.color": MUTED,
        "ytick.color": MUTED,
        "xtick.labelcolor": INK,
        "ytick.labelcolor": INK,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "legend.fontsize": 7.5,
        "legend.frameon": False,
        "legend.handlelength": 1.8,
        "lines.linewidth": 1.3,
        "lines.markersize": 4,
        "figure.titlesize": 10.5,
        "figure.titleweight": "bold",
    })


def save(fig, path: str) -> None:
    """Write ``fig`` to ``path`` at the style's resolution and close it."""
    import matplotlib.pyplot as plt

    fig.savefig(path, dpi=DPI)
    plt.close(fig)


def field(ax, x: npt.ArrayLike, y: npt.ArrayLike, color: str = INK, label: Optional[str] = None, **kw):
    """An element field through its shape functions (``NaN`` breaks between elements)."""
    return ax.plot(x, y, "-", color=color, lw=kw.pop("lw", 1.3), label=label, **kw)


def nodal(ax, x, y, color: str = INK, label: Optional[str] = None, **kw):
    """Explicit nodal values: solid triangles."""
    return ax.plot(x, y, "^", ls="none", color=color, mfc=color, mec=color, ms=kw.pop("ms", 3.6),
                   label=label, **kw)


def gauss(ax, x, y, color: str = INK, label: Optional[str] = None, **kw):
    """Gauss point values: small solid circles."""
    return ax.plot(x, y, "o", ls="none", color=color, mfc=color, mec=color, ms=kw.pop("ms", 2.6),
                   label=label, **kw)


def projected(ax, x, y, color: str = INK, label: Optional[str] = None, **kw):
    """Nodal values projected from Gauss points: hollow circles."""
    return ax.plot(x, y, "o", ls="none", color=color, mfc="white", mec=color, ms=kw.pop("ms", 4.2),
                   mew=0.9, label=label, **kw)


def solves(ax, x, y, color: str = INK, label: Optional[str] = None, marker: str = "o", **kw):
    """One value per solve (load step, speed, length): filled markers, not joined."""
    return ax.plot(x, y, marker, ls="none", color=color, mfc=color, mec="white", mew=0.5,
                   ms=kw.pop("ms", 5.0), label=label, **kw)


def analytic(ax, x, y, label: Optional[str] = "analytic", color: str = INK, **kw):
    """An analytic reference: black dashed line."""
    return ax.plot(x, y, "--", color=color, lw=kw.pop("lw", 1.0), label=label, **kw)


def reference(ax, x, y, label: Optional[str] = "reference", color: str = INK, **kw):
    """A computed reference (shell, solid, finite strip): hollow black diamonds."""
    return ax.plot(x, y, "D", ls="none", color=color, mfc="white", mec=color, mew=0.9,
                   ms=kw.pop("ms", 4.6), label=label, **kw)


_KEY = {
    "field": dict(marker="", ls="-", lw=1.3, label="element field (shape functions)"),
    "nodal": dict(marker="^", ls="none", ms=3.6, label="nodal value"),
    "gauss": dict(marker="o", ls="none", ms=2.6, label="Gauss point"),
    "projected": dict(marker="o", ls="none", ms=4.2, mfc="white", label="projected to nodes"),
    "solves": dict(marker="o", ls="none", ms=5.0, mec="white", label="one solve"),
    "analytic": dict(marker="", ls="--", lw=1.0, label="analytic"),
    "reference": dict(marker="D", ls="none", ms=4.6, mfc="white", label="reference model"),
}


def marks_legend(ax, kinds: Iterable[str], loc: str = "lower right", labels: Optional[dict] = None, **kw):
    """A grey legend of what the marks mean, added beside the series legend."""
    from matplotlib.lines import Line2D

    labels = labels or {}
    handles = []
    for k in kinds:
        spec = dict(_KEY[k])
        spec["label"] = labels.get(k, spec["label"])
        spec.setdefault("mfc", MUTED)
        handles.append(Line2D([], [], color=MUTED, mec=MUTED, **spec))
    from matplotlib.legend import Legend

    leg = Legend(ax, handles, [h.get_label() for h in handles], loc=loc, fontsize=7, labelcolor=MUTED,
                 frameon=False, **kw)
    ax.add_artist(leg)
    return leg


def panel_title(ax, title: str, note: Optional[str] = None) -> None:
    """A bold panel title with an optional grey note under it."""
    ax.set_title(title, pad=(4 + 10 * (note.count("\n") + 1)) if note else 6)
    if note:
        ax.text(0.0, 1.015, note, transform=ax.transAxes, fontsize=7.5, color=MUTED, va="bottom", ha="left")


def finite(values: Sequence[float]) -> npt.NDArray[np.float64]:
    a = np.asarray(values, dtype=np.float64)
    return a[np.isfinite(a)]

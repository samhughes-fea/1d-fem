# post_processing/graphical_visualisers/thin_walled_drawing.py
"""Drawings of a thin-walled section and of wall fields around its contour.

``draw_section`` draws every wall of a ``ThinWalledSectionAnalysis`` at its true
laminate thickness about its mid-line, ply by ply, coloured by material and shaded by
fibre angle (0° plies darkest). The nodes of the section's wall elements
(``n_sub`` cubic elements per wall) are marked on the mid-line.

``contour_field`` draws a wall field, such as ``N_xx`` or the shear flow ``N_xs``, as a
diagram normal to each wall. The value is taken at points inside every wall element
from the section's own solution (``SaintVenantRecovery.at``, which evaluates the
element's cubic shape functions), so the diagram is the finite element field and
nothing is fitted between points. The section's Gauss points are marked on it.

``unrolled_field`` draws the same field against the distance around the contour.

Public API
----------
draw_section, contour_field, unrolled_field, wall_samples, gauss_s.
"""

from __future__ import annotations

from typing import Callable, Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import numpy.typing as npt

from onedfem.post_processing.graphical_visualisers import figure_style as fs

__all__ = ["draw_section", "contour_field", "unrolled_field", "wall_samples", "gauss_s"]

Array = npt.NDArray[np.float64]
ValueAt = Callable[[int, float], float]


def _mix(c: str, w: float) -> Tuple[float, float, float]:
    """``c`` mixed with white by ``w`` (0 keeps ``c``)."""
    from matplotlib.colors import to_rgb

    r = np.array(to_rgb(c))
    return tuple(r + (1.0 - r) * w)


def _material_colour(mat, materials: Sequence[Tuple[str, object]], seen: Dict[int, str]) -> str:
    for name, m in materials:
        if m == mat:
            return fs.MATERIAL_COLOURS.get(name, fs.colour(name, len(seen)))
    key = id(mat)
    if key not in seen:
        seen[key] = fs.colour("", len(seen))
    return seen[key]


def gauss_s(analysis, wall_id: int) -> Array:
    """Distances along a wall of the section analysis's Gauss points."""
    return analysis.gauss_s(wall_id)


def wall_samples(analysis, wall_id: int, per_element: int = 9) -> List[Array]:
    """Distances along a wall, ``per_element`` points inside each wall element (ends included)."""
    sec = analysis.section
    L = sec.length(sec.wall(wall_id))
    le = L / analysis.n_sub
    return [le * (k + np.linspace(0.0, 1.0, per_element)) for k in range(analysis.n_sub)]


def draw_section(ax, analysis, materials: Sequence[Tuple[str, object]] = (), nodes: bool = True,
                 angle_shading: bool = True, edge: bool = True, walls: Optional[Iterable[int]] = None,
                 legend: bool = True, zorder: float = 2.0, ply_lines: bool = False):
    """Draw the walls at true thickness, ply by ply; returns the legend handles by material.

    ``materials`` pairs a label with a ``PlyMaterial`` (``[("CFRP", cfrp), ...]``); the label
    picks the colour from ``MATERIAL_COLOURS``.
    """
    from matplotlib.patches import Patch, Polygon

    sec = analysis.section
    wanted = set(walls) if walls is not None else None
    seen: Dict[int, str] = {}
    used: Dict[str, str] = {}
    for w in sec.walls:
        if wanted is not None and w.wall_id not in wanted:
            continue
        a, b = sec.endpoints(w)
        nrm = sec.normal(w)
        lam = analysis.laminates[w.laminate_id]
        z = lam.interfaces
        for k, ply in enumerate(lam.plies):
            c = _material_colour(ply.material, materials, seen)
            shade = 0.7 * min(abs(ply.angle_deg), 90.0) / 45.0 * 0.5 if angle_shading else 0.0
            poly = np.array([a + nrm * z[k], b + nrm * z[k], b + nrm * z[k + 1], a + nrm * z[k + 1]])
            ax.add_patch(Polygon(poly, closed=True, fc=_mix(c, shade), ec="none", zorder=zorder))
            label = next((n for n, m in materials if m == ply.material), None) or f"material {len(used) + 1}"
            used.setdefault(label, c)
        if ply_lines:
            for zz in z[1:-1]:
                p0, p1 = a + nrm * zz, b + nrm * zz
                ax.plot([p0[0], p1[0]], [p0[1], p1[1]], color="white", lw=0.25, zorder=zorder + 0.05)
        if edge:
            for zz in (z[0], z[-1]):
                p0, p1 = a + nrm * zz, b + nrm * zz
                ax.plot([p0[0], p1[0]], [p0[1], p1[1]], color=fs.INK, lw=0.35, zorder=zorder + 0.1)
        if nodes:
            t = sec.tangent(w)
            s_nodes = analysis.element_ends(w.wall_id)
            pts = np.array([a + s * t for s in s_nodes])
            ax.plot(pts[:, 0], pts[:, 1], "|", color="white", ms=3, mew=0.7, zorder=zorder + 0.2)
    ax.set_aspect("equal")
    ax.autoscale_view()
    handles = [Patch(fc=c, ec="none", label=n) for n, c in used.items()]
    if angle_shading and handles:
        handles.append(Patch(fc="none", ec="none", label="lighter: off-axis plies"))
    if legend and handles:
        ax.legend(handles=handles, loc="upper right", fontsize=7)
    return handles


def contour_field(ax, analysis, value_at: ValueAt, scale: float, per_element: int = 9,
                  pos_colour: str = "#d4570b", neg_colour: str = "#0b6fb0", alpha: float = 0.35,
                  mark_gauss: bool = True, walls: Optional[Iterable[int]] = None, zorder: float = 3.0):
    """Diagram of ``value_at(wall_id, s)`` normal to each wall, ``scale`` length per unit value.

    Positive values are drawn on the ``+n`` side of the wall (``n = x × s``) in ``pos_colour``.
    Returns ``{wall_id: (s, values)}`` of the sampled field.
    """
    from matplotlib.collections import PolyCollection

    sec = analysis.section
    wanted = set(walls) if walls is not None else None
    out: Dict[int, Tuple[Array, Array]] = {}
    quads, colours = [], []
    for w in sec.walls:
        if wanted is not None and w.wall_id not in wanted:
            continue
        a, _ = sec.endpoints(w)
        t, nrm = sec.tangent(w), sec.normal(w)
        s_all, v_all = [], []
        for s_el in wall_samples(analysis, w.wall_id, per_element):
            v = np.array([value_at(w.wall_id, float(s)) for s in s_el])
            base = a + np.outer(s_el, t)
            tip = base + np.outer(v * scale, nrm)
            for i in range(len(s_el) - 1):
                quads.append([base[i], base[i + 1], tip[i + 1], tip[i]])
                colours.append(pos_colour if v[i] + v[i + 1] >= 0.0 else neg_colour)
            ax.plot(tip[:, 0], tip[:, 1], color=fs.INK, lw=0.6, zorder=zorder + 0.1)
            s_all.append(s_el)
            v_all.append(v)
        if mark_gauss:
            sg = gauss_s(analysis, w.wall_id)
            vg = np.array([value_at(w.wall_id, float(s)) for s in sg])
            pg = a + np.outer(sg, t) + np.outer(vg * scale, nrm)
            ax.plot(pg[:, 0], pg[:, 1], "o", ms=1.6, color=fs.INK, zorder=zorder + 0.2)
        out[w.wall_id] = (np.concatenate(s_all), np.concatenate(v_all))
    ax.add_collection(PolyCollection(quads, facecolors=colours, edgecolors=colours, linewidths=0.3, alpha=alpha,
                                     zorder=zorder))
    ax.set_aspect("equal")
    ax.autoscale_view()
    return out


def unrolled_field(ax, analysis, value_at: ValueAt, wall_order: Sequence[int], color: str = fs.INK,
                   label: Optional[str] = None, per_element: int = 9, mark_gauss: bool = True,
                   wall_labels: Optional[Dict[int, str]] = None, shade: bool = True,
                   groups: Optional[Dict[str, Sequence[int]]] = None) -> Array:
    """``value_at`` against the distance around the contour, walls in ``wall_order``.

    Each wall element is drawn on its own through its shape functions; Gauss points are
    marked. Returns the contour distance at the start of each wall (and the end).
    """
    sec = analysis.section
    starts = [0.0]
    for wid in wall_order:
        starts.append(starts[-1] + sec.length(sec.wall(wid)))
    xs, vs, xg, vg = [], [], [], []
    for i, wid in enumerate(wall_order):
        for s_el in wall_samples(analysis, wid, per_element):
            xs.extend(starts[i] + s_el)
            vs.extend(value_at(wid, float(s)) for s in s_el)
            xs.append(np.nan)
            vs.append(np.nan)
        if mark_gauss:
            sg = gauss_s(analysis, wid)
            xg.extend(starts[i] + sg)
            vg.extend(value_at(wid, float(s)) for s in sg)
    fs.field(ax, xs, vs, color=color, label=label, lw=1.1)
    if mark_gauss:
        fs.gauss(ax, xg, vg, color=color, ms=1.8)
    if shade and groups:
        index = {wid: i for i, wid in enumerate(wall_order)}
        for j, (name, ws) in enumerate(groups.items()):
            ids = [index[w] for w in ws if w in index]
            if not ids:
                continue
            a, b = starts[min(ids)], starts[max(ids) + 1]
            if j % 2 == 1:
                ax.axvspan(a, b, color=fs.GRID, alpha=0.45, lw=0, zorder=0)
            ax.text(0.5 * (a + b), 0.99, name, fontsize=6.5, color=fs.MUTED, ha="center", va="top",
                    transform=ax.get_xaxis_transform())
    elif shade:
        for i, wid in enumerate(wall_order):
            if i % 2 == 1:
                ax.axvspan(starts[i], starts[i + 1], color=fs.GRID, alpha=0.35, lw=0, zorder=0)
        if wall_labels:
            for i, wid in enumerate(wall_order):
                if wid in wall_labels:
                    ax.text(0.5 * (starts[i] + starts[i + 1]), 1.0, wall_labels[wid], rotation=0, fontsize=6.5,
                            color=fs.MUTED, ha="center", va="bottom", transform=ax.get_xaxis_transform())
    return np.array(starts)

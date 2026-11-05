# post_processing/stations.py
"""Where along the beam to write section and wall stresses.

By default the stress files hold every element's mid-point. ``stress_stations.txt`` in
the job directory picks stations instead::

    [Stress_Station]
    [x]
    12.5
    30.0

A row with one value is the plane ``X = x``: the station is where the beam axis crosses it
(the first element, in element order, whose axis does). A row with three values
``[x] [y] [z]`` is the nearest point of the beam axis to that point.

An optional ``[Wall_Points]`` block picks where along each wall of a thin-walled section
``wall_stresses.csv`` is written (by default 5 evenly spaced points, or the wall ends on
sections of more than 64 walls)::

    [Wall_Points]
    nodes gauss

``nodes`` are the wall elements' end nodes, ``gauss`` the section's Gauss points (where the
section stiffness was integrated), and ``uniform n`` ``n`` evenly spaced points. Several
kinds on one line are all written; the ``point_kind`` column says which each row is. The
file may hold only this block, and then the stations are the element mid-points.

Public API
----------
read_stress_stations, read_wall_points, wall_point_list, axis_point, locate_station, element_points.
"""

from __future__ import annotations

import os
from typing import List, Optional, Sequence, Tuple

import numpy as np

FILENAME = "stress_stations.txt"


BLOCKS = ("[stress_station]", "[wall_points]")
WALL_POINT_KINDS = ("nodes", "gauss", "uniform")


def _blocks(job_dir: str) -> Optional[dict]:
    """Data lines of ``stress_stations.txt`` by block (lines before any block are stations),
    and under ``"seen"`` the blocks the file names."""
    path = os.path.join(job_dir, FILENAME)
    if not os.path.isfile(path):
        return None
    out = {b: [] for b in BLOCKS}
    out["seen"] = set()
    block = "[stress_station]"
    with open(path, encoding="utf-8") as f:
        for line in f:
            s = line.split("#", 1)[0].strip()
            if not s:
                continue
            if s.lower() in BLOCKS:
                block = s.lower()
                out["seen"].add(block)
            elif not s.startswith("["):                 # else a column header such as [x] [y] [z]
                out[block].append(s)
    return out


def read_stress_stations(job_dir: str) -> Optional[List[np.ndarray]]:
    """Stations of ``stress_stations.txt`` (each (1,) or (3,)), or ``None`` without the file or
    when it holds no stations (only a ``[Wall_Points]`` block)."""
    blocks = _blocks(job_dir)
    if blocks is None:
        return None
    out = []
    for s in blocks["[stress_station]"]:
        v = np.array([float(t) for t in s.replace(",", " ").split()], dtype=np.float64)
        if v.size not in (1, 3):
            raise ValueError(f"{FILENAME}: a station is [x] or [x] [y] [z], got {s!r}")
        out.append(v)
    if not out and "[wall_points]" in blocks["seen"]:
        return None
    return out


def read_wall_points(job_dir: str) -> Optional[List[Tuple[str, int]]]:
    """The ``[Wall_Points]`` block of ``stress_stations.txt`` as ``[(kind, n), ...]``, or ``None``.

    ``kind`` is ``"nodes"``, ``"gauss"`` or ``"uniform"``; ``n`` is the point count of
    ``uniform`` (0 for the others).
    """
    blocks = _blocks(job_dir)
    if blocks is None or not blocks["[wall_points]"]:
        return None
    words = " ".join(blocks["[wall_points]"]).replace(",", " ").split()
    spec, k = [], 0
    while k < len(words):
        w = words[k].lower()
        if w not in WALL_POINT_KINDS:
            raise ValueError(f"{FILENAME} [Wall_Points]: {words[k]!r} is not one of {', '.join(WALL_POINT_KINDS)}")
        if w == "uniform":
            if k + 1 >= len(words) or not words[k + 1].isdigit() or int(words[k + 1]) < 2:
                raise ValueError(f"{FILENAME} [Wall_Points]: 'uniform' needs a point count of 2 or more")
            spec.append(("uniform", int(words[k + 1])))
            k += 2
        else:
            spec.append((w, 0))
            k += 1
    return spec


def wall_point_list(analysis, wall_id: int, spec: Sequence[Tuple[str, int]]) -> List[Tuple[float, str]]:
    """``[(s, point_kind), ...]`` along wall ``wall_id`` for a ``read_wall_points`` spec.

    ``point_kind`` is ``"wall_node"``, ``"wall_gauss"`` or ``"uniform"``.
    """
    out: List[Tuple[float, str]] = []
    for kind, n in spec:
        if kind == "nodes":
            out += [(float(s), "wall_node") for s in analysis.element_ends(wall_id)]
        elif kind == "gauss":
            out += [(float(s), "wall_gauss") for s in analysis.gauss_s(wall_id)]
        else:
            L = analysis.section.length(analysis.section.wall(wall_id))
            out += [(float(s), "uniform") for s in np.linspace(0.0, L, n)]
    return out


def axis_point(element, xi: float) -> np.ndarray:
    """Global position (3,) of the beam axis at ``xi`` of ``element``."""
    t = 0.5 * (float(xi) + 1.0)
    if getattr(element, "_frame_mesh", False):
        p0 = np.asarray(element._frame_p0, dtype=np.float64)
        arc = getattr(element, "_frame_arc", None)
        if arc is None:
            L = float(element.grid_array[1, 0])
            return p0 + t * L * np.asarray(element._frame_R0, dtype=np.float64)[:, 0]
        p_mid = p0 - arc.position(-0.5 * arc.length, np.zeros(3))
        return arc.position((t - 0.5) * arc.length, p_mid)
    X = getattr(element, "grid_array", None)
    if X is None:
        X = np.array([[element.x_start, 0.0, 0.0], [element.x_end, 0.0, 0.0]])
    X = np.asarray(X, dtype=np.float64)
    return (1.0 - t) * X[0] + t * X[1]


def _ends(element) -> Tuple[np.ndarray, np.ndarray]:
    return axis_point(element, -1.0), axis_point(element, 1.0)


def _refine(element, xi: float, f) -> float:
    """Polish ``xi`` on a curved element, where the chord is only a first guess."""
    if getattr(element, "_frame_arc", None) is None:
        return xi
    from scipy.optimize import minimize_scalar

    r = minimize_scalar(lambda x: f(axis_point(element, x)), bounds=(-1.0, 1.0), method="bounded",
                        options={"xatol": 1e-10})
    return float(r.x)


def locate_station(elements: Sequence, station: np.ndarray) -> Tuple[int, float]:
    """``(element index, xi)`` of a station (see the module notes)."""
    station = np.asarray(station, dtype=np.float64)
    if station.size == 1:
        x = float(station[0])
        for i, el in enumerate(elements):
            a, b = _ends(el)
            lo, hi = min(a[0], b[0]), max(a[0], b[0])
            tol = 1e-9 * max(abs(hi - lo), 1.0)
            if lo - tol <= x <= hi + tol:
                xi = -1.0 if hi == lo else 2.0 * (x - a[0]) / (b[0] - a[0]) - 1.0
                xi = _refine(el, float(np.clip(xi, -1.0, 1.0)), lambda p: (p[0] - x) ** 2)
                return i, xi
        raise ValueError(f"stress station X = {x:g} is not on the beam")
    best = (np.inf, 0, 0.0)
    for i, el in enumerate(elements):
        a, b = _ends(el)
        d = b - a
        t = float(np.clip((station - a) @ d / (d @ d), 0.0, 1.0))
        xi = _refine(el, 2.0 * t - 1.0, lambda p: float((p - station) @ (p - station)))
        dist = float(np.linalg.norm(axis_point(el, xi) - station))
        if dist < best[0]:
            best = (dist, i, xi)
    return best[1], best[2]


def element_points(elements: Sequence, stations: Optional[Sequence]) -> List[Tuple[int, float]]:
    """``(element index, xi)`` to write: every element's mid-point, or the stations."""
    if stations is None:
        return [(i, 0.0) for i in range(len(elements))]
    return [locate_station(elements, s) for s in stations]

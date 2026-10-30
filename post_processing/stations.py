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

Public API
----------
read_stress_stations, axis_point, locate_station, element_points.
"""

from __future__ import annotations

import os
from typing import List, Optional, Sequence, Tuple

import numpy as np

FILENAME = "stress_stations.txt"


def read_stress_stations(job_dir: str) -> Optional[List[np.ndarray]]:
    """Stations of ``stress_stations.txt`` (each (1,) or (3,)), or ``None`` without the file."""
    path = os.path.join(job_dir, FILENAME)
    if not os.path.isfile(path):
        return None
    out = []
    with open(path, encoding="utf-8") as f:
        for line in f:
            s = line.split("#", 1)[0].strip()
            if not s or s.startswith("["):
                continue
            v = np.array([float(t) for t in s.replace(",", " ").split()], dtype=np.float64)
            if v.size not in (1, 3):
                raise ValueError(f"{FILENAME}: a station is [x] or [x] [y] [z], got {s!r}")
            out.append(v)
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

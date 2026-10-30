# post_processing/thin_walled/wall_stresses.py
"""Write wall and ply stresses of a solved job with thin-walled sections.

For every element with a thin-walled section, at its mid-point (``xi = 0``, or at the
stations of ``stress_stations.txt``, see ``post_processing/stations.py``) and at
``points_per_wall`` points along each wall, the recovered state is written to
``secondary_results/wall_stresses/wall_stresses.csv``: wall resultants and, at the
bottom and top of every ply, the wall-axis, ply-axis and interlaminar stresses.
See docs/element_library/thin_walled_stress_recovery.md.

Public API
----------
WALL_STRESS_COLUMNS, element_wall_rows, write_wall_stresses.
"""

from __future__ import annotations

import csv
import os
from typing import Iterable, List, Optional, Sequence

import numpy as np

from pre_processing.section_integration.thin_walled.element_recovery import element_wall_stress

MANY_WALLS = 64

WALL_STRESS_COLUMNS = (
    "element_id", "x", "wall_id", "s", "y", "z",
    "N_xx", "N_ss", "N_xs", "M_xx", "M_ss", "M_xs",
    "ply", "n", "sigma_xx", "sigma_ss", "tau_xs", "sigma_11", "sigma_22", "tau_12",
    "tau_xn", "tau_sn", "sigma_nn", "xi", "mode",
)


def element_wall_rows(element, U_e, analysis, points_per_wall: int = 5, xi: float = 0.0, mode="") -> List[list]:
    """CSV rows (``WALL_STRESS_COLUMNS``) for one element; ``x`` is the global ``X`` of the axis at ``xi``."""
    from post_processing.stations import axis_point

    rows = []
    x = float(axis_point(element, xi)[0])
    for w in analysis.section.walls:
        L = analysis.section.length(w)
        for s in np.linspace(0.0, L, points_per_wall):
            p = element_wall_stress(element, U_e, xi, analysis, w.wall_id, float(s), points_per_ply=2)
            for q in p.plies:
                rows.append([
                    element.element_id, x, w.wall_id, float(s), *p.position, *p.resultants,
                    q.ply, q.n, *q.stress_wall, *q.stress_ply, q.tau_xn, q.tau_sn, q.sigma_nn, xi, mode,
                ])
    return rows


def write_wall_stresses(
    job_results_dir: str,
    elements: Sequence,
    analyses: Sequence,
    U_global: np.ndarray,
    points_per_wall: Optional[int] = None,
    stations: Optional[Sequence] = None,
    modes: Optional[np.ndarray] = None,
) -> str:
    """Write ``wall_stresses.csv`` for the elements that have a thin-walled analysis; return its path.

    ``points_per_wall`` defaults to 5, or to 2 (the wall ends) on sections of more than
    ``MANY_WALLS`` walls, such as a blade cut whose curved skins are many short walls.
    ``stations`` (``post_processing.stations``) replaces the element mid-points. ``modes``
    (n_dof, n_modes) writes one block of rows per mode shape instead of ``U_global``.
    """
    from post_processing.stations import element_points

    out_dir = os.path.join(job_results_dir, "secondary_results", "wall_stresses")
    os.makedirs(out_dir, exist_ok=True)
    path = os.path.join(out_dir, "wall_stresses.csv")
    states = [("", np.asarray(U_global, dtype=np.float64).ravel())] if modes is None else \
        [(k + 1, np.asarray(modes, dtype=np.float64)[:, k]) for k in range(np.asarray(modes).shape[1])]
    points = element_points(elements, stations)
    with open(path, "w", newline="", encoding="utf-8") as f:
        wr = csv.writer(f)
        wr.writerow(WALL_STRESS_COLUMNS)
        for mode, U in states:
            for i, xi in points:
                el, an = elements[i], analyses[i]
                if an is None:
                    continue
                U_e = U[el.assemble_global_dof_indices()]
                n = points_per_wall or (2 if len(an.section.walls) > MANY_WALLS else 5)
                try:
                    wr.writerows(element_wall_rows(el, U_e, an, n, xi=xi, mode=mode))
                except NotImplementedError:
                    continue  # element type without wall recovery
    return path

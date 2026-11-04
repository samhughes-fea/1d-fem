# post_processing/solid/section_stresses.py
"""Write 3D section stresses of a solved job with sections meshed in 2D.

For every element with a 2D section, at its mid-point (``xi = 0``, or at the stations of
``stress_stations.txt``, see ``post_processing/stations.py``), the 3D stresses at every
quadrature point of the section are written to
``secondary_results/section_stresses/section_stresses.csv``: in section axes, in material
(ply or core) axes, and the interlaminar stresses in wall axes. The interlaminar stresses
on the ply interfaces, recovered by ``SolidSectionAnalysis.interface_stresses``, go to
``interface_stresses.csv`` next to it.
See docs/element_library/solid_section_analysis.md.

Public API
----------
SECTION_STRESS_COLUMNS, INTERFACE_STRESS_COLUMNS, element_section_rows, element_interface_rows,
write_section_stresses.
"""

from __future__ import annotations

import csv
import os
from typing import List, Optional, Sequence

import numpy as np

from pre_processing.section_integration.solid.element_recovery import element_section_state, element_section_stresses

SECTION_STRESS_COLUMNS = (
    "element_id", "x", "section_element", "y", "z",
    "sigma_xx", "sigma_yy", "sigma_zz", "tau_yz", "tau_xz", "tau_xy",
    "sigma_11", "sigma_22", "sigma_33", "tau_23", "tau_13", "tau_12",
    "sigma_nn", "tau_sn", "tau_xn", "xi", "mode",
)

INTERFACE_STRESS_COLUMNS = (
    "element_id", "x", "below", "above", "y", "z", "sigma_nn", "tau_sn", "tau_xn", "xi", "mode",
)


def _x(element, xi: float) -> float:
    from post_processing.stations import axis_point

    return float(axis_point(element, xi)[0])


def element_section_rows(element, U_e, analysis, xi: float = 0.0, mode="") -> List[list]:
    """CSV rows (``SECTION_STRESS_COLUMNS``) for one element; ``x`` is the global ``X`` of the axis at ``xi``."""
    pts = element_section_stresses(element, U_e, xi, analysis)
    x = _x(element, xi)
    return [[element.element_id, x, p.element, p.y, p.z, *p.stress, *p.stress_material, *p.interlaminar, xi, mode]
            for p in pts]


def element_interface_rows(element, U_e, analysis, xi: float = 0.0, mode="") -> List[list]:
    """CSV rows (``INTERFACE_STRESS_COLUMNS``): the interlaminar stresses on every ply interface of
    the section (``SolidSectionAnalysis.interface_stresses``), at the middle of each interface edge."""
    forces, benscoter = element_section_state(element, U_e, xi, analysis)
    x = _x(element, xi)
    return [[element.element_id, x, e, f, y, z, *s, xi, mode]
            for e, f, y, z, s in analysis.interface_stresses(forces=forces, benscoter=benscoter)]


def write_section_stresses(job_results_dir: str, elements: Sequence, analyses: Sequence, U_global,
                           stations: Optional[Sequence] = None, modes: Optional[np.ndarray] = None) -> str:
    """Write ``section_stresses.csv`` and ``interface_stresses.csv`` for the elements that have a 2D
    section; return the path of the first.

    ``stations`` (``post_processing.stations``) replaces the element mid-points. ``modes``
    (n_dof, n_modes) writes one block of rows per mode shape instead of ``U_global``.
    """
    from post_processing.stations import element_points

    out_dir = os.path.join(job_results_dir, "secondary_results", "section_stresses")
    os.makedirs(out_dir, exist_ok=True)
    path = os.path.join(out_dir, "section_stresses.csv")
    states = [("", np.asarray(U_global, dtype=np.float64).ravel())] if modes is None else \
        [(k + 1, np.asarray(modes, dtype=np.float64)[:, k]) for k in range(np.asarray(modes).shape[1])]
    points = element_points(elements, stations)
    with open(path, "w", newline="", encoding="utf-8") as f, \
            open(os.path.join(out_dir, "interface_stresses.csv"), "w", newline="", encoding="utf-8") as g:
        wr, wi = csv.writer(f), csv.writer(g)
        wr.writerow(SECTION_STRESS_COLUMNS)
        wi.writerow(INTERFACE_STRESS_COLUMNS)
        for mode, U in states:
            for i, xi in points:
                el, an = elements[i], analyses[i]
                if an is None:
                    continue
                U_e = U[el.assemble_global_dof_indices()]
                try:
                    wr.writerows(element_section_rows(el, U_e, an, xi=xi, mode=mode))
                    if hasattr(an, "interface_stresses"):
                        wi.writerows(element_interface_rows(el, U_e, an, xi=xi, mode=mode))
                except NotImplementedError:
                    continue  # element type without section recovery
    return path

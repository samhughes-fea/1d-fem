# post_processing/solid/section_stresses.py
"""Write 3D section stresses of a solved job with sections meshed in 2D.

For every element with a 2D section, at its mid-point (``xi = 0``), the 3D stresses at
every quadrature point of the section are written to
``secondary_results/section_stresses/section_stresses.csv``: in section axes, in material
(ply or core) axes, and the interlaminar stresses in wall axes.
See docs/element_library/solid_section_analysis.md.

Public API
----------
SECTION_STRESS_COLUMNS, element_section_rows, write_section_stresses.
"""

from __future__ import annotations

import csv
import os
from typing import List, Sequence

import numpy as np

from pre_processing.section_integration.solid.element_recovery import element_section_stresses

SECTION_STRESS_COLUMNS = (
    "element_id", "x", "section_element", "y", "z",
    "sigma_xx", "sigma_yy", "sigma_zz", "tau_yz", "tau_xz", "tau_xy",
    "sigma_11", "sigma_22", "sigma_33", "tau_23", "tau_13", "tau_12",
    "sigma_nn", "tau_sn", "tau_xn",
)


def element_section_rows(element, U_e, analysis, xi: float = 0.0) -> List[list]:
    """CSV rows (``SECTION_STRESS_COLUMNS``) for one element."""
    x = 0.5 * (1 - xi) * element.x_start + 0.5 * (1 + xi) * element.x_end
    return [[element.element_id, x, p.element, p.y, p.z, *p.stress, *p.stress_material, *p.interlaminar]
            for p in element_section_stresses(element, U_e, xi, analysis)]


def write_section_stresses(job_results_dir: str, elements: Sequence, analyses: Sequence, U_global) -> str:
    """Write ``section_stresses.csv`` for the elements that have a 2D section; return its path."""
    out_dir = os.path.join(job_results_dir, "secondary_results", "section_stresses")
    os.makedirs(out_dir, exist_ok=True)
    path = os.path.join(out_dir, "section_stresses.csv")
    U = np.asarray(U_global, dtype=np.float64).ravel()
    with open(path, "w", newline="", encoding="utf-8") as f:
        wr = csv.writer(f)
        wr.writerow(SECTION_STRESS_COLUMNS)
        for el, an in zip(elements, analyses):
            if an is None:
                continue
            try:
                wr.writerows(element_section_rows(el, U[el.assemble_global_dof_indices()], an))
            except NotImplementedError:
                continue  # element type without section recovery
    return path

# pre_processing/parsing/precurvature_parser.py
"""Parse optional per-element reference curvature / twist for straight beam elements."""

from __future__ import annotations

import logging
import os
import re
from typing import Any, Dict, Optional

import numpy as np
import numpy.typing as npt

logger = logging.getLogger(__name__)

SECTION_PATTERN = re.compile(r"^\[Precurvature\]$", re.IGNORECASE)


def parse_precurvature(
    file_path: Optional[str],
    element_ids: npt.NDArray[np.int64],
) -> npt.NDArray[np.float64]:
    """
    Read ``precurvature.txt`` and return ``(N_e, 3)`` columns ``[k_x0, k_y0, k_z0]`` (1/m),
    row order matching *element_ids* (same order as ``element_dictionary['ids']``).

    If *file_path* is missing or the file does not exist, returns zeros.

    File format::

        [Precurvature]
        [element_id] [k_x0] [k_y0] [k_z0]
        1  0.0  0.0  0.0

    Lines starting with ``#`` and blank lines are ignored. Inline ``#`` comments are stripped.
    """
    n_e = int(element_ids.size)
    out = np.zeros((n_e, 3), dtype=np.float64)
    if not file_path or not os.path.isfile(file_path):
        return out

    id_to_row: Dict[int, npt.NDArray[np.float64]] = {}
    in_section = False
    header_skipped = False

    with open(file_path, encoding="utf-8") as fh:
        for line_no, raw in enumerate(fh, 1):
            line = raw.split("#")[0].strip()
            if not line:
                continue
            if SECTION_PATTERN.match(line):
                in_section = True
                continue
            if not in_section:
                continue

            parts = line.split()
            if not header_skipped:
                # Skip subheader row like [element_id] [k_x0] ...
                if any(p.startswith("[") for p in parts):
                    header_skipped = True
                    continue
                # First data row (no bracket subheader in file)
                header_skipped = True

            if len(parts) != 4:
                logger.warning(
                    "precurvature.txt line %d: expected 4 tokens (element_id k_x0 k_y0 k_z0), got %d; skipping",
                    line_no,
                    len(parts),
                )
                continue
            try:
                eid = int(parts[0])
                kx0, ky0, kz0 = float(parts[1]), float(parts[2]), float(parts[3])
            except ValueError:
                logger.warning("precurvature.txt line %d: invalid numbers; skipping", line_no)
                continue
            id_to_row[eid] = np.array([kx0, ky0, kz0], dtype=np.float64)

    # Fill rows in element_ids order
    ids_list = element_ids.tolist()
    missing: list[int] = []
    for i, eid in enumerate(ids_list):
        row = id_to_row.get(int(eid))
        if row is None:
            missing.append(int(eid))
            continue
        out[i, :] = row

    if missing:
        logger.info(
            "precurvature.txt: %d element id(s) had no row (using zeros): %s",
            len(missing),
            missing[:20] + (["..."] if len(missing) > 20 else []),
        )

    extra = set(id_to_row.keys()) - set(ids_list)
    if extra:
        logger.warning("precurvature.txt: ignoring %d unknown element id(s)", len(extra))

    return out


def element_reference_strain_voigt(
    element_dictionary: Dict[str, Any],
    element_id: int,
) -> npt.NDArray[np.float64]:
    """
    Return ``(6,)`` Voigt ``E_0`` for *element_id* from ``precurvature_per_element`` if present, else zeros.
    """
    idx = int(np.where(np.asarray(element_dictionary["ids"], dtype=np.int64) == int(element_id))[0][0])
    prec = element_dictionary.get("precurvature_per_element")
    if prec is None:
        return np.zeros(6, dtype=np.float64)
    row = np.asarray(prec, dtype=np.float64)[idx]
    return reference_strain_voigt(row)


def reference_strain_voigt(k_xyz: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
    """
    Build 6-vector ``E_0`` matching Voigt ``[eps_x, kappa_y, kappa_z, gamma_xy, gamma_xz, phi_x]``.

    Parameters
    ----------
    k_xyz
        ``(3,)`` with ``[k_x0, k_y0, k_z0]`` (twist rate, bending curvatures).
    """
    k = np.asarray(k_xyz, dtype=np.float64).reshape(3)
    kx0, ky0, kz0 = float(k[0]), float(k[1]), float(k[2])
    return np.array([0.0, ky0, kz0, 0.0, 0.0, kx0], dtype=np.float64)

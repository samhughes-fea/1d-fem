# post_processing/failure/section_force_checks.py
"""Ply checks under section forces given from outside the job, and their envelope.

With ``section_forces.txt`` (``pre_processing/parsing/section_forces_parser.py``) and
``ply_strength.txt`` in the job, every row's forces are applied to the section of its element
and the plies checked by the section provider's ``SectionStrength``: the four criteria of
``ply_failure`` at sampling points over the section (thin-walled: ``points_per_wall`` along
every wall and the bottom and top of every ply; 2D: every quadrature point). Each section is
solved once for six unit forces, so a load case costs a matrix product and the criteria.

Two files go to ``secondary_results/section_force_checks/``:

- ``section_force_checks.csv``: for every row and criterion (and ``exposure``, the largest of
  them), the most exposed point: exposure, margin, mode, wall, ply, height through the wall
  ``n``, position ``(y, z)`` in the section, 2D element and material;
- ``section_force_envelope.csv``: for every element and station ``x`` and criterion, the
  largest exposure over all load cases, the load case that sets it and where.

The exposure ``f`` scales with the load, the reserve factor is ``1 / f`` and the margin
``1 / f - 1`` (negative: failed). Table sections (``section_matrix.txt``) have no plies and
cannot be checked.

Public API
----------
CHECK_COLUMNS, ENVELOPE_COLUMNS, section_strength_checks, write_section_force_checks.
"""

from __future__ import annotations

import csv
import os
from typing import Dict, List, Mapping, Optional, Sequence

import numpy as np

from onedfem.sections import section_import

CHECK_COLUMNS = ("load_case", "element_id", "x", "criterion", "exposure", "margin", "mode", "wall_id", "ply", "n",
                 "y", "z", "section_element", "material_id")
ENVELOPE_COLUMNS = ("element_id", "x", "criterion", "exposure", "margin", "mode", "load_case", "wall_id", "ply", "n",
                    "y", "z", "section_element", "material_id")


def _material_lookup(job_dir: str, strengths):
    """``PlyMaterial -> PlyStrength`` and ``PlyMaterial -> material id`` for thin-walled plies."""
    from onedfem.post_processing.failure.ply_failure import _ply_material_keys

    keys = _ply_material_keys(job_dir, strengths)
    return (lambda mat: keys[mat][1] if mat in keys else None), {mat: mid for mat, (mid, _) in keys.items()}


def section_strength_checks(job_dir: str, element_ids: Sequence[int], thin_walled: Optional[Sequence],
                            solid: Optional[Sequence], strengths: Mapping[int, object],
                            points_per_wall: Optional[int] = None) -> Dict[int, tuple]:
    """``{element_id: (SectionStrength, material id per sampling point)}`` for every element with a
    thin-walled or 2D section; elements sharing a section share its check."""
    SectionStrength = section_import("section_strength", "SectionStrength")
    from onedfem.post_processing.thin_walled.wall_stresses import MANY_WALLS

    lookup, mid_of = _material_lookup(job_dir, strengths)
    cache: Dict[int, tuple] = {}
    out: Dict[int, tuple] = {}
    for i, eid in enumerate(element_ids):
        tw = thin_walled[i] if thin_walled is not None else None
        so = solid[i] if solid is not None else None
        analysis = tw if tw is not None else so
        if analysis is None:
            continue
        if id(analysis) not in cache:
            if tw is not None:
                n = points_per_wall or (2 if len(tw.section.walls) > MANY_WALLS else 5)
                check = SectionStrength.from_thin_walled(tw, lookup, points_per_wall=n, points_per_ply=2)
                mids = []
                for w, p in zip(check.wall, check.ply):
                    lam = tw.laminates[tw.section.wall(int(w)).laminate_id]
                    mids.append(mid_of.get(lam.plies[int(p)].material, ""))
            else:
                check = SectionStrength.from_solid(so, dict(strengths))
                mids = [int(so.mesh.material_ids[int(e)]) for e in check.element]
            cache[id(analysis)] = (analysis, check, mids)
        out[int(eid)] = cache[id(analysis)][1:]
    return out


def _row(point, mids):
    i = point.index
    return [point.wall if point.wall is not None else "", point.ply, point.n if point.wall is not None else "",
            point.position[0], point.position[1], point.element if point.element is not None else "", mids[i]]


def write_section_force_checks(job_dir: str, job_results_dir: str, elements: Sequence,
                               thin_walled: Optional[Sequence], solid: Optional[Sequence],
                               strengths: Optional[Mapping[int, object]] = None,
                               points_per_wall: Optional[int] = None) -> Optional[str]:
    """Write ``section_force_checks.csv`` and ``section_force_envelope.csv``; return the first path,
    or ``None`` when the job has no ``section_forces.txt``."""
    from onedfem.pre_processing.parsing.ply_strength_parser import PLY_STRENGTH_FILE, parse_ply_strength_file
    from onedfem.pre_processing.parsing.section_forces_parser import (
        SECTION_FORCES_FILE,
        parse_section_forces_file,
        to_engine_forces,
    )

    forces_path = os.path.join(job_dir, SECTION_FORCES_FILE)
    if not os.path.exists(forces_path):
        return None
    rows = parse_section_forces_file(forces_path)
    if strengths is None:
        strengths = parse_ply_strength_file(os.path.join(job_dir, PLY_STRENGTH_FILE))
    if not strengths:
        raise ValueError(f"{SECTION_FORCES_FILE} needs {PLY_STRENGTH_FILE}: the plies to check have no strengths")
    CRITERIA = section_import("ply_failure", "CRITERIA")
    ids = [int(el.element_id) for el in elements]
    unknown = sorted({eid for _, eid, *_ in rows} - set(ids))
    if unknown:
        raise ValueError(f"{SECTION_FORCES_FILE} names unknown elements {unknown}")
    checks = section_strength_checks(job_dir, ids, thin_walled, solid, strengths, points_per_wall)
    unchecked = sorted({eid for _, eid, *_ in rows} - set(checks))
    if unchecked:
        raise ValueError(f"{SECTION_FORCES_FILE}: elements {unchecked} have no thin-walled or 2D section with plies "
                         "to check (table sections carry no plies)")

    criteria = ("exposure",) + tuple(CRITERIA)
    out_dir = os.path.join(job_results_dir, "secondary_results", "section_force_checks")
    os.makedirs(out_dir, exist_ok=True)
    path = os.path.join(out_dir, "section_force_checks.csv")
    envelope: Dict[tuple, list] = {}
    with open(path, "w", newline="", encoding="utf-8") as f:
        wr = csv.writer(f)
        wr.writerow(CHECK_COLUMNS)
        for case, eid, x, frame, forces in rows:
            check, mids = checks[eid]
            result = check.evaluate(to_engine_forces(forces, frame), order="analysis")
            for c in criteria:
                p = result.critical(c)
                margin = p.margin
                wr.writerow([case, eid, x, c, p.exposure, margin, p.mode, *_row(p, mids)])
                key = (eid, x, c)
                if key not in envelope or p.exposure > envelope[key][0]:
                    envelope[key] = [p.exposure, margin, p.mode, case, *_row(p, mids)]
    with open(os.path.join(out_dir, "section_force_envelope.csv"), "w", newline="", encoding="utf-8") as f:
        wr = csv.writer(f)
        wr.writerow(ENVELOPE_COLUMNS)
        order = {c: k for k, c in enumerate(criteria)}
        for (eid, x, c) in sorted(envelope, key=lambda k: (k[0], k[1], order[k[2]])):
            wr.writerow([eid, x, c, *envelope[(eid, x, c)]])
    return path

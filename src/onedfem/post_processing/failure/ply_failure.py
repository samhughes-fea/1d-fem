# post_processing/failure/ply_failure.py
"""Ply failure exposures and margins at the beam elements' Gauss points.

With ``ply_strength.txt`` in the job (``pre_processing/parsing/ply_strength_parser.py``),
a static job writes ``secondary_results/ply_failure/ply_failure.csv``. For every element
with a thin-walled or 2D section, at each Gauss point of its stiffness integration (the
``xi`` of the cached element object's ``gauss_data``), the section stresses are recovered
and every ply point with a strength is checked by the four criteria of
``pre_processing/section_integration/ply_failure.py``:

- thin-walled sections: the bottom and top of every ply, at ``points_per_wall`` points along
  each wall (as ``wall_stresses.csv``), with the transverse shears turned to ply axes;
- 2D sections: every quadrature point of the section mesh, in material axes.

``ply_failure_summary.csv`` gives the smallest margin of each element and criterion and
where it is. The exposure ``f`` scales with the load, the reserve factor is ``1 / f`` and
the margin ``1 / f - 1`` (negative: failed).

Public API
----------
PLY_FAILURE_COLUMNS, SUMMARY_COLUMNS, write_ply_failure.
"""

from __future__ import annotations

import csv
import os
from typing import Dict, List, Mapping, Optional, Sequence

import numpy as np

from onedfem.pre_processing.section_integration.ply_failure import CRITERIA, PlyStrength, exposures, ply_stress_from_wall

PLY_FAILURE_COLUMNS = (
    "element_id", "gauss_point", "xi", "x", "wall_id", "s", "section_element", "ply", "n", "y", "z",
    "material_id", "sigma_11", "sigma_22", "sigma_33", "tau_23", "tau_13", "tau_12",
    "max_stress", "max_stress_mode", "tsai_wu", "hashin", "hashin_mode", "puck", "puck_mode", "puck_angle",
    "exposure", "margin", "critical",
)
SUMMARY_COLUMNS = ("element_id", "criterion", "exposure", "margin", "mode", "gauss_point", "xi", "x", "wall_id", "s",
                   "section_element", "ply", "material_id")


def _ply_material_keys(job_dir: str, strengths: Mapping[int, PlyStrength]) -> Dict[object, tuple]:
    """``{PlyMaterial: (material_id, PlyStrength)}`` from ``laminate.txt`` for the thin-walled sections."""
    path = os.path.join(job_dir, "laminate.txt")
    if not os.path.exists(path):
        return {}
    from onedfem.pre_processing.parsing.laminate_parser import parse_ply_material_file

    out: Dict[object, tuple] = {}
    for mid, mat in parse_ply_material_file(path).items():
        if mid not in strengths:
            continue
        if mat in out and out[mat][1] != strengths[mid]:
            raise ValueError(f"ply_strength.txt: materials {out[mat][0]} and {mid} have the same properties but "
                             "different strengths, so a thin-walled ply cannot tell them apart")
        out.setdefault(mat, (mid, strengths[mid]))
    return out


def _gauss_xi(element_object) -> List[float]:
    data = getattr(element_object, "gauss_data", None) or []
    return [float(g.xi) for g in data]


def _thin_walled_points(element, U_e, xi, analysis, keys, points_per_wall):
    from onedfem.pre_processing.section_integration.thin_walled.element_recovery import element_wall_stress

    for w in analysis.section.walls:
        lam = analysis.laminates[w.laminate_id]
        L = analysis.section.length(w)
        for s in np.linspace(0.0, L, points_per_wall):
            p = element_wall_stress(element, U_e, xi, analysis, w.wall_id, float(s), points_per_ply=2)
            for q in p.plies:
                ply = lam.plies[q.ply]
                key = keys.get(ply.material)
                if key is None:
                    continue
                stress = ply_stress_from_wall(*q.stress_ply, q.tau_xn, q.tau_sn, q.sigma_nn, ply.angle_deg)
                yield (w.wall_id, float(s), "", q.ply, q.n, *p.position), key, stress


def _solid_points(element, U_e, xi, analysis, strengths):
    from onedfem.pre_processing.section_integration.solid.element_recovery import element_section_stresses

    mesh = analysis.mesh
    for p in element_section_stresses(element, U_e, xi, analysis):
        mid = int(mesh.material_ids[p.element])
        st = strengths.get(mid)
        if st is None:
            continue
        yield ("", "", p.element, "", "", p.y, p.z), (mid, st), np.asarray(p.stress_material, dtype=np.float64)


def write_ply_failure(job_dir: str, job_results_dir: str, elements: Sequence, element_objects: Sequence,
                      U_global, thin_walled: Optional[Sequence] = None, solid: Optional[Sequence] = None,
                      strengths: Optional[Mapping[int, PlyStrength]] = None,
                      points_per_wall: Optional[int] = None) -> Optional[str]:
    """Write ``ply_failure.csv`` and ``ply_failure_summary.csv``; return the first path, or
    ``None`` when the job has no ``ply_strength.txt`` (and no ``strengths`` are given)."""
    from onedfem.post_processing.stations import axis_point
    from onedfem.post_processing.thin_walled.wall_stresses import MANY_WALLS

    if strengths is None:
        from onedfem.pre_processing.parsing.ply_strength_parser import PLY_STRENGTH_FILE, parse_ply_strength_file

        strengths = parse_ply_strength_file(os.path.join(job_dir, PLY_STRENGTH_FILE))
    if not strengths:
        return None
    keys = _ply_material_keys(job_dir, strengths) if thin_walled is not None else {}
    U = np.asarray(U_global, dtype=np.float64).ravel()
    out_dir = os.path.join(job_results_dir, "secondary_results", "ply_failure")
    os.makedirs(out_dir, exist_ok=True)
    path = os.path.join(out_dir, "ply_failure.csv")
    summary: List[list] = []
    with open(path, "w", newline="", encoding="utf-8") as f:
        wr = csv.writer(f)
        wr.writerow(PLY_FAILURE_COLUMNS)
        for i, el in enumerate(elements):
            tw = thin_walled[i] if thin_walled is not None else None
            so = solid[i] if solid is not None else None
            if tw is None and so is None:
                continue
            U_e = U[el.assemble_global_dof_indices()]
            locs, ids, stresses, where = [], [], [], []
            try:
                for g, xi in enumerate(_gauss_xi(element_objects[i])):
                    x = float(axis_point(el, xi)[0])
                    if tw is not None:
                        n = points_per_wall or (2 if len(tw.section.walls) > MANY_WALLS else 5)
                        pts = _thin_walled_points(el, U_e, xi, tw, keys, n)
                    else:
                        pts = _solid_points(el, U_e, xi, so, strengths)
                    for loc, (mid, st), s in pts:
                        locs.append((g, xi, x, *loc))
                        ids.append((mid, st))
                        stresses.append(s)
            except NotImplementedError:
                continue  # element type without section recovery
            if not stresses:
                continue
            S = np.array(stresses)
            res = {c: np.zeros(len(S)) for c in CRITERIA}
            modes = {k: np.empty(len(S), dtype=object) for k in ("max_stress_mode", "hashin_mode", "puck_mode")}
            angle = np.zeros(len(S))
            # evaluate per material, vectorised
            by_mid: Dict[int, list] = {}
            for k, (mid, _) in enumerate(ids):
                by_mid.setdefault(mid, []).append(k)
            for mid, rows in by_mid.items():
                rows = np.array(rows)
                e = exposures(S[rows], ids[rows[0]][1])
                for c in CRITERIA:
                    res[c][rows] = e[c]
                for k in modes:
                    modes[k][rows] = e[k]
                angle[rows] = e["puck_angle"]
            F = np.stack([res[c] for c in CRITERIA])
            f = F.max(axis=0)
            crit = np.array(CRITERIA)[F.argmax(axis=0)]
            with np.errstate(divide="ignore"):
                margin = np.where(f > 0, 1.0 / np.where(f > 0, f, 1.0) - 1.0, np.inf)
            for k in range(len(S)):
                g, xi, x, wall, s, sec_el, ply, n, y, z = locs[k]
                wr.writerow([el.element_id, g, xi, x, wall, s, sec_el, ply, n, y, z, ids[k][0], *S[k],
                             res["max_stress"][k], modes["max_stress_mode"][k], res["tsai_wu"][k],
                             res["hashin"][k], modes["hashin_mode"][k], res["puck"][k], modes["puck_mode"][k],
                             angle[k], f[k], margin[k], crit[k]])
            for c in CRITERIA:
                k = int(np.argmax(res[c]))
                g, xi, x, wall, s, sec_el, ply, n, y, z = locs[k]
                mode = {"max_stress": modes["max_stress_mode"][k], "hashin": modes["hashin_mode"][k],
                        "puck": modes["puck_mode"][k]}.get(c, "")
                fc = res[c][k]
                summary.append([el.element_id, c, fc, (1.0 / fc - 1.0) if fc > 0 else np.inf, mode, g, xi, x, wall, s,
                                sec_el, ply, ids[k][0]])
    with open(os.path.join(out_dir, "ply_failure_summary.csv"), "w", newline="", encoding="utf-8") as f:
        wr = csv.writer(f)
        wr.writerow(SUMMARY_COLUMNS)
        wr.writerows(summary)
    return path

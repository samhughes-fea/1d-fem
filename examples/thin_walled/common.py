# examples/thin_walled/common.py
"""Shared helpers for the thin-walled section examples.

Each example describes a cantilever with a thin-walled section in Python, writes
its 1d-fem job files, runs it with ``process_job`` and reads back the
displacements and ``wall_stresses.csv``.
"""

from __future__ import annotations

import csv
import os
import shutil
import sys
from dataclasses import dataclass, field
from typing import Dict, List, Sequence, Tuple

import numpy as np

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
if REPO not in sys.path:
    sys.path.insert(0, REPO)

STANDARD_DOFS = ["UX", "UY", "UZ", "RX", "RY", "RZ"]


@dataclass
class Material:
    E1: float
    E2: float
    G12: float
    nu12: float
    G13: float
    G23: float
    rho: float


@dataclass
class SectionSpec:
    """Points ``{id: (y, z)}``, walls ``[(wall_id, start, end, laminate_id)]`` and laminates."""

    points: Dict[int, Tuple[float, float]]
    walls: List[Tuple[int, int, int, int]]
    materials: Dict[int, Material]
    laminates: Dict[int, List[Tuple[int, float, float]]]  # laminate_id -> [(material_id, angle, t)]
    n_sub: int = 1


@dataclass
class CantileverSpec:
    """A straight cantilever along ``x``, clamped at ``x = 0`` with all extra DOFs fixed."""

    length: float
    n_el: int
    element_type: str
    section: SectionSpec
    tip_load: Sequence[float]  # [F_x, F_y, F_z, M_x, M_y, M_z]
    modes: int = 0
    shear_modes: int = 0
    warping: int = 0
    free_root_extras: Sequence[int] = field(default_factory=tuple)  # extra DOF numbers (1-based Q) left free


def write_job(job_dir: str, spec: CantileverSpec) -> None:
    """Write the 1d-fem input files of ``spec`` into ``job_dir`` (replaced)."""
    if os.path.isdir(job_dir):
        shutil.rmtree(job_dir)
    os.makedirs(job_dir)
    n, L = spec.n_el, spec.length

    def w(name, text):
        with open(os.path.join(job_dir, name), "w", encoding="utf-8") as f:
            f.write(text)

    w("grid.txt", "[Grid]\n[node_id] [x] [y] [z]\n" + "".join(f"{i} {L * i / n:.9g} 0 0\n" for i in range(n + 1)))
    w("element.txt",
      "[Element]\n[element_id] [node1] [node2] [element_type] [axial_order] [bending_y_order] "
      "[bending_z_order] [shear_y_order] [shear_z_order] [torsion_order] [load_order] [warping]\n"
      + "".join(f"{e} {e} {e + 1} {spec.element_type} 2 2 2 2 2 2 2 {spec.warping}\n" for e in range(n)))
    # classical inputs are required by the parsers but unused with a thin-walled section
    w("material.txt", "[Material]\n[element_id] [E] [G] [nu] [rho]\n" + "".join(f"{e} 7e10 2.7e10 0.3 0\n" for e in range(n)))
    w("section.txt", "[Section]\n[element_id] [A] [I_x] [I_y] [I_z] [J_t] [kappa] [alpha]\n"
      + "".join(f"{e} 1e-3 0.0 1e-8 1e-8 1e-9 0.8333 0.0\n" for e in range(n)))
    w("point_load.txt", "[Point load]\n[x] [y] [z] [F_x] [F_y] [F_z] [M_x] [M_y] [M_z]\n"
      + f"{L} 0 0 " + " ".join(f"{v:.9g}" for v in spec.tip_load) + "\n")
    w("simulation_settings.txt",
      "[Simulation]\n[Type]\nstatic\n\n[Solver]\ntype = direct\ntolerance = 1e-12\nmax_iterations = 500\n")

    sec = spec.section
    mats = "".join(f"{i} {m.E1:.9g} {m.E2:.9g} {m.G12:.9g} {m.nu12:.9g} {m.G13:.9g} {m.G23:.9g} {m.rho:.9g}\n"
                   for i, m in sec.materials.items())
    lams = "".join(f"{lid} {k + 1} {mid} {ang:.9g} {t:.9g}\n"
                   for lid, plies in sec.laminates.items() for k, (mid, ang, t) in enumerate(plies))
    w("laminate.txt", "[Ply_Material]\n[material_id] [E1] [E2] [G12] [nu12] [G13] [G23] [rho]\n" + mats
      + "[Laminate]\n[laminate_id] [ply] [material_id] [angle_deg] [thickness]\n" + lams)
    w("section_geometry.txt",
      "[Section_Point]\n[section_id] [point_id] [y] [z]\n"
      + "".join(f"1 {p} {y:.9g} {z:.9g}\n" for p, (y, z) in sec.points.items())
      + "[Section_Wall]\n[section_id] [wall_id] [point_start] [point_end] [laminate_id]\n"
      + "".join(f"1 {wid} {a} {b} {lid}\n" for wid, a, b, lid in sec.walls)
      + f"[Section_Mesh]\n[section_id] [n_sub]\n1 {sec.n_sub}\n")
    cols, row = "[element_id] [section_id]", "1"
    if spec.modes or spec.shear_modes:
        cols += " [modes] [shear_modes]"
        row += f" {spec.modes} {spec.shear_modes}"
    w("element_section.txt", f"[Element_Section]\n{cols}\n" + "".join(f"{e} {row}\n" for e in range(n)))


def n_extra_dofs(spec: CantileverSpec) -> int:
    """Extra nodal DOFs of the job (from the section, as ``attach_section_analyses`` sets them)."""
    if "GBT" not in spec.element_type:
        return 1 if spec.warping else 0
    from pre_processing.section_integration.thin_walled.cross_section_modes import CrossSectionModes

    return CrossSectionModes(build_analysis(spec.section)).n_sv_warping + spec.shear_modes + 2 * spec.modes


def write_supports(job_dir: str, spec: CantileverSpec) -> None:
    dofs = STANDARD_DOFS + [f"Q{k}" for k in range(1, n_extra_dofs(spec) + 1) if k not in spec.free_root_extras]
    with open(os.path.join(job_dir, "prescribed_displacement.txt"), "w", encoding="utf-8") as f:
        f.write("[Prescribed Displacement]\n[id] [node_id] [dof] [value] [type] [comment]\n")
        f.writelines(f"{i} 0 {d} 0.0 displacement\n" for i, d in enumerate(dofs))


def build_analysis(sec: SectionSpec):
    from pre_processing.section_integration.thin_walled import (
        ThinWalledSection,
        ThinWalledSectionAnalysis,
        Wall,
        WallLaminate,
    )
    from pre_processing.section_integration.thin_walled.laminate import Ply, PlyMaterial

    mats = {i: PlyMaterial(m.E1, m.E2, m.G12, m.nu12, m.G13, m.G23, m.rho) for i, m in sec.materials.items()}
    lams = {lid: WallLaminate([Ply(mats[mid], a, t) for mid, a, t in plies]) for lid, plies in sec.laminates.items()}
    walls = [Wall(wid, a, b, lid) for wid, a, b, lid in sec.walls]
    return ThinWalledSectionAnalysis(ThinWalledSection(dict(sec.points), walls), lams, n_sub=sec.n_sub)


def results_dir(base: str) -> str:
    """A fresh job results directory with the standard sub-directories."""
    if os.path.isdir(base):
        shutil.rmtree(base)
    for sub in ("element_stiffness_matrices", "element_force_vectors", "primary_results", "secondary_results",
                "logs", "maps", "diagnostics"):
        os.makedirs(os.path.join(base, sub))
    return base


@dataclass
class JobResult:
    U: np.ndarray            # (n_nodes, stride)
    wall_rows: List[dict]    # rows of wall_stresses.csv


def run(spec: CantileverSpec, work_dir: str, name: str) -> JobResult:
    """Write, run and read one job under ``work_dir/name``."""
    import logging

    from workflow_orchestrator.run_job import process_job

    job_dir = os.path.join(work_dir, "jobs", name)
    write_job(job_dir, spec)
    write_supports(job_dir, spec)
    res = results_dir(os.path.join(work_dir, "results", name))
    logging.getLogger().setLevel(logging.WARNING)
    process_job(job_dir, res, {}, {}, force_serial=True, max_processes_per_job=1)
    with open(os.path.join(res, "primary_results", "global", "U_global.csv"), newline="", encoding="utf-8") as f:
        U = np.array([float(r[1]) for r in list(csv.reader(f))[1:]])
    with open(os.path.join(res, "secondary_results", "wall_stresses", "wall_stresses.csv"), newline="",
              encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    return JobResult(U.reshape(spec.n_el + 1, -1), rows)


def contour(rows: List[dict], element_id: int, face: str = "mid") -> List[dict]:
    """Rows of one element along the contour, one per wall point (ply 0 bottom face)."""
    out = [r for r in rows if int(r["element_id"]) == element_id and int(r["ply"]) == 0]
    seen, uniq = set(), []
    for r in out:
        key = (r["wall_id"], r["s"])
        if key not in seen:
            seen.add(key)
            uniq.append(r)
    return uniq


def bredt_multicell(analysis, torque: float) -> Dict[int, float]:
    """Cell shear flows under a torque by the classical multi-cell Bredt-Batho method.

    Each wall is a membrane with shear stiffness ``1 / inv(A)[2, 2]`` (no plate
    twisting). Returns the shear flow in each wall, positive along the wall.
    """
    sec = analysis.section
    walls = list(sec.walls)
    pts = sec.points
    # cells from the wall graph: fundamental cycles of a spanning tree

    adj: Dict[int, List[Tuple[int, int, int]]] = {}
    for w in walls:
        adj.setdefault(w.start, []).append((w.end, w.wall_id, +1))
        adj.setdefault(w.end, []).append((w.start, w.wall_id, -1))
    root = walls[0].start
    parent = {root: None}
    order = [root]
    for p in order:
        for q, wid, sgn in adj[p]:
            if q not in parent:
                parent[q] = (p, wid, sgn)
                order.append(q)
    tree = {v[1] for v in parent.values() if v is not None}

    def path_to_root(p):
        out = []
        while parent[p] is not None:
            q, wid, sgn = parent[p]
            out.append((wid, -sgn))  # from p toward q
            p = q
        return out

    cycles = []
    for w in walls:
        if w.wall_id in tree:
            continue
        # cycle: start -> end along w, end -> root, root -> start
        cyc = {w.wall_id: 1.0}
        for wid, s in path_to_root(w.end):
            cyc[wid] = cyc.get(wid, 0.0) + s
        for wid, s in path_to_root(w.start):
            cyc[wid] = cyc.get(wid, 0.0) - s
        cycles.append({k: v for k, v in cyc.items() if v != 0.0})
    wall_of = {w.wall_id: w for w in walls}

    def area(cyc):
        a = 0.0
        for wid, s in cyc.items():
            w = wall_of[wid]
            (y0, z0), (y1, z1) = (pts[w.start], pts[w.end]) if s > 0 else (pts[w.end], pts[w.start])
            a += 0.5 * (y0 * z1 - y1 * z0)
        return a

    areas = [area(c) for c in cycles]
    for i, a in enumerate(areas):  # orient counter-clockwise
        if a < 0:
            cycles[i] = {k: -v for k, v in cycles[i].items()}
            areas[i] = -a
    flex = {w.wall_id: sec.length(w) * np.linalg.inv(analysis.laminates[w.laminate_id].A())[2, 2] for w in walls}
    nc = len(cycles)
    # unknowns q_1..q_n and twist rate; rows: compatibility per cell, torque
    M = np.zeros((nc + 1, nc + 1))
    rhs = np.zeros(nc + 1)
    for i, ci in enumerate(cycles):
        for j, cj in enumerate(cycles):
            M[i, j] = sum(ci[w] * cj[w] * flex[w] for w in ci if w in cj)
        M[i, nc] = -2.0 * areas[i]
    M[nc, :nc] = 2.0 * np.array(areas)
    rhs[nc] = torque
    sol = np.linalg.solve(M, rhs)
    q = {w.wall_id: 0.0 for w in walls}
    for i, c in enumerate(cycles):
        for wid, s in c.items():
            q[wid] += s * sol[i]
    return q

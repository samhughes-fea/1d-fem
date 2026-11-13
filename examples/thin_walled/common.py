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
    offsets: Dict[int, float] = field(default_factory=dict)  # laminate_id -> mid-plane above the wall line


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
      + "[Laminate]\n[laminate_id] [ply] [material_id] [angle_deg] [thickness]\n" + lams
      + ("[Laminate_Offset]\n[laminate_id] [offset]\n"
         + "".join(f"{lid} {e:.9g}\n" for lid, e in sec.offsets.items()) if sec.offsets else ""))
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
    from onedfem.pre_processing.element_library.element_factory import ElementFactory

    fixed = int(ElementFactory.capability(spec.element_type, "extra_dofs", 0) or 0)
    if fixed:
        return fixed
    if not ElementFactory.capability(spec.element_type, "cross_section_modes", False):
        return 1 if spec.warping else 0
    from anisotropic_section_model.thin_walled.cross_section_modes import CrossSectionModes

    return CrossSectionModes(build_analysis(spec.section)).n_sv_warping + spec.shear_modes + 2 * spec.modes


def write_supports(job_dir: str, spec: CantileverSpec) -> None:
    dofs = STANDARD_DOFS + [f"Q{k}" for k in range(1, n_extra_dofs(spec) + 1) if k not in spec.free_root_extras]
    with open(os.path.join(job_dir, "prescribed_displacement.txt"), "w", encoding="utf-8") as f:
        f.write("[Prescribed Displacement]\n[id] [node_id] [dof] [value] [type] [comment]\n")
        f.writelines(f"{i} 0 {d} 0.0 displacement\n" for i, d in enumerate(dofs))


def build_analysis(sec: SectionSpec):
    from anisotropic_section_model.thin_walled import (
        ThinWalledSection,
        ThinWalledSectionAnalysis,
        Wall,
        WallLaminate,
    )
    from anisotropic_section_model.thin_walled.laminate import Ply, PlyMaterial

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

    from onedfem.workflow_orchestrator.run_job import process_job

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


def by_ply(points) -> List[list]:
    """Stress recovery points grouped by ply, in order through the thickness, for drawing
    each ply on its own (the recovered stresses are piecewise in ``n``, broken at ply faces)."""
    out: Dict[int, list] = {}
    for p in points:
        out.setdefault(p.ply, []).append(p)
    return [out[k] for k in sorted(out)]


def ply_materials(sec: SectionSpec, names: Dict[int, str]):
    """``[(label, PlyMaterial)]`` of the section's materials, for ``draw_section``."""
    from anisotropic_section_model.thin_walled.laminate import PlyMaterial

    return [(names[i], PlyMaterial(m.E1, m.E2, m.G12, m.nu12, m.G13, m.G23, m.rho)) for i, m in sec.materials.items()]


def _centres(ax, an):
    from onedfem.post_processing.graphical_visualisers import figure_style as fs

    ax.plot(*an.shear_centre, "+", color=fs.INK, ms=9, mew=1.2, label="shear centre", zorder=6)
    ax.plot(*an.tension_centre, "o", color=fs.INK, mfc="white", ms=5, label="tension centre", zorder=6)
    ax.plot(0, 0, "x", color=fs.INK, ms=6, label="reference axis", zorder=6)


def plot_section(an, mats, path: str, title: str, note: str, zoom, zoom_title: str, figsize=(10, 4.4),
                 width_ratios=(2.6, 1)) -> None:
    """The section at true thickness, ply by ply, and a close-up ``zoom = (y0, z0, dy, dz)``."""
    from onedfem.post_processing.graphical_visualisers import figure_style as fs
    from onedfem.post_processing.graphical_visualisers.thin_walled_drawing import draw_section

    fs.use()
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle

    fig = plt.figure(figsize=figsize)
    gs = fig.add_gridspec(1, 2, width_ratios=list(width_ratios))
    ax = fig.add_subplot(gs[0])
    handles = draw_section(ax, an, mats, legend=False)
    _centres(ax, an)
    ax.add_patch(Rectangle(zoom[:2], zoom[2], zoom[3], fill=False, ec=fs.MUTED, lw=0.7, ls=":"))
    ax.set(xlabel="y (m)", ylabel="z (m)")
    ax.legend(handles=handles + ax.get_legend_handles_labels()[0], loc="upper center", bbox_to_anchor=(0.5, -0.2),
              ncol=6, fontsize=7)
    fs.panel_title(ax, title, note)
    az = fig.add_subplot(gs[1])
    draw_section(az, an, mats, legend=False, nodes=False, ply_lines=True)
    az.set(xlim=(zoom[0], zoom[0] + zoom[2]), ylim=(zoom[1], zoom[1] + zoom[3]), xlabel="y (m)")
    az.grid(False)
    fs.panel_title(az, zoom_title, "every ply drawn; lighter = off-axis")
    fig.tight_layout()
    fs.save(fig, path)


def plot_section_fields(an, rec, mats, fields, path: str, size: float, figsize, ylim=None, labels=None,
                        xlim=None, across: bool = False) -> None:
    """Wall fields drawn normal to the walls, one panel per ``(title, forces, component, unit, unit name)``.

    ``size`` is the drawn length of the largest value. ``labels`` maps a panel index to
    ``[(wall_id, s, (dy, dz), text)]`` annotations.
    """
    from onedfem.post_processing.graphical_visualisers import figure_style as fs
    from onedfem.post_processing.graphical_visualisers.thin_walled_drawing import (
        contour_field,
        draw_section,
        gauss_s,
    )

    fs.use()
    import matplotlib.pyplot as plt

    shape = (1, len(fields)) if across else (len(fields), 1)
    fig, axes = plt.subplots(*shape, figsize=figsize, sharex=True, sharey=True)
    axes = np.atleast_1d(axes)
    for i, (ax, (title, F, comp, unit, uname)) in enumerate(zip(axes, fields)):
        draw_section(ax, an, mats, legend=False, nodes=False, angle_shading=False, edge=False)
        vmax = max(abs(rec.at(w.wall_id, float(s), F).resultants[comp]) for w in an.section.walls
                   for s in gauss_s(an, w.wall_id))
        contour_field(ax, an, lambda w, s, F=F, c=comp: rec.at(w, s, F).resultants[c], size / vmax)
        _centres(ax, an)
        ax.set(ylabel="z (m)" if (i == 0 or not across) else "")
        fs.panel_title(ax, title, f"drawn normal to the walls, {size * 1e3:.0f} mm = {vmax / unit:.1f} {uname}; "
                                  + ("\n" if across else "") + "orange +, blue −; dots: Gauss points")
        for w, s, (dy, dz), text in (labels or {}).get(i, []):
            p = an.point_at(w, s).position
            ax.annotate(text, p + np.array([dy, dz]), fontsize=6.5, ha="center", va="center", color=fs.INK,
                        bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.85), zorder=8)
    if ylim:
        axes[0].set_ylim(*ylim)
    if xlim:
        axes[0].set_xlim(*xlim)
    for ax in (axes if across else axes[-1:]):
        ax.set_xlabel("y (m)")
    fig.tight_layout()
    fs.save(fig, path)


def plot_deflection(results, element_types, length: float, n_el: int, exact, path: str,
                    twist_note: str = "") -> None:
    """Flapwise deflection and twist along the beam: each element's shape-function field, the
    nodal values, and the exact Saint-Venant tip ``exact = (w, twist)``."""
    from onedfem.post_processing.graphical_visualisers import figure_style as fs
    from onedfem.post_processing.graphical_visualisers.beam_field_sampling import beam_curve

    fs.use()
    import matplotlib.pyplot as plt

    x = np.linspace(0, length, n_el + 1)
    fig, ax = plt.subplots(1, 2, figsize=(10, 3.8))
    for i, (k, r) in enumerate(results.items()):
        c = fs.colour(k, i)
        for a, comp, scale in ((ax[0], 2, 1e3), (ax[1], 3, np.rad2deg(1.0))):
            xs, v = beam_curve(x, r.U, comp, element_types[k])
            fs.field(a, xs, v * scale, color=c, label=k, ls="--" if "free" in k else "-")
            fs.nodal(a, x, r.U[:, comp] * scale, color=c)
    fs.reference(ax[0], [length], [exact[0] * 1e3], label="exact Saint-Venant tip")
    fs.reference(ax[1], [length], [np.rad2deg(exact[1])], label="exact Saint-Venant tip")
    ax[0].set(xlabel="x (m)", ylabel="flapwise deflection w (mm)")
    ax[1].set(xlabel="x (m)", ylabel="twist (deg)")
    fs.panel_title(ax[0], "Flapwise deflection", f"{n_el} two-node elements")
    fs.panel_title(ax[1], "Twist", twist_note or None)
    ax[0].legend(loc="upper left")
    fs.marks_legend(ax[1], ("field", "nodal", "reference"), loc="upper left")
    fig.tight_layout()
    fs.save(fig, path)


def plot_interlaminar(ax, points, title: str = "Mid-span, through the thickness") -> None:
    """``tau_xn`` and ``tau_sn`` through each wall of ``points = [(WallPointStress, label)]``, ply by ply."""
    from onedfem.post_processing.graphical_visualisers import figure_style as fs

    for i, (q, lab) in enumerate(points):
        for j, ply in enumerate(by_ply(q.plies)):
            n = [p.n * 1e3 for p in ply]
            fs.field(ax, [p.tau_xn / 1e3 for p in ply], n, color=fs.colour("", 2 * i),
                     label=f"τ_xn, {lab}" if j == 0 else None, lw=1.1)
            fs.field(ax, [p.tau_sn / 1e3 for p in ply], n, color=fs.colour("", 2 * i + 1), ls="--",
                     label=f"τ_sn, {lab}" if j == 0 else None, lw=1.1)
    ax.axvline(0.0, color=fs.MUTED, lw=0.6)
    ax.set(xlabel="interlaminar shear (kPa)", ylabel="n (mm)")
    fs.panel_title(ax, title, "each ply drawn on its own; zero at both faces")
    ax.legend(loc="lower right")


def plot_ply_stresses(laminate, point, path: str, title: str, names: Dict[int, str]) -> None:
    """Ply-axis stresses through one wall, each ply on its own, ply bands shaded by material.

    ``names`` maps a material's ``E1`` in GPa (rounded) to its name.
    """
    from onedfem.post_processing.graphical_visualisers import figure_style as fs

    fs.use()
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(1, 3, figsize=(11, 3.9), sharey=True)
    z = laminate.interfaces * 1e3
    for i, lab in enumerate(("σ_11", "σ_22", "τ_12")):
        for k, ply in enumerate(laminate.plies):
            name = names.get(round(ply.material.E1 / 1e9), "")
            c = fs.MATERIAL_COLOURS.get(name, fs.MUTED)
            ax[i].axhspan(z[k], z[k + 1], color=c, alpha=0.07 + 0.12 * abs(ply.angle_deg) / 90.0, lw=0)
            ax[i].axhline(z[k], color="white", lw=0.4)
        if i == 0:   # one label per run of equal plies
            k = 0
            plies = laminate.plies
            while k < len(plies):
                j = k
                while j + 1 < len(plies) and plies[j + 1].angle_deg == plies[k].angle_deg \
                        and plies[j + 1].material == plies[k].material:
                    j += 1
                name = names.get(round(plies[k].material.E1 / 1e9), "")
                a = plies[k].angle_deg
                text = f"{name} {a:+.0f}°" if a else f"{name} 0°"
                ax[i].text(0.99, 0.5 * (z[k] + z[j + 1]), text + (f" ×{j - k + 1}" if j > k else ""),
                           transform=ax[i].get_yaxis_transform(), fontsize=6, color=fs.MUTED, ha="right",
                           va="center")
                k = j + 1
        for ply in by_ply(point.plies):
            fs.field(ax[i], [p.stress_ply[i] / 1e6 for p in ply], [p.n * 1e3 for p in ply], color=fs.INK)
        ax[i].set(xlabel=f"{lab} (MPa)")
        ax[i].grid(axis="y", visible=False)
    ax[0].set_ylabel("n through the wall (mm)")
    fs.panel_title(ax[0], "Fibre direction", "bands: plies, darker = off-axis")
    fs.panel_title(ax[1], "Transverse")
    fs.panel_title(ax[2], "In-plane shear")
    fig.suptitle(title, x=0.01, ha="left")
    fig.tight_layout()
    fs.save(fig, path)

# examples/thin_walled/gbt_buckling/run.py
"""Local, distortional and global buckling of a lipped channel column with GBT elements.

Run from the repository root::

    python examples/thin_walled/gbt_buckling/run.py

An aluminium lipped channel (web 100 mm, flanges 60 mm, lips 10 mm, walls 2 mm) is
compressed at its centroid between simply supported ends (``v = w = theta_x = 0`` and a
diaphragm, ``a_k = 0``, at both ends; warping free). 1d-fem's linear buckling with
``LinearGBTBeamElement3D`` is run for lengths from 40 mm to 3 m. The section carries its
distortional and local plate modes (``n_sub = 2``), so one model buckles locally in short
columns, distortionally at intermediate lengths and globally (flexural-torsional) in long
ones.

The reference is a finite strip model of the same column (``finite_strip.py``: full plate
theory, no GBT hypotheses, lowest over the number of half-waves). Benscoter beam elements,
which keep the section shape, give the global buckling alone.

Writes ``jobs/`` and ``results/`` (not tracked), and ``summary.json``, ``signature.png``
and ``mode_shapes.png`` next to this file.
"""

from __future__ import annotations

import glob
import json
import logging
import os
import sys
from typing import Dict, List

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
for p in (REPO, os.path.dirname(HERE), HERE):
    if p not in sys.path:
        sys.path.insert(0, p)

from common import CantileverSpec, Material, SectionSpec, build_analysis, results_dir, write_job  # noqa: E402
from finite_strip import column_stress  # noqa: E402

E, NU = 70e9, 0.3
ALU = Material(E1=E, E2=E, G12=E / (2 * (1 + NU)), nu12=NU, G13=E / (2 * (1 + NU)), G23=E / (2 * (1 + NU)), rho=2700.0)
H, B, C, T = 0.100, 0.060, 0.010, 2e-3        # web, flange, lip, wall thickness (centre lines), m
N_SUB, MODES = 2, 10
LENGTHS = (0.04, 0.06, 0.08, 0.12, 0.2, 0.3, 0.4, 0.6, 0.9, 1.4, 2.0, 3.0)
SHOWN = (0.08, 0.4, 2.0)                       # lengths whose buckled sections are drawn
EL_LENGTH = 0.01                               # target element length, m (local half-waves are ~80 mm)
CLASSES = ("global", "distortional", "local")
COLOURS = {"global": "#2a78d6", "distortional": "#eb6834", "local": "#1baf7a", "fsm": "#52514e"}


def points() -> Dict[int, tuple]:
    """Centre-line points with the centroid at the origin (the beam axis)."""
    raw = {1: (B, H / 2 - C), 2: (B, H / 2), 3: (0.0, H / 2), 4: (0.0, -H / 2), 5: (B, -H / 2), 6: (B, -H / 2 + C)}
    y_c = (2 * B * (B / 2) + 2 * C * B) / (H + 2 * B + 2 * C)
    return {k: (y - y_c, z) for k, (y, z) in raw.items()}


WALLS = [(1, 1, 2), (2, 2, 3), (3, 3, 4), (4, 4, 5), (5, 5, 6)]


def section(n_sub: int = N_SUB) -> SectionSpec:
    return SectionSpec(points=points(), walls=[(w, a, b, 1) for w, a, b in WALLS], materials={1: ALU},
                       laminates={1: [(1, 0.0, T)]}, n_sub=n_sub)


def n_elements(length: float) -> int:
    n = int(np.clip(np.ceil(length / EL_LENGTH), 16, 120))
    return n + n % 2                            # even: a node at midspan


def write_column_job(job: str, length: float, element_type: str, modes: int, n_warping: int) -> int:
    """A pinned column, compressed by a unit force at the centroid. Returns the element count."""
    n_el = n_elements(length)
    write_job(job, CantileverSpec(length, n_el, element_type, section(), [0.0] * 6, modes=modes,
                                  warping=1 if "Benscoter" in element_type else 0))
    with open(os.path.join(job, "point_load.txt"), "w", encoding="utf-8") as f:
        f.write("[Point load]\n[x] [y] [z] [F_x] [F_y] [F_z] [M_x] [M_y] [M_z]\n")
        f.write(f"{length:.9g} 0 0 -1 0 0 0 0 0\n")
    diaphragm = [f"Q{n_warping + 1 + 2 * k}" for k in range(modes)]          # a_k = 0
    rows = [(0, d) for d in ["UX", "UY", "UZ", "RX"] + diaphragm] + [(n_el, d) for d in ["UY", "UZ", "RX"] + diaphragm]
    with open(os.path.join(job, "prescribed_displacement.txt"), "w", encoding="utf-8") as f:
        f.write("[Prescribed Displacement]\n[id] [node_id] [dof] [value] [type] [comment]\n")
        f.writelines(f"{i} {n} {d} 0.0 displacement\n" for i, (n, d) in enumerate(rows))
    with open(os.path.join(job, "simulation_settings.txt"), "w", encoding="utf-8") as f:
        f.write("[Simulation]\n[Type]\nbuckling\n\n[Solver]\ntype = direct\ntolerance = 1e-12\n\n"
                "[Buckling]\nenabled = true\nnum_modes = 4\nbuckling_prestress = linear_static\n"
                "buckling_load_factor = 1.0\n")
    return n_el


def run_buckling(work: str, name: str, length: float, element_type: str, modes: int, n_warping: int):
    """Load factors (ascending) and mode shapes (n_nodes, stride, n_modes) of one job."""
    from onedfem.workflow_orchestrator.run_job import process_job

    job = os.path.join(work, "jobs", name)
    n_el = write_column_job(job, length, element_type, modes, n_warping)
    res = results_dir(os.path.join(work, "results", name))
    process_job(job, res, {}, {}, force_serial=True, max_processes_per_job=1)
    f = glob.glob(os.path.join(res, "primary_results", "modal_results", "*_buckling_load_factors.txt"))[0]
    lam = np.atleast_1d(np.loadtxt(f))
    V = np.loadtxt(f.replace("load_factors", "mode_shapes")).reshape(-1, lam.size)
    order = np.argsort(lam)
    return lam[order], V[:, order].reshape(n_el + 1, -1, lam.size)


# --------------------------------------------------------------------------- #
# Buckled section shapes, from the section modes and the frame shape functions
# --------------------------------------------------------------------------- #
class SectionShape:
    """In-plane displacement of the section for nodal DOFs ``[v, w, theta_x, a(m)]``.

    Each mode ``a_k`` displaces the section frame by its nodal vector; the beam DOFs move it
    rigidly. Along every frame element the displacement is that of the frame's own shape
    functions (tangential linear, normal cubic Hermite), so the drawn shape is the exact
    field of the GBT section model.
    """

    def __init__(self, modes: int):
        from anisotropic_section_model.thin_walled.cross_section_modes import CrossSectionModes

        self.cm = CrossSectionModes(build_analysis(section()))
        self.X = self.cm._modes[1][:, :modes]
        ends = [pid for _, a, b in WALLS for pid in (a, b)]
        folds = [points()[pid] for pid in set(ends) if ends.count(pid) > 1]
        corners = [i for i, xy in enumerate(self.cm.node_xy) if any(np.allclose(xy, p) for p in folds)]
        # a local mode leaves the folds (where walls meet) in place; a distortional one moves them
        self.kind = []
        for k in range(modes):
            m = self.cm.mode(k)
            d = np.hypot(m["dy"], m["dz"])
            self.kind.append("local" if d[corners].max() < 0.05 * d.max() else "distortional")

    def frame_vector(self, q: np.ndarray, which: str = "all") -> np.ndarray:
        """Frame nodal vector ``[dy, dz, beta, u]`` per node for ``q = [v, w, theta, a]``."""
        y, z = self.cm.node_xy[:, 0], self.cm.node_xy[:, 1]
        out = np.zeros(self.cm.n_frame_dof)
        if which in ("all", "global"):
            out[0::4] += q[0] - z * q[2]
            out[1::4] += q[1] + y * q[2]
            out[2::4] += q[2]
        for k in range(self.X.shape[1]):
            if which in ("all", self.kind[k]):
                out += q[3 + k] * self.X[:, k]
        return out

    def curve(self, f: np.ndarray, n: int = 9) -> List[np.ndarray]:
        """Displaced points (n per frame element) of every frame element for frame vector ``f``."""
        out = []
        for e in self.cm._elements:
            d = self.cm._element_dofs(e)
            pts = []
            for xi in np.linspace(0.0, 1.0, n):
                fl = self.cm._fields(e, xi)
                pts.append(e["a"] + xi * e["le"] * e["t"] + (fl["v"] @ f[d]) * e["t"] + (fl["w"] @ f[d]) * e["n"])
            out.append(np.array(pts))
        return out

    def share(self, Q: np.ndarray) -> Dict[str, float]:
        """Share of each class in the in-plane displacement, summed over the column nodes.

        ``Q`` (n_nodes, 3 + m). The measure is the frame nodal displacement squared of each
        class, as a fraction of the sum over the classes.
        """
        tot = {c: 0.0 for c in CLASSES}
        for q in Q:
            for c in CLASSES:
                f = self.frame_vector(q, c)
                tot[c] += float(f[0::4] @ f[0::4] + f[1::4] @ f[1::4])
        s = sum(tot.values()) or 1.0
        return {c: v / s for c, v in tot.items()}


def in_plane_dofs(U: np.ndarray, n_warping: int, modes: int) -> np.ndarray:
    """``[v, w, theta_x, a(m)]`` per node from nodal DOF rows."""
    a = U[:, 6 + n_warping:6 + n_warping + 2 * modes:2]
    return np.column_stack([U[:, 1], U[:, 2], U[:, 3], a])


# --------------------------------------------------------------------------- #
def main(plot: bool = True, work_dir: str = HERE, out_dir: str = HERE, lengths=LENGTHS) -> dict:
    from anisotropic_section_model.thin_walled.cross_section_modes import CrossSectionModes

    area = T * (H + 2 * B + 2 * C)
    p = CrossSectionModes(build_analysis(section())).n_sv_warping
    shape = SectionShape(MODES)
    walls = [(points()[a], points()[b]) for _, a, b in WALLS]
    rows, shapes = [], {}
    logging.disable(logging.WARNING)
    try:
        for L in lengths:
            lam, V = run_buckling(work_dir, f"gbt_{L:g}", L, "LinearGBTBeamElement3D", MODES, p)
            lam_b, _ = run_buckling(work_dir, f"benscoter_{L:g}", L, "LinearBenscoterBeamElement3D", 0, 1)
            Q = in_plane_dofs(V[:, :, 0], p, MODES)
            share = shape.share(Q)
            s_fsm, waves = column_stress(walls, T, E, NU, L)
            rows.append({"length": L, "elements": V.shape[0] - 1, "gbt": float(lam[0] / area),
                         "fsm": s_fsm, "fsm_half_waves": waves, "benscoter": float(lam_b[0] / area),
                         "gbt_over_fsm": float(lam[0] / area / s_fsm), "share": share,
                         "governs": max(share, key=share.get)})
            if any(np.isclose(L, s) for s in SHOWN):
                shapes[L] = Q
    finally:
        logging.disable(logging.NOTSET)
    result = {"section": {"web": H, "flange": B, "lip": C, "t": T, "E": E, "nu": NU, "n_sub": N_SUB,
                          "modes": MODES, "mode_kinds": shape.kind},
              "rows": rows,
              "max_gbt_error": max(abs(r["gbt_over_fsm"] - 1) for r in rows)}
    with open(os.path.join(out_dir, "summary.json"), "w", encoding="utf-8") as f:
        json.dump(result, f, indent=1)
    if plot:
        figures(rows, shapes, shape, out_dir)
    return result


def figures(rows, shapes, shape: SectionShape, out_dir: str) -> None:
    """Buckling stress against length (one marker per solved length), the share of each mode
    class in the buckling mode, and the buckled sections drawn through the frame shape functions."""
    from onedfem.post_processing.graphical_visualisers import figure_style as fs

    fs.use()
    import matplotlib.pyplot as plt

    L = np.array([r["length"] for r in rows])
    fig, (ax, axs) = plt.subplots(2, 1, figsize=(7.5, 5.6), sharex=True, gridspec_kw={"height_ratios": [3, 1]})
    fs.reference(ax, L, [r["fsm"] / 1e6 for r in rows], label="finite strip (reference)", ms=7)
    for c in CLASSES:
        sel = [r for r in rows if r["governs"] == c]
        if sel:
            fs.solves(ax, [r["length"] for r in sel], [r["gbt"] / 1e6 for r in sel], color=COLOURS[c],
                      label=f"1d-fem GBT, {c}")
    ax.plot(L, [r["benscoter"] / 1e6 for r in rows], "x", ls="none", ms=6, color=fs.MUTED,
            label="1d-fem Benscoter (rigid section)")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_ylabel("buckling stress P / A (MPa)")
    fs.panel_title(ax, "Lipped channel column: lowest buckling stress at each length",
                   "one marker per column length solved; GBT within 1.3 % of the finite strip model")
    ax.legend(loc="upper right")
    bottom = np.zeros(len(rows))
    for c in ("local", "distortional", "global"):
        v = np.array([r["share"][c] for r in rows])
        axs.bar(L, v, bottom=bottom, width=L * 0.18, color=COLOURS[c], label=c, lw=0)
        bottom += v
    axs.set(xlabel="column length (m)", ylabel="mode share", ylim=(0, 1))
    from matplotlib.ticker import FuncFormatter

    for a in (ax, axs):
        a.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
    ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
    axs.grid(axis="x", visible=False)
    fs.panel_title(axs, "Share of the in-plane displacement by mode class", None)
    axs.legend(loc="center left", bbox_to_anchor=(1.0, 0.5))
    fig.tight_layout()
    fs.save(fig, os.path.join(out_dir, "signature.png"))

    fig, axes = plt.subplots(1, len(shapes), figsize=(3.3 * len(shapes), 3.6))
    by_len = {r["length"]: r for r in rows}
    for ax, (Lk, Q) in zip(np.atleast_1d(axes), sorted(shapes.items())):
        mid = Q.shape[0] // 2
        q = Q[mid] if np.abs(Q[mid]).max() > 0.2 * np.abs(Q).max() else Q[np.argmax(np.abs(Q).max(axis=1))]
        f = shape.frame_vector(q)
        peak = max(np.hypot(f[0::4], f[1::4]).max(), 1e-30)
        f = f * (0.25 * H / peak)
        for seg in shape.curve(np.zeros_like(f)):
            ax.plot(seg[:, 0] * 1e3, seg[:, 1] * 1e3, color=fs.GRID, lw=3.0)
        r = by_len[Lk]
        for seg in shape.curve(f):
            fs.field(ax, seg[:, 0] * 1e3, seg[:, 1] * 1e3, color=COLOURS[r["governs"]], lw=1.5)
        ax.set_aspect("equal")
        ax.set_xlabel("y (mm)")
        fs.panel_title(ax, f"L = {Lk:g} m: {r['governs']}",
                       f"{r['gbt'] / 1e6:.1f} MPa (finite strip {r['fsm'] / 1e6:.1f})")
    np.atleast_1d(axes)[0].set_ylabel("z (mm)")
    fig.suptitle("Buckled section at a node near mid-span, through the frame shape functions (grey: undeformed)",
                 x=0.01, ha="left")
    fig.tight_layout()
    fs.save(fig, os.path.join(out_dir, "mode_shapes.png"))


if __name__ == "__main__":
    out = main()
    for r in out["rows"]:
        print(f"L = {r['length']:5.2f} m  GBT {r['gbt'] / 1e6:7.2f} MPa  FSM {r['fsm'] / 1e6:7.2f} MPa "
              f"({r['fsm_half_waves']} half-waves)  ratio {r['gbt_over_fsm']:.4f}  Benscoter {r['benscoter'] / 1e6:8.2f}  "
              f"{r['governs']}")

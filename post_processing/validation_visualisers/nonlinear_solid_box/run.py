"""Corotational Timoshenko beam against the CalculiX ``NLGEOM`` solid box.

Runs the 1d-fem jobs of the CAS composite box (``make_reference.py``) at every load level of
the solid reference, with two sections:

- ``thin_walled``: the thin-walled section of the shell validation (walls on their mid-lines);
- ``2d``: the 2D section mesh that the solid model is extruded from, so the beam and the
  solid have the same section stiffness and the comparison isolates the large-rotation
  kinematics.

Each section also runs once with the linear Timoshenko element at load ``P``. Writes
``summary.json`` and ``nonlinear_box.png`` here. The beams take about a minute::

    python post_processing/validation_visualisers/nonlinear_solid_box/run.py
"""

from __future__ import annotations

import json
import os
import sys
import tempfile
from typing import Dict

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
if REPO not in sys.path:
    sys.path.insert(0, REPO)

N_EL = 24
INCREMENTS_PER_P = 8


def reference() -> dict:
    with open(os.path.join(HERE, "reference.json"), encoding="utf-8") as f:
        return json.load(f)


def _beam(level: float, element: str) -> Dict:
    from post_processing.validation_visualisers.abaqus.shell_validation.cases import BOX_L, box_tip_force

    P = box_tip_force()
    return dict(kind="nonlinear", element=element, length=BOX_L, n_el=N_EL, warping=0, supports="clamp",
                loads=[(BOX_L, [0, 0, P, 0, 0, 0])], levels=(level,),
                increments=max(1, int(round(INCREMENTS_PER_P * level))))


def run_job(work: str, section: str, level: float, linear: bool = False) -> np.ndarray:
    """Nodal displacements (rows: beam nodes; columns: UX, UY, UZ, RX, RY, RZ) of one job."""
    import csv

    from post_processing.validation_visualisers.abaqus.shell_validation import onedfem
    from post_processing.validation_visualisers.abaqus.shell_validation.cases import box_section
    from post_processing.validation_visualisers.nonlinear_solid_box.make_reference import box_mesh
    from pre_processing.section_integration.solid.from_geometry import write_section_mesh_2d

    element = "LinearTimoshenkoBeamElement3D" if linear else "CorotationalTimoshenkoBeamElement3D"
    beam = _beam(level, element)
    name = f"{section}_{'linear' if linear else f'{level:g}'}"
    job = os.path.join(work, "jobs", name)
    res = os.path.join(work, "results", name)
    P = beam["loads"][0][1][2]
    settings = (onedfem.LINEAR_SETTINGS if linear else
                onedfem.NONLINEAR_SETTINGS.format(tol=1e-8 * max(P * level, 1.0), inc=beam["increments"]))
    onedfem.write_beam_job(job, box_section(), beam, settings, scale=level)
    if section == "2d":
        os.remove(os.path.join(job, "section_geometry.txt"))
        os.remove(os.path.join(job, "laminate.txt"))
        write_section_mesh_2d(os.path.join(job, "section_mesh_2d.txt"), {1: box_mesh()})
    onedfem._process(job, res)
    if not linear:
        with open(os.path.join(res, "logs", "newton_history.csv"), newline="") as fh:
            last = list(csv.DictReader(fh))[-1]
        if float(last["load_factor"]) < 1.0 - 1e-12 or last["residual_ok"] != "1":
            raise RuntimeError(f"{name}: Newton did not reach the full load (last row {last})")
    return onedfem._read_U(res, N_EL)[:, :6]


def compare(beams: Dict[str, Dict], ref: dict) -> dict:
    """Tip errors against the solid, as fractions of the solid's tip deflection ``|u|`` at
    each level (and of the twist for the rotation about x)."""
    L = ref["length"]
    out = {}
    nl = ref["nlgeom"]
    for sec, b in beams.items():
        rows = []
        for k, lv in enumerate(nl["levels"]):
            s = np.array(nl["tip"][k])
            t = np.array(b["tip"][k])
            scale = np.linalg.norm(s[:3])
            rows.append({"level": lv, "tip_deflection_over_span": float(s[2] / L),
                         "u_z": float((t[2] - s[2]) / s[2]), "u_x": float((t[0] - s[0]) / scale),
                         "theta_x": float((t[3] - s[3]) / s[3])})
        lin_s = np.array(ref["linear"]["tip"][0])
        lin_b = np.array(b["linear_tip"])
        out[sec] = {"levels": rows,
                    "linear": {"u_z": float((lin_b[2] - lin_s[2]) / lin_s[2]),
                               "theta_x": float((lin_b[3] - lin_s[3]) / lin_s[3])},
                    "max_abs_u_z": max(abs(r["u_z"]) for r in rows),
                    "max_abs_u_x": max(abs(r["u_x"]) for r in rows),
                    "max_abs_theta_x": max(abs(r["theta_x"]) for r in rows)}
    return out


def _marker(sec: str, colour: str) -> dict:
    # thin-walled: open squares; 2D: small filled triangles drawn on top
    if sec == "2d":
        return dict(marker="^", ls="none", ms=4, color=colour)
    return dict(marker="s", ls="none", ms=8, mfc="none", mew=1.2, color=colour)


def figures(beams: Dict[str, Dict], ref: dict, out_dir: str) -> None:
    """Load paths at the tip and the deformed axis at the top load.

    Everything is discrete: solid load levels and beam load steps are markers, the dashed
    lines are the linear solutions (exactly linear in the load), and the deformed axis is
    shown at the 1d-fem nodes and the solid's printed stations."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    L = ref["length"]
    nl = ref["nlgeom"]
    lv = np.array(nl["levels"])
    s_tip = np.array(nl["tip"])
    colours = {"2d": "#2a78d6", "thin_walled": "#d6762a"}
    labels = {"2d": "beam, 2D section", "thin_walled": "beam, thin-walled section"}
    fig, axes = plt.subplots(1, 4, figsize=(17, 4.2))
    panels = [(2, "flapwise tip deflection  $u_z / L$", 1.0 / L), (0, "tip shortening  $-u_x / L$", -1.0 / L),
              (3, "tip twist  $\\theta_x$ (rad)", 1.0)]
    for ax, (j, title, f) in zip(axes, panels):
        lin = ref["linear"]["tip"][0][j] * f
        ax.plot([0, lv[-1]], [0, lv[-1] * lin], "--", color="#c3c2b7", lw=1.0, label="solid, linear")
        ax.plot(lv, s_tip[:, j] * f, "o", ms=7, mfc="none", color="#222222", label="solid C3D20R, NLGEOM")
        for sec in ("thin_walled", "2d"):
            if sec in beams:
                ax.plot(lv, np.array(beams[sec]["tip"])[:, j] * f, **_marker(sec, colours[sec]), label=labels[sec])
        ax.set_title(title, fontsize=10)
        ax.set_xlabel("load level  $F / P$")
        ax.grid(alpha=0.3)
    axes[0].legend(fontsize=8, loc="upper left")
    ax = axes[3]
    depth = 0.1
    for k, m in zip(sorted(nl["flattening"]), ("o", "D", "v")):
        ax.plot(lv, 100 * np.array(nl["flattening"][k]) / depth, m, ms=6, mfc="none", color="#222222",
                label=f"solid, x = {float(k):g} L")
    ax.set_title("flange flattening (Brazier), % of depth", fontsize=10)
    ax.set_xlabel("load level  $F / P$")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8, loc="upper left")
    fig.suptitle("CAS composite box, 6 m cantilever, dead flapwise tip force "
                 f"(P = {ref['P']:.0f} N: 20 % of span in linear theory)", fontsize=10)
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, "nonlinear_box.png"), dpi=150)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6.5, 4.2))
    x = np.linspace(0, L, N_EL + 1)
    top = len(lv) - 1
    for sec in ("thin_walled", "2d"):
        if sec in beams:
            U = np.array(beams[sec]["top_nodes"])
            ax.plot(x + U[:, 0], U[:, 2], **_marker(sec, colours[sec]), label=f"{labels[sec]} (nodes)")
    st = nl["stations"]
    pts = [(0.0, 0.0, 0.0)] + [(float(k) * L + st[k][top][0], st[k][top][2], 0) for k in st] \
        + [(L + s_tip[top, 0], s_tip[top, 2], 0)]
    ax.plot([p[0] for p in pts], [p[1] for p in pts], "o", ms=7, mfc="none", color="#222222",
            label="solid (mean section motion)")
    ax.set_aspect("equal")
    ax.set_xlabel("x (m)")
    ax.set_ylabel("z (m)")
    ax.set_title(f"Deformed axis at F = {lv[-1]:g} P", fontsize=10)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, "deformed_axis.png"), dpi=150)
    plt.close(fig)


def main(plot: bool = True, work_dir: str = "", out_dir: str = HERE, sections=("2d", "thin_walled"),
         plot_only: bool = False) -> dict:
    """Run the beams and compare; ``plot_only`` redraws the figures from ``summary.json``."""
    ref = reference()
    if plot_only:
        with open(os.path.join(out_dir, "summary.json"), encoding="utf-8") as f:
            summary = json.load(f)
        figures(summary["beams"], ref, out_dir)
        return summary
    work = work_dir or tempfile.mkdtemp(prefix="nl_box_beam_")
    beams: Dict[str, Dict] = {}
    for sec in sections:
        tips, top = [], None
        for lv in ref["nlgeom"]["levels"]:
            U = run_job(work, sec, lv)
            tips.append(U[-1].tolist())
            top = U
        lin = run_job(work, sec, 1.0, linear=True)
        beams[sec] = {"tip": tips, "linear_tip": lin[-1].tolist(), "top_nodes": top.tolist()}
    summary = {"reference": {k: ref[k] for k in ("element", "n_along", "elements_along_wall", "length", "P")},
               "beam": {"element": "CorotationalTimoshenkoBeamElement3D", "n_el": N_EL},
               "comparison": compare(beams, ref), "beams": beams}
    with open(os.path.join(out_dir, "summary.json"), "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=1)
    if plot:
        figures(beams, ref, out_dir)
    return summary


if __name__ == "__main__":
    s = main(plot_only="--plot-only" in sys.argv)
    for sec, c in s["comparison"].items():
        print(sec, "linear", c["linear"])
        for r in c["levels"]:
            print(f"  {r['level']:5.2f}  w/L={r['tip_deflection_over_span']:.3f}  u_z {100 * r['u_z']:+.2f} %  "
                  f"u_x {100 * r['u_x']:+.2f} %  theta_x {100 * r['theta_x']:+.2f} %")

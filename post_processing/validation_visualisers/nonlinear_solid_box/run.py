"""Corotational Timoshenko beam against the CalculiX ``NLGEOM`` solid box.

Runs the 1d-fem jobs of the CAS composite box (``make_reference.py``) at every load level of
the solid reference, with two sections:

- ``thin_walled``: the thin-walled section of the shell validation (walls on their mid-lines);
- ``2d``: the 2D section mesh that the solid model is extruded from, so the beam and the
  solid have the same section stiffness and the comparison isolates the large-rotation
  kinematics.

and ``gbt``: ``CorotationalGBTBeamElement3D`` on the thin-walled section with ``GBT_MODES``
distortional modes, whose Brazier term lets the flanges flatten. Its flattening is read
from the nodal mode amplitudes at the solid's stations.

Each Timoshenko section also runs once with the linear Timoshenko element at load ``P``.
Writes ``summary.json`` and ``nonlinear_box.png`` here. The Timoshenko beams take about a
minute, the GBT beam about half an hour. ``--sections=gbt`` reruns one beam and keeps the
others from ``summary.json``::

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
GBT_MODES = 8
STATIONS = ("0.25", "0.5", "0.75")


def reference() -> dict:
    with open(os.path.join(HERE, "reference.json"), encoding="utf-8") as f:
        return json.load(f)


def _beam(level: float, element: str, warping: int = 0) -> Dict:
    from post_processing.validation_visualisers.abaqus.shell_validation.cases import BOX_L, box_tip_force

    P = box_tip_force()
    return dict(kind="nonlinear", element=element, length=BOX_L, n_el=N_EL, warping=warping, supports="clamp",
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


def gbt_modes():
    """The thin-walled section's GBT modes, as the 1d-fem job builds them (hoop restrained)."""
    from post_processing.validation_visualisers.abaqus.shell_validation import onedfem
    from post_processing.validation_visualisers.abaqus.shell_validation.cases import box_section
    from pre_processing.section_integration.thin_walled.cross_section_modes import CrossSectionModes
    from pre_processing.section_integration.thin_walled.section_analysis import ThinWalledSectionAnalysis

    tw = onedfem.section_analysis(box_section())
    return CrossSectionModes(ThinWalledSectionAnalysis(tw.section, tw.laminates, hoop="restrained", n_sub=tw.n_sub))


def flattening(modes, amplitudes: np.ndarray) -> float:
    """How much closer the flanges' mid-points come than their corners, from the mode
    amplitudes (the modes' in-plane motion at the wall nodes; the definition of the solid's
    ``make_reference.flattening``)."""
    from post_processing.validation_visualisers.abaqus.shell_validation.cases import BOX_H, BOX_W

    xy = modes.node_xy

    def node(y, z):
        return int(np.argmin(np.hypot(xy[:, 0] - y, xy[:, 1] - z)))

    hw, hd = 0.5 * BOX_W, 0.5 * BOX_H
    ct, cb, mt, mb = node(hw, hd), node(hw, -hd), node(0.0, hd), node(0.0, -hd)
    dz = sum(a * modes.mode(k)["dz"] for k, a in enumerate(amplitudes))
    return float((dz[ct] - dz[cb]) - (dz[mt] - dz[mb]))


def run_gbt_job(work: str, level: float, modes) -> np.ndarray:
    """Nodal DOFs of one corotational GBT job (rows: beam nodes; columns: the six beam DOFs,
    the Saint-Venant warping, then each mode's amplitude and its derivative)."""
    import csv

    from post_processing.validation_visualisers.abaqus.shell_validation import onedfem
    from post_processing.validation_visualisers.abaqus.shell_validation.cases import box_section

    beam = _beam(level, "CorotationalGBTBeamElement3D", warping=1)
    name = f"gbt_{level:g}"
    job = os.path.join(work, "jobs", name)
    res = os.path.join(work, "results", name)
    P = beam["loads"][0][1][2]
    onedfem.write_beam_job(job, box_section(), beam,
                           onedfem.NONLINEAR_SETTINGS.format(tol=1e-8 * max(P * level, 1.0), inc=beam["increments"]),
                           scale=level)
    with open(os.path.join(job, "element_section.txt"), "w", encoding="utf-8") as f:
        f.write("[Element_Section]\n[element_id] [section_id] [modes] [shear_modes]\n"
                + "".join(f"{e} 1 {GBT_MODES} 0\n" for e in range(N_EL)))
    n_field = modes.n_sv_warping + 2 * GBT_MODES
    fixed = [d for d in ("UX", "UY", "UZ", "RX", "RY", "RZ")] + [f"Q{q + 1}" for q in range(n_field)]
    with open(os.path.join(job, "prescribed_displacement.txt"), "w", encoding="utf-8") as f:
        f.write("[Prescribed Displacement]\n[id] [node_id] [dof] [value] [type] [comment]\n"
                + "".join(f"{i} 0 {d} 0.0 displacement\n" for i, d in enumerate(fixed)))
    onedfem._process(job, res)
    with open(os.path.join(res, "logs", "newton_history.csv"), newline="") as fh:
        last = list(csv.DictReader(fh))[-1]
    if float(last["load_factor"]) < 1.0 - 1e-12 or last["residual_ok"] != "1":
        raise RuntimeError(f"{name}: Newton did not reach the full load (last row {last})")
    with open(os.path.join(res, "primary_results", "global", "U_global.csv"), newline="", encoding="utf-8") as f:
        U = np.array([float(r[1]) for r in list(csv.reader(f))[1:]])
    return U.reshape(N_EL + 1, -1)


def run_gbt(work: str, levels) -> Dict:
    modes = gbt_modes()
    nw = modes.n_sv_warping
    tips, top, flat = [], None, {st: [] for st in STATIONS}
    for lv in levels:
        U = run_gbt_job(work, lv, modes)
        tips.append(U[-1, :6].tolist())
        top = U[:, :6]
        for st in STATIONS:
            a = U[int(round(float(st) * N_EL)), 6 + nw::2][:GBT_MODES]
            flat[st].append(flattening(modes, a))
    return {"tip": tips, "top_nodes": top.tolist(), "flattening": flat}


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
        out[sec] = {"levels": rows,
                    "max_abs_u_z": max(abs(r["u_z"]) for r in rows),
                    "max_abs_u_x": max(abs(r["u_x"]) for r in rows),
                    "max_abs_theta_x": max(abs(r["theta_x"]) for r in rows)}
        if "linear_tip" in b:
            lin_s = np.array(ref["linear"]["tip"][0])
            lin_b = np.array(b["linear_tip"])
            out[sec]["linear"] = {"u_z": float((lin_b[2] - lin_s[2]) / lin_s[2]),
                                  "theta_x": float((lin_b[3] - lin_s[3]) / lin_s[3])}
        if "flattening" in b:
            out[sec]["flattening_ratio"] = {st: [float(fb / fs) for fb, fs in zip(b["flattening"][st], nl["flattening"][st])]
                                            for st in b["flattening"]}
    return out


def _marker(sec: str, colour: str) -> dict:
    # thin-walled: open squares; 2D: small filled triangles drawn on top; GBT: crosses
    if sec == "2d":
        return dict(marker="^", ls="none", ms=4, color=colour)
    if sec == "gbt":
        return dict(marker="x", ls="none", ms=7, mew=1.4, color=colour)
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
    colours = {"2d": "#2a78d6", "thin_walled": "#d6762a", "gbt": "#1a9e5a"}
    labels = {"2d": "Timoshenko, 2D section", "thin_walled": "Timoshenko, thin-walled section",
              "gbt": f"GBT, {GBT_MODES} modes (Brazier)"}
    order = ("thin_walled", "2d", "gbt")
    fig, axes = plt.subplots(1, 4, figsize=(17, 4.2))
    panels = [(2, "flapwise tip deflection  $u_z / L$", 1.0 / L), (0, "tip shortening  $-u_x / L$", -1.0 / L),
              (3, "tip twist  $\\theta_x$ (rad)", 1.0)]
    for ax, (j, title, f) in zip(axes, panels):
        lin = ref["linear"]["tip"][0][j] * f
        ax.plot([0, lv[-1]], [0, lv[-1] * lin], "--", color="#c3c2b7", lw=1.0, label="solid, linear")
        ax.plot(lv, s_tip[:, j] * f, "o", ms=7, mfc="none", color="#222222", label="solid C3D20R, NLGEOM")
        for sec in order:
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
        if "gbt" in beams:
            ax.plot(lv, 100 * np.array(beams["gbt"]["flattening"][k]) / depth, m, ms=4, color=colours["gbt"],
                    label=f"GBT, x = {float(k):g} L")
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
    for sec in order:
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


def main(plot: bool = True, work_dir: str = "", out_dir: str = HERE, sections=("2d", "thin_walled", "gbt"),
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
    old = os.path.join(out_dir, "summary.json")
    if os.path.exists(old):                     # keep the beams that are not rerun
        with open(old, encoding="utf-8") as f:
            beams.update({k: v for k, v in json.load(f)["beams"].items() if k not in sections})
    for sec in sections:
        if sec == "gbt":
            beams[sec] = run_gbt(work, ref["nlgeom"]["levels"])
            continue
        tips, top = [], None
        for lv in ref["nlgeom"]["levels"]:
            U = run_job(work, sec, lv)
            tips.append(U[-1].tolist())
            top = U
        lin = run_job(work, sec, 1.0, linear=True)
        beams[sec] = {"tip": tips, "linear_tip": lin[-1].tolist(), "top_nodes": top.tolist()}
    summary = {"reference": {k: ref[k] for k in ("element", "n_along", "elements_along_wall", "length", "P")},
               "beam": {"element": "CorotationalTimoshenkoBeamElement3D", "gbt_element": "CorotationalGBTBeamElement3D",
                        "gbt_modes": GBT_MODES, "n_el": N_EL},
               "comparison": compare(beams, ref), "beams": beams}
    with open(os.path.join(out_dir, "summary.json"), "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=1)
    if plot:
        figures(beams, ref, out_dir)
    return summary


if __name__ == "__main__":
    only = [a.split("=", 1)[1].split(",") for a in sys.argv if a.startswith("--sections=")]
    s = main(plot_only="--plot-only" in sys.argv, **({"sections": tuple(only[0])} if only else {}))
    for sec, c in s["comparison"].items():
        print(sec, "linear", c.get("linear"), "flattening / solid", c.get("flattening_ratio"))
        for r in c["levels"]:
            print(f"  {r['level']:5.2f}  w/L={r['tip_deflection_over_span']:.3f}  u_z {100 * r['u_z']:+.2f} %  "
                  f"u_x {100 * r['u_x']:+.2f} %  theta_x {100 * r['theta_x']:+.2f} %")

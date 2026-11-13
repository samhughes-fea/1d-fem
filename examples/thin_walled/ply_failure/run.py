# examples/thin_walled/ply_failure/run.py
"""Ply failure margins of a composite box spar at the beam's Gauss points.

Run from the repository root::

    python examples/thin_walled/ply_failure/run.py

The box spar of ``box_spar`` (0.6 x 0.3 m, 12-ply carbon caps, +-45 webs, 6 m) carries a
flapwise tip force and a torque. With ``ply_strength.txt`` in the job, 1d-fem checks every
ply at every Gauss point of every element by four criteria (maximum stress, Tsai-Wu,
Hashin, Puck) and writes ``ply_failure.csv``. The run is linear, so the smallest reserve
factor ``1 / exposure`` is the first-ply-failure load factor. The same spar is solved with
thin-walled sections and with the section meshed in 2D.

Writes ``jobs/`` and ``results/`` (not tracked), and the figures and ``summary.json`` next
to this file.
"""

from __future__ import annotations

import csv
import json
import logging
import os
import sys
from typing import Dict, List

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
for p in (REPO, os.path.dirname(HERE)):
    if p not in sys.path:
        sys.path.insert(0, p)

from common import CantileverSpec, Material, SectionSpec, build_analysis, results_dir, write_job, write_supports  # noqa: E402

from onedfem.pre_processing.parsing.ply_strength_parser import format_ply_strength_file  # noqa: E402
from anisotropic_section_model.ply_failure import CRITERIA, PlyStrength  # noqa: E402

CFRP = Material(E1=135e9, E2=10e9, G12=5e9, nu12=0.3, G13=5e9, G23=3.8e9, rho=1550.0)
# carbon-epoxy strengths (T300/914-like) and Puck's carbon parameters
STRENGTH = PlyStrength(Xt=1500e6, Xc=1200e6, Yt=50e6, Yc=200e6, S12=70e6, p_TL_t=0.35, p_TL_c=0.30, p_TT=0.25)
B, H, T_PLY = 0.6, 0.3, 0.5e-3
L, N_EL = 6.0, 24
LOAD = [0.0, 0.0, 20e3, 5e3, 0.0, 0.0]      # flapwise tip force (z) and torque, N and N m
NAMES = {"max_stress": "maximum stress", "tsai_wu": "Tsai-Wu", "hashin": "Hashin", "puck": "Puck"}
COLOURS = {"max_stress": "#8b8a83", "tsai_wu": "#2a78d6", "hashin": "#eb6834", "puck": "#1baf7a"}
WALLS = {"1": "bottom cap", "2": "rear web", "3": "top cap", "4": "front web"}


def section() -> SectionSpec:
    cap = [(1, a, T_PLY) for a in (45, -45, 0, 0, 0, 0, 0, 0, 0, 0, -45, 45)]   # 6 mm
    web = [(1, a, T_PLY) for a in (45, -45, -45, 45)]                          # 2 mm
    return SectionSpec(
        points={1: (-B / 2, -H / 2), 2: (B / 2, -H / 2), 3: (B / 2, H / 2), 4: (-B / 2, H / 2)},
        walls=[(1, 1, 2, 1), (2, 2, 3, 2), (3, 3, 4, 1), (4, 4, 1, 2)],
        materials={1: CFRP}, laminates={1: cap, 2: web}, n_sub=4)


def run_case(work: str, kind: str) -> List[dict]:
    """Solve the spar with ``kind`` sections ("thin_walled" or "2d"); rows of ``ply_failure.csv``."""
    from onedfem.workflow_orchestrator.run_job import process_job

    spec = CantileverSpec(L, N_EL, "LinearTimoshenkoBeamElement3D", section(), LOAD)
    job = os.path.join(work, "jobs", kind)
    write_job(job, spec)
    write_supports(job, spec)
    if kind == "2d":
        from onedfem.interfaces.blade.solid_sections import write_section_mesh_2d
        from anisotropic_section_model.solid.mesh import SectionMesh

        an = build_analysis(spec.section)
        mesh = SectionMesh.from_thin_walled(an.section, an.laminates, elements_along=12)
        write_section_mesh_2d(os.path.join(job, "section_mesh_2d.txt"), {1: mesh})
        os.remove(os.path.join(job, "section_geometry.txt"))
    with open(os.path.join(job, "ply_strength.txt"), "w", encoding="utf-8") as f:
        f.write(format_ply_strength_file({1: STRENGTH}))
    res = results_dir(os.path.join(work, "results", kind))
    process_job(job, res, {}, {}, force_serial=True, max_processes_per_job=1)
    with open(os.path.join(res, "secondary_results", "ply_failure", "ply_failure.csv"), newline="",
              encoding="utf-8") as f:
        return list(csv.DictReader(f))


def first_ply_failure(rows: List[dict]) -> Dict[str, dict]:
    """Per criterion: the largest exposure, its load factor (``1 / exposure``), mode and place."""
    out = {}
    for c in CRITERIA:
        r = max(rows, key=lambda r: float(r[c]))
        f = float(r[c])
        mode = r.get(f"{c}_mode", "") or ""
        place = WALLS.get(r["wall_id"], f"section element {r['section_element']}") if r["wall_id"] else \
            f"section element {r['section_element']}"
        out[c] = {"exposure": f, "load_factor": 1.0 / f, "mode": mode, "x": float(r["x"]), "place": place,
                  "ply": r["ply"], "y": float(r["y"]), "z": float(r["z"])}
    return out


def along_span(rows: List[dict]) -> Dict[str, np.ndarray]:
    """Per Gauss point (sorted by x): the largest exposure over the section, per criterion."""
    by_x: Dict[float, Dict[str, float]] = {}
    for r in rows:
        d = by_x.setdefault(float(r["x"]), {c: 0.0 for c in CRITERIA})
        for c in CRITERIA:
            d[c] = max(d[c], float(r[c]))
    xs = np.array(sorted(by_x))
    return {"x": xs, **{c: np.array([by_x[x][c] for x in xs]) for c in CRITERIA}}


def contour_coordinates(r: dict) -> tuple:
    """``(contour distance, n)`` of a row: from its wall and ``s`` (thin-walled rows), or from
    the nearest wall mid-line of the box to its ``(y, z)`` (2D section rows)."""
    if r["wall_id"]:
        w = r["wall_id"]
        return _OFFSETS[w] + float(r["s"]), float(r["n"])
    y, z = float(r["y"]), float(r["z"])
    # walls: 1 bottom (+y), 2 rear web (+z), 3 top (-y), 4 front web (-z); n = x cross s points inward
    d = {"1": abs(z + H / 2), "2": abs(y - B / 2), "3": abs(z - H / 2), "4": abs(y + B / 2)}
    w = min(d, key=d.get)
    s_, n = {"1": (y + B / 2, z + H / 2), "2": (z + H / 2, B / 2 - y),
             "3": (B / 2 - y, H / 2 - z), "4": (H / 2 - z, y + B / 2)}[w]
    return _OFFSETS[w] + s_, n


_OFFSETS = {"1": 0.0, "2": B, "3": B + H, "4": 2 * B + H}


def figures(tw: List[dict], so: List[dict], out_dir: str) -> None:
    """Exposure along the span at every Gauss point, and maps of the exposure over every ply
    face (thin-walled) or quadrature point (2D) of the root Gauss point's section, unrolled."""
    from onedfem.post_processing.graphical_visualisers import figure_style as fs

    fs.use()
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(10.5, 3.8), sharey=True)
    for ax, rows, title in ((axes[0], tw, "Thin-walled sections"), (axes[1], so, "2D sections")):
        a = along_span(rows)
        for c in CRITERIA:
            fs.gauss(ax, a["x"], a[c], color=COLOURS[c], label=NAMES[c], ms=3.4)
        ax.set(xlabel="span x at the beam Gauss points (m)")
        fs.panel_title(ax, title, "largest exposure over the section at each Gauss point")
    axes[0].set_ylabel("exposure (1 / first-ply-failure load factor)")
    axes[0].legend(loc="upper right")
    fs.marks_legend(axes[1], ("gauss",), loc="upper right", labels={"gauss": "beam Gauss point (2 per element)"})
    fig.tight_layout()
    fs.save(fig, os.path.join(out_dir, "exposure_along_span.png"))

    # the root Gauss point: exposure on every ply face / quadrature point, around the unrolled contour
    fig, axes = plt.subplots(2, 2, figsize=(12, 5.6), sharex=True, sharey="row")
    for i, (rows, label) in enumerate(((tw, "thin-walled: every ply face at points along each wall"),
                                       (so, "2D section: every quadrature point of the section mesh"))):
        g0 = min(float(r["x"]) for r in rows)
        root = [r for r in rows if float(r["x"]) == g0]
        cn = np.array([contour_coordinates(r) for r in root])
        for j, c in enumerate(("tsai_wu", "puck")):
            ax = axes[i, j]
            e = np.array([float(r[c]) for r in root])
            order = np.argsort(e)
            sc = ax.scatter(cn[order, 0], cn[order, 1] * 1e3, c=e[order], s=6 if i == 0 else 3, cmap="viridis",
                            vmin=0, vmax=max(e.max(), 1e-12), marker="s", lw=0)
            k = int(np.argmax(e))
            ax.plot(cn[k, 0], cn[k, 1] * 1e3, "o", mfc="none", mec="#d4570b", ms=9, mew=1.4)
            ax.annotate(f"{e[k]:.3f}", (cn[k, 0], cn[k, 1] * 1e3), xytext=(8, 6), textcoords="offset points",
                        fontsize=7, color="#d4570b")
            for w, o in _OFFSETS.items():
                ax.axvline(o, color=fs.MUTED, lw=0.6)
                if i == 0:
                    ax.text(o + 0.02, 1.0, WALLS[w], transform=ax.get_xaxis_transform(), fontsize=6.5,
                            color=fs.MUTED, va="bottom")
            ax.grid(False)
            fig.colorbar(sc, ax=ax, pad=0.01, label=f"{NAMES[c]} exposure")
            if i == 1:
                ax.set_xlabel("contour distance (m)")
            if j == 0:
                ax.set_ylabel("n through the wall (mm)")
            ax.set_title(f"{NAMES[c]}, {label}", fontsize=8, pad=14 if i == 0 else 6)
    fig.suptitle(f"Root Gauss point (x = {min(float(r['x']) for r in tw):.3f} m): exposure over the unrolled "
                 "section; circle = worst point", x=0.01, ha="left")
    fig.tight_layout()
    fs.save(fig, os.path.join(out_dir, "root_contour.png"))


def main(plot: bool = True, work_dir: str = HERE, out_dir: str = HERE) -> dict:
    logging.disable(logging.WARNING)
    try:
        tw = run_case(work_dir, "thin_walled")
        so = run_case(work_dir, "2d")
    finally:
        logging.disable(logging.NOTSET)
    result = {"load": LOAD, "strength": STRENGTH.__dict__,
              "thin_walled": first_ply_failure(tw), "2d": first_ply_failure(so)}
    if plot:
        figures(tw, so, out_dir)
    with open(os.path.join(out_dir, "summary.json"), "w", encoding="utf-8") as f:
        json.dump(result, f, indent=1)
    return result


if __name__ == "__main__":
    print(json.dumps(main(), indent=1))

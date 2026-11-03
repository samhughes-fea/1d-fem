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

from pre_processing.parsing.ply_strength_parser import format_ply_strength_file  # noqa: E402
from pre_processing.section_integration.ply_failure import CRITERIA, PlyStrength  # noqa: E402

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
    from workflow_orchestrator.run_job import process_job

    spec = CantileverSpec(L, N_EL, "LinearTimoshenkoBeamElement3D", section(), LOAD)
    job = os.path.join(work, "jobs", kind)
    write_job(job, spec)
    write_supports(job, spec)
    if kind == "2d":
        from pre_processing.section_integration.solid.from_geometry import write_section_mesh_2d
        from pre_processing.section_integration.solid.mesh import SectionMesh

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


def figures(tw: List[dict], so: List[dict], out_dir: str) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(10, 3.6), sharey=True)
    for ax, rows, title in ((axes[0], tw, "thin-walled sections"), (axes[1], so, "2D sections")):
        a = along_span(rows)
        for c in CRITERIA:
            ax.plot(a["x"], a[c], "o", ms=3.5, color=COLOURS[c], label=NAMES[c])
        ax.set_xlabel("span x at the Gauss points, m", fontsize=8)
        ax.set_title(title, fontsize=9, loc="left")
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        ax.tick_params(labelsize=7)
    axes[0].set_ylabel("largest exposure in the section", fontsize=8)
    axes[0].legend(fontsize=7, frameon=False)
    fig.suptitle("Ply failure exposure at every Gauss point (1 / exposure = first-ply-failure load factor)",
                 fontsize=9, x=0.01, ha="left")
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, "exposure_along_span.png"), dpi=150)
    plt.close(fig)

    # the root Gauss point of the thin-walled model, around the contour, per ply
    g0 = min(float(r["x"]) for r in tw)
    root = [r for r in tw if float(r["x"]) == g0]
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.6), sharey=True)
    offsets = {"1": 0.0, "2": B, "3": B + H, "4": 2 * B + H}
    for ax, c in zip(axes, ("tsai_wu", "puck")):
        for r in root:
            s = offsets[r["wall_id"]] + float(r["s"])
            ax.plot(s, float(r[c]), "o", ms=2.5, color=COLOURS[c], alpha=0.6)
        for w, o in offsets.items():
            ax.axvline(o, color="#e9e8e3", lw=0.8, zorder=0)
            ax.text(o + 0.02, 0.98, WALLS[w], transform=ax.get_xaxis_transform(), fontsize=6, color="#52514e",
                    va="top")
        ax.set_xlabel("contour coordinate (bottom cap, rear web, top cap, front web), m", fontsize=8)
        ax.set_title(f"{NAMES[c]}, every ply face at the root Gauss point", fontsize=9, loc="left")
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        ax.tick_params(labelsize=7)
    axes[0].set_ylabel("exposure", fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, "root_contour.png"), dpi=150)
    plt.close(fig)


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

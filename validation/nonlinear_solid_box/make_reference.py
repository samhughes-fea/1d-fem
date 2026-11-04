"""CalculiX ``NLGEOM`` reference for the corotational beams: the CAS composite box in solids.

The bend-twist coupled CFRP box of the shell validation (``abaqus/shell_validation/cases.py``,
case ``composite_box``: 0.1 m square, eight 0.5 mm plies per wall, 6 m long) is meshed in 2D
(``SectionMesh.from_thin_walled``, one ``Q8`` row per ply, walls tied at the corners) and
extruded into ``C3D20R`` bricks (``iea15mw_sections/prism.py``). The root face is clamped and
the tip face is tied to a reference node (``*RIGID BODY``) that carries a dead flapwise force.

The force ``P`` deflects the tip 20 % of the span in linear theory. One ``NLGEOM`` step takes
it to ``LEVELS[-1] P``, with the tip motion printed at every level; a second deck without
``NLGEOM`` gives the linear solution at ``P``. Rotations of the reference node are total
rotation vectors, as in the corotational beams.

Solid models need no shell offsets or rotational junction ties, which is what stopped the
``NLGEOM`` shell models converging. Writes ``reference.json`` here. Needs ``ccx`` on the PATH
(a few minutes)::

    python validation/nonlinear_solid_box/make_reference.py
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
import tempfile
from typing import Dict, List

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
for p in (REPO, os.path.join(REPO, "validation", "iea15mw_sections")):
    if p not in sys.path:
        sys.path.insert(0, p)

N_ALONG = 60                       # brick layers along the 6 m box
ELEMENTS_ALONG_WALL = 6            # Q8 elements along each wall of the section
LEVELS = (0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0)    # load levels, multiples of P
STATIONS = (0.25, 0.5, 0.75)       # span fractions where the mean section motion is printed


def box_mesh():
    """The composite box section as a conforming layered 2D mesh with mitred corners."""
    from validation.abaqus.shell_validation.cases import box_section
    from validation.abaqus.shell_validation.onedfem import section_analysis
    from onedfem.pre_processing.section_integration.solid.mesh import SectionMesh

    tw = section_analysis(box_section())
    return SectionMesh.from_thin_walled(tw.section, tw.laminates, elements_along=ELEMENTS_ALONG_WALL,
                                       junctions="mitre")


def _nset(name: str, ids: List[int]) -> List[str]:
    out = [f"*NSET, NSET={name}"]
    for a in range(0, len(ids), 12):
        out.append(", ".join(map(str, ids[a:a + 12])) + ",")
    return out


def _steps(P: float, nlgeom: bool, station_sets: List[str]) -> List[str]:
    top = LEVELS[-1] if nlgeom else 1.0
    out = []
    if nlgeom:
        out += ["*TIME POINTS, NAME=TP", ", ".join(f"{lv / top:.10g}" for lv in LEVELS)]
        out += ["*STEP, NLGEOM, INC=400", "*STATIC", "0.02, 1.0, 1e-5, 0.05"]
    else:
        out += ["*STEP", "*STATIC"]
    out += ["*CLOAD", f"REF, 3, {top * P:.10e}"]
    tp = ", TIME POINTS=TP" if nlgeom else ""
    out += [f"*NODE PRINT, NSET=REF{tp}", "U", f"*NODE PRINT, NSET=ROT{tp}", "U"]
    for s in station_sets:
        out += [f"*NODE PRINT, NSET={s}{tp}", "U"]
    out += ["*END STEP"]
    return out


def _read_blocks(dat: str, nset: str) -> List[Dict[int, np.ndarray]]:
    """Every ``displacements ... for set NSET`` block of a ``.dat`` file, in order."""
    blocks = []
    lines = dat.splitlines()
    for i, ln in enumerate(lines):
        if ln.strip().startswith("displacements") and f"set {nset}" in ln:
            vals = {}
            for row in lines[i + 2:]:
                parts = row.split()
                if len(parts) != 4:
                    break
                vals[int(parts[0])] = np.array([float(v) for v in parts[1:]])
            blocks.append(vals)
    return blocks


def flattening(mesh, nid: Dict, k: int, block: Dict[int, np.ndarray], half_width: float,
               half_depth: float) -> float:
    """Brazier flattening of the section at level ``k``: how much closer the flanges' mid-points
    come than their corners, after the section's rigid motion is taken out (a best-fit
    rotation of all its nodes)."""
    ids = [(i, nid[(k, i)]) for i in range(mesh.nodes.shape[0]) if nid[(k, i)] in block]
    yz = np.array([mesh.nodes[i] for i, _ in ids])
    P = np.c_[np.zeros(len(ids)), yz]
    x = P + np.array([block[n] for _, n in ids])
    c0, c1 = P.mean(0), x.mean(0)
    u, _, vt = np.linalg.svd((P - c0).T @ (x - c1))
    loc = (x - c1) @ (vt.T @ u.T)
    y, z = yz[:, 0], yz[:, 1]
    tol = 1e-3 * half_depth
    on_flange = np.abs(np.abs(z) - half_depth) < 0.05 * half_depth
    mid = on_flange & (np.abs(y) < tol)
    corner = on_flange & (np.abs(np.abs(y) - half_width) < tol)

    def gap(sel):
        return loc[sel & (z > 0), 2].mean() - loc[sel & (z < 0), 2].mean()

    return float(gap(corner) - gap(mid))


def run_deck(mesh, P: float, nlgeom: bool, work_dir: str) -> Dict:
    import prism

    from validation.abaqus.shell_validation.cases import BOX_H, BOX_L, BOX_W

    name = "nlgeom" if nlgeom else "linear"
    deck = os.path.join(work_dir, f"{name}.inp")
    # station node sets: every node of the section at the station level (even levels only)
    probe = prism.write_deck(mesh, BOX_L, N_ALONG, deck, (0, 1), steps=[])
    nid = probe["nid"]
    station_sets, extra = [], []
    for s in STATIONS:
        k = 2 * int(round(s * N_ALONG))
        ids = [nid[(k, i)] for i in range(probe["n2"]) if i not in mesh.ties]
        station_sets.append(f"S{int(round(100 * s))}")
        extra += _nset(station_sets[-1], ids)
    info = prism.write_deck(mesh, BOX_L, N_ALONG, deck, (0, 1), steps=extra + _steps(P, nlgeom, station_sets))
    subprocess.run(["ccx", name], cwd=work_dir, check=True, stdout=subprocess.DEVNULL)
    dat = open(os.path.join(work_dir, f"{name}.dat"), encoding="utf-8").read()
    tip_u = [b[info["ref"]] for b in _read_blocks(dat, "REF")]
    tip_r = [b[info["rot"]] for b in _read_blocks(dat, "ROT")]
    levels = list(LEVELS) if nlgeom else [1.0]
    if len(tip_u) != len(levels):
        raise RuntimeError(f"{name}: CalculiX printed {len(tip_u)} of {len(levels)} load levels")
    out = {"levels": levels, "tip": [list(map(float, np.r_[u, r])) for u, r in zip(tip_u, tip_r)],
           "stations": {}}
    out["flattening"] = {}
    for s, name_s in zip(STATIONS, station_sets):
        blocks = _read_blocks(dat, name_s)
        out["stations"][f"{s:g}"] = [list(map(float, np.mean(list(b.values()), axis=0))) for b in blocks]
        k = 2 * int(round(s * N_ALONG))
        out["flattening"][f"{s:g}"] = [flattening(mesh, nid, k, b, 0.5 * BOX_W, 0.5 * BOX_H) for b in blocks]
    return out


def main(work_dir: str = "") -> dict:
    from validation.abaqus.shell_validation.cases import BOX_L, box_tip_force

    mesh = box_mesh()
    P = box_tip_force()
    work_dir = work_dir or tempfile.mkdtemp(prefix="nl_box_")
    out = {"element": "C3D20R", "n_along": N_ALONG, "elements_along_wall": ELEMENTS_ALONG_WALL,
           "length": BOX_L, "P": P,
           "note": "tip: [u_x, u_y, u_z, theta_x, theta_y, theta_z] of the rigid tip reference node "
                   "(total rotation vector); stations: mean [u_x, u_y, u_z] of the section nodes",
           "linear": run_deck(mesh, P, False, work_dir),
           "nlgeom": run_deck(mesh, P, True, work_dir)}
    with open(os.path.join(HERE, "reference.json"), "w", encoding="utf-8") as f:
        json.dump(out, f, indent=1)
    return out


if __name__ == "__main__":
    res = main(sys.argv[1] if len(sys.argv) > 1 else "")
    print("linear tip", res["linear"]["tip"][0])
    for lv, t in zip(res["nlgeom"]["levels"], res["nlgeom"]["tip"]):
        print(lv, t)

"""Shell reference for the blade_beam example: the IEA outboard segment as a
composite-laminate-model shell blade.

Needs composite-laminate-model, implicit-geometry-kernel and CalculiX (``ccx``), none of
which 1d-fem depends on, so the output, ``data/shell_reference.json``, is committed. To
regenerate it (about 15 minutes)::

    python examples/blade_beam/make_shell_reference.py

The blade is meshed from the same IGK loft as the beams (``examples/iea15mw/segment.py``)
with curved Quad8 FSDT shells on the walls' mid-surfaces, 40 elements along the span, and
clamped at the root. Three analyses:

- **static** (composite-laminate-model's own solver): 10 kN flap and edge tip forces and a
  10 kN m torque through a rigid tip (RBE2 on the translations), for the tapered,
  twisted segment and for a prismatic, untwisted one;
- **modal** (CalculiX ``*FREQUENCY`` of the same shell deck): eight frequencies and the
  tip motion of each mode;
- **nonlinear** (CalculiX ``NLGEOM``): a flapwise tip force of up to 1.8 MN spread over the
  web tip nodes, in 5 % increments, with the mean tip motion at each.
"""

from __future__ import annotations

import json
import os
import re
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, os.path.join(REPO, "examples", "iea15mw"))
import segment as G  # noqa: E402
import stations as S  # noqa: E402

N_SPAN = 40
STATIC = {"flap": (1, 1e4), "edge": (0, 1e4), "torque": (5, 1e4)}
P_NONLINEAR = 1.8e6


def build(prismatic: bool = False):
    from composite_laminate_model import BladeBuilder, FSDTModel, Laminate, Ply
    from composite_laminate_model.laminate_materials import MaterialCard
    from composite_laminate_model.laminate_materials_3d import ElasticProperties3D, StrengthProperties3D

    big = 1e12      # strengths are not used
    cards = {n: MaterialCard(n, ElasticProperties3D.transversely_isotropic(p["E1"], p["E2"], p["nu12"], p["G12"], S.NU23),
                             StrengthProperties3D(big, big, big, big, big), density=p["rho"])
             for n, p in S.DATA["materials"].items()}
    lams = {k: Laminate([Ply(cards[n], t, 0.0) for n, t in v]) for k, v in G.MODEL.laminae.items()}
    return BladeBuilder.from_geometry(G.geometry(prismatic), n_span=N_SPAN,
                                      element_size=0.03 * G.TIP["chord"]).build(lams, FSDTModel())


def static(prismatic: bool = False) -> dict:
    from composite_laminate_model.laminate_connectors import add_reference_node, rigid_link
    from composite_laminate_model.laminate_solver_pipeline import BoundaryConditions, LoadCase, StructuralGlobalSolver

    blade = build(prismatic)
    mesh = blade.mesh.copy()
    bcs = BoundaryConditions().clamp(blade.nodes("root"))
    tip = add_reference_node(mesh, (0.0, 0.0, G.LENGTH))
    rigid_link(mesh, bcs, tip, blade.nodes("tip"), dofs=(0, 1, 2))
    cases = [LoadCase(k).add_nodal_load(tip, d, v) for k, (d, v) in STATIC.items()]
    sols = StructuralGlobalSolver().solve_many(mesh, bcs, cases)
    return {"tip": {k: s.node_displacement(tip).tolist() for k, s in zip(STATIC, sols)},
            "loads": {k: list(v) for k, v in STATIC.items()}, "n_elements": len(blade.mesh.element_ids)}


def _ccx(mode: str, work: Path) -> dict:
    from composite_laminate_model.kernels.base import AnalysisModel
    from composite_laminate_model.kernels.calculix import CalculiXDeckWriter
    from composite_laminate_model.laminate_solver_pipeline import BoundaryConditions, LoadCase

    blade = build()
    bcs = BoundaryConditions().clamp(blade.nodes("root"))
    tipn = sorted(blade.nodes("tip"))
    webs = [w for w in blade.node_sets if w.startswith("shear_web") or w.startswith("web")]
    web_tip = sorted(set(tipn) & {n for w in webs for n in blade.nodes(w)})
    load = LoadCase("lc")
    if mode == "nonlinear":
        for n in web_tip:
            load.add_nodal_load(n, 1, P_NONLINEAR / len(web_tip))
    else:
        load.add_nodal_load(tipn[0], 1, 0.0)
    work.mkdir(parents=True, exist_ok=True)
    deck = CalculiXDeckWriter().write(AnalysisModel(blade.mesh, bcs, load, title=f"seg_{mode}"), work, "seg")
    head, step = deck.path.read_text().split("*STEP\n*STATIC\n", 1)
    head += "*NSET, NSET=TIPN\n" + "\n".join(", ".join(map(str, tipn[i:i + 12])) + "," for i in range(0, len(tipn), 12)) + "\n"
    cloads = step.split("*NODE PRINT")[0]
    if mode == "modal":
        step = "*STEP\n*FREQUENCY\n8\n*NODE PRINT, NSET=TIPN\nU\n*END STEP\n"
    else:
        step = ("*STEP, NLGEOM, INC=400\n*STATIC\n0.05, 1.0, 1e-6, 0.05\n" + cloads +
                "*NODE PRINT, NSET=TIPN, FREQUENCY=1\nU\n*END STEP\n")
    deck.path.write_text(head + step)
    r = subprocess.run(["ccx", "seg"], cwd=work, capture_output=True, text=True,
                       env=dict(os.environ, OMP_NUM_THREADS=str(os.cpu_count() or 1)))
    (work / "ccx.log").write_text(r.stdout + r.stderr)
    dat = (work / "seg.dat").read_text()
    X = np.array([blade.mesh.node_coordinates(n) for n in tipn])
    series = []
    for b in re.split(r"\n\s*displacements \(vx,vy,vz\) for set TIPN and time\s+", dat)[1:]:
        lines = b.strip().splitlines()
        U = {}
        for ln in lines[1:]:
            p = ln.split()
            if len(p) == 4:
                try:
                    U[int(p[0])] = [float(v) for v in p[1:]]
                    continue
                except ValueError:
                    pass
            if U:
                break
        U = np.array([U[n] for n in tipn])
        r_, u_ = X[:, :2] - X[:, :2].mean(0), U[:, :2] - U[:, :2].mean(0)
        twist = float(np.sum(r_[:, 0] * u_[:, 1] - r_[:, 1] * u_[:, 0]) / np.sum(r_ ** 2))
        series.append([float(lines[0].split()[0]), *U.mean(0).tolist(), twist])
    out = {"tip_series": series, "ccx_returncode": r.returncode}
    if mode == "modal":
        freqs = []
        m = re.search(r"E I G E N V A L U E   O U T P U T(.*?)P A R T I C", dat, re.S)
        for ln in (m.group(1) if m else "").splitlines():
            p = ln.split()
            if len(p) == 5 and p[0].isdigit():
                freqs.append(float(p[3]))
        out["frequencies_hz"] = freqs
    else:
        out["P"] = P_NONLINEAR
    return out


def main():
    t0 = time.time()
    work = Path(HERE) / "output"
    ref = {"length": G.LENGTH, "n_span": N_SPAN,
           "static": static(), "static_prismatic": static(prismatic=True),
           "modal": _ccx("modal", work / "modal"), "nonlinear": _ccx("nonlinear", work / "nonlinear")}
    ref["time"] = time.time() - t0
    os.makedirs(os.path.join(HERE, "data"), exist_ok=True)
    with open(os.path.join(HERE, "data", "shell_reference.json"), "w", encoding="utf-8") as f:
        json.dump(ref, f, indent=1)
    print("done in", round(ref["time"]), "s")


if __name__ == "__main__":
    main()

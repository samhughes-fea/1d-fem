"""A 3D solid check of a 2D section: the section extruded into a CalculiX ``C3D20R`` prism.

The 2D section mesh (``Q8``) is extruded along ``x`` into 20-node bricks, each with its
ply's material in its ply axes. The root face is clamped and the tip face is tied to a
reference point at the section origin (``*RIGID BODY``). Unit tip loads in ``FORCES`` order
(``N``, ``M_y``, ``M_z``, ``V_y``, ``V_z``, ``T``) and their 15 pairs are applied one per
step, and CalculiX prints the strain energy of the middle part of the prism, away from
the end zones, for each.

The energy of a beam segment ``[x1, x2]`` is ``U = 1/2 int F(x)^T S F(x) dx`` with the
section forces ``F(x) = A(x) P`` of the tip load ``P`` (the shear loads carry a moment that
grows towards the root). The 21 energies fix the 21 terms of the flexibility ``S``, and
its inverse is the 6 x 6 stiffness, couplings included. It is the stiffness the 2D
section analysis computes, so the two can be compared term by term. (Fitting a rigid
motion to the nodes of two sections instead gives shear terms 10-45 % off, since shear
warps the section.)
"""

from __future__ import annotations

import os
import subprocess
from typing import Dict, Tuple

import numpy as np

from pre_processing.section_integration.solid.material import material_axes
from pre_processing.section_integration.solid.mesh import SectionMesh

LOADS = np.eye(6)          # unit load cases in FORCES order: N, M_y, M_z, V_y, V_z, T


def load_cases() -> np.ndarray:
    """The six unit loads, then the 15 pairs of them (for the energy fit)."""
    cases = [LOADS[i] for i in range(6)]
    cases += [LOADS[i] + LOADS[j] for i in range(6) for j in range(i + 1, 6)]
    return np.array(cases)


def _fmt(v: float) -> str:
    return f"{v:.10e}"


def write_deck(mesh: SectionMesh, length: float, n_along: int, path: str, stations: Tuple[int, int]) -> Dict:
    """Write the ``.inp`` deck; return the node bookkeeping needed to read the results.

    ``stations`` are two element layer indices (0..n_along) bounding the elements whose
    energy is printed.
    """
    nodes2d = mesh.nodes
    n2 = nodes2d.shape[0]
    corner = np.zeros(n2, dtype=bool)
    for conn in mesh.elements:
        corner[conn[:4]] = True
    # node ids: level k (0..2 n_along); even levels hold every 2D node, odd levels only corners
    nid = {}
    lines = []
    next_id = 1
    for k in range(2 * n_along + 1):
        x = length * k / (2 * n_along)
        for i in range(n2):
            if k % 2 == 1 and not corner[i]:
                continue
            nid[(k, i)] = next_id
            y, z = nodes2d[i]
            lines.append(f"{next_id}, {_fmt(x)}, {_fmt(y)}, {_fmt(z)}")
            next_id += 1
    ref, rot = next_id, next_id + 1
    lines.append(f"{ref}, {_fmt(length)}, 0.0, 0.0")
    lines.append(f"{rot}, {_fmt(length)}, 0.0, 0.0")

    out = ["*HEADING", "2D section extruded into a C3D20R prism", "*NODE"] + lines
    # elements grouped by (material, theta, alpha)
    groups: Dict[tuple, list] = {}
    eid = 1
    elem_lines = []
    for e, conn in enumerate(mesh.elements):
        key = (mesh.material_ids[e], round(float(mesh.theta[e]), 6), round(float(mesh.alpha[e]), 6))
        for j in range(n_along):
            k0, k1, km = 2 * j, 2 * j + 2, 2 * j + 1
            c = [nid[(k0, i)] for i in conn[:4]] + [nid[(k1, i)] for i in conn[:4]] \
                + [nid[(k0, i)] for i in conn[4:8]] + [nid[(k1, i)] for i in conn[4:8]] \
                + [nid[(km, i)] for i in conn[:4]]
            elem_lines.append((key, eid, c))
            groups.setdefault(key, []).append(eid)
            eid += 1
    out.append("*ELEMENT, TYPE=C3D20R, ELSET=EALL")
    for _, e_, c in elem_lines:
        out.append(f"{e_}, " + ", ".join(map(str, c[:15])) + ",")
        out.append(", ".join(map(str, c[15:])))
    for g, (key, ids) in enumerate(groups.items()):
        out.append(f"*ELSET, ELSET=G{g}")
        for a in range(0, len(ids), 12):
            out.append(", ".join(map(str, ids[a:a + 12])) + ",")
    for mid, m in mesh.materials.items():
        out += [f"*MATERIAL, NAME=M{mid}", "*ELASTIC, TYPE=ENGINEERING CONSTANTS",
                ", ".join(_fmt(v) for v in (m.E1, m.E2, m.E3, m.nu12, m.nu13, m.nu23, m.G12, m.G13)) + ",",
                _fmt(m.G23)]
    for g, key in enumerate(groups):
        mid, th, al = key
        R = material_axes(th, al)
        a, b = R[:, 0], R[:, 1]
        out += [f"*ORIENTATION, NAME=O{g}, SYSTEM=RECTANGULAR",
                ", ".join(_fmt(v) for v in (*a, *b))]
        out += [f"*SOLID SECTION, ELSET=G{g}, MATERIAL=M{mid}, ORIENTATION=O{g}"]
    root = [nid[(0, i)] for i in range(n2)]
    tip = [nid[(2 * n_along, i)] for i in range(n2)]

    def nset(name, ids):
        res = [f"*NSET, NSET={name}"]
        for a in range(0, len(ids), 12):
            res.append(", ".join(map(str, ids[a:a + 12])) + ",")
        return res

    out += nset("ROOT", root) + nset("TIP", tip)
    out += ["*RIGID BODY, NSET=TIP, REF NODE=%d, ROT NODE=%d" % (ref, rot)]
    out += ["*BOUNDARY", "ROOT, 1, 3, 0.0"]
    mid = [e_ for key, e_, c in elem_lines if stations[0] <= (e_ - 1) % n_along < stations[1]]
    out.append("*ELSET, ELSET=MID")
    for a in range(0, len(mid), 12):
        out.append(", ".join(map(str, mid[a:a + 12])) + ",")
    for case, F in enumerate(load_cases()):
        N, My, Mz, Vy, Vz, T = F
        out += ["*STEP", "*STATIC", "*CLOAD"]
        out += [f"{ref}, 1, {_fmt(N)}", f"{ref}, 2, {_fmt(Vy)}", f"{ref}, 3, {_fmt(Vz)}",
                f"{rot}, 1, {_fmt(T)}", f"{rot}, 2, {_fmt(My)}", f"{rot}, 3, {_fmt(Mz)}"]
        out += ["*EL PRINT, ELSET=MID, TOTALS=ONLY", "ELSE"]
        out += ["*END STEP"]
    with open(path, "w", encoding="utf-8") as f:
        f.write("\n".join(out) + "\n")
    return {"nid": nid, "n2": n2}


def read_energies(path: str) -> np.ndarray:
    """Total strain energy of the printed element set, per step, from a CalculiX ``.dat``."""
    out, grab = [], False
    with open(path, encoding="utf-8", errors="replace") as f:
        for line in f:
            s = line.strip()
            if s.lower().startswith("total internal energy"):
                grab = True
                continue
            if grab and s:
                out.append(float(s.split()[0]))
                grab = False
    return np.array(out)


def energy_flexibility(cases: np.ndarray, energies: np.ndarray, length: float, x1: float, x2: float) -> np.ndarray:
    """Flexibility ``S`` (6 x 6) from the strain energy of ``[x1, x2]`` under tip loads ``cases``:
    ``U = 1/2 int F(x)^T S F(x) dx`` with ``F(x) = A(x) P``, the shear loads carrying their moment."""
    iu = np.triu_indices(6)
    xg, wg = np.polynomial.legendre.leggauss(3)
    xs = 0.5 * (x1 + x2) + 0.5 * (x2 - x1) * xg
    ws = 0.5 * (x2 - x1) * wg
    rows = []
    for P in cases:
        r = np.zeros(len(iu[0]))
        for x, w in zip(xs, ws):
            A = np.eye(6)
            A[2, 3] = length - x
            A[1, 4] = -(length - x)
            F = A @ P
            Q = np.outer(F, F)
            Q = Q + Q.T - np.diag(np.diag(Q))          # off-diagonal unknowns appear twice
            r += 0.5 * w * Q[iu]
        rows.append(r)
    sol = np.linalg.lstsq(np.array(rows), energies, rcond=None)[0]
    S = np.zeros((6, 6))
    S[iu] = sol
    return S + S.T - np.diag(np.diag(S))


def section_stiffness(mesh: SectionMesh, length: float, n_along: int, work_dir: str, ccx: str = "ccx",
                      stations=None) -> Dict:
    """Run the prism and return its 6 x 6 stiffness (``FORCES`` / ``STRAINS`` order) and details."""
    os.makedirs(work_dir, exist_ok=True)
    s1, s2 = stations or (int(round(n_along / 3)), int(round(2 * n_along / 3)))
    inp = os.path.join(work_dir, "prism.inp")
    book = write_deck(mesh, length, n_along, inp, (s1, s2))
    env = dict(os.environ, OMP_NUM_THREADS=str(os.cpu_count() or 1))
    with open(os.path.join(work_dir, "ccx.log"), "w", encoding="utf-8") as log:
        subprocess.run([ccx, "-i", "prism"], cwd=work_dir, check=True, env=env, stdout=log, stderr=subprocess.STDOUT)
    cases = load_cases()
    energies = read_energies(os.path.join(work_dir, "prism.dat"))
    if energies.size != len(cases):
        raise RuntimeError(f"expected {len(cases)} energies in prism.dat, found {energies.size}")
    x1, x2 = length * s1 / n_along, length * s2 / n_along
    S = energy_flexibility(cases, energies, length, x1, x2)
    return {"stiffness": np.linalg.inv(S), "flexibility": S, "x_range": (x1, x2),
            "n_nodes": len(book["nid"]), "n_elements": len(mesh.elements) * n_along}

"""Does the shell reference see the sandwich cores' shear? A strip of the blade skin, in CalculiX.

Needs CalculiX (``ccx``), so the output, ``data/sandwich_strip.json``, is committed. To
regenerate it (under a minute)::

    python examples/blade_beam/sandwich_strip_check.py

The shell reference (``make_shell_reference.py``) writes composite ``S8R`` shells, which
CalculiX expands into one ``C3D20R`` brick through the wall: the displacement through the
whole sandwich is one quadratic, so the faces cannot slide over the core and the core's
transverse shear costs nothing. The 2D sections mesh every ply, so their panels bend with
the core's shear. The check bends a strip of the skin laminate at 70 % span (gelcoat,
triax, foam, triax) across its width, simply supported over spans of 0.5, 1 and 2 m, and
compares the first bending frequency of

- the composite ``S8R`` strip, as the shell reference has it;
- ``C3D20`` bricks, two per ply through the thickness;
- a sandwich beam with and without the core's shear: ``f = f_b / sqrt(1 + (D / S) (pi / L)^2)``,
  ``S = G_c d^2 / c``.

The strip is free across its width (a narrow beam) in all three.
"""

from __future__ import annotations

import json
import os
import re
import subprocess
import sys
import tempfile

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
for p in (REPO, os.path.join(REPO, "examples", "iea15mw")):
    if p not in sys.path:
        sys.path.insert(0, p)
import stations as S  # noqa: E402

SPANS = (0.5, 1.0, 2.0)
N_ALONG = 40
BRICKS_PER_PLY = 2
OUT = os.path.join(HERE, "data", "sandwich_strip.json")


def _layup():
    """``[(thickness, name)]`` bottom ply first, and the materials, for the skin at 70 % span."""
    plies = S.station_model(20).laminae["skin"]
    return [(t, m) for m, t in plies], S.ply_materials()


def _materials(mats) -> str:
    """CalculiX materials; the strip runs along global ``x`` = the ply's direction 2 (chordwise),
    ``y`` = direction 1 (the fibres, spanwise), ``z`` through the wall."""
    out = ""
    for name, m in mats.items():
        nu21 = m.nu12 * m.E2 / m.E1
        nu23 = m.E2 / (2.0 * m.G23) - 1.0
        out += (f"*MATERIAL, NAME={name}\n*ELASTIC, TYPE=ENGINEERING CONSTANTS\n"
                f"{m.E2:.8g},{m.E1:.8g},{m.E2:.8g},{nu21:.8g},{nu23:.8g},{nu21:.8g},{m.G12:.8g},{m.G23:.8g},\n"
                f"{m.G13:.8g}\n*DENSITY\n{m.rho:.8g}\n")
    return out


def _bending(work, name, deck, n=6):
    """Frequencies of the modes that move the strip through the wall: their participation
    factor in ``z`` dominates (the symmetric bending modes, the first among them)."""
    with open(os.path.join(work, name + ".inp"), "w", encoding="utf-8") as f:
        f.write(deck + f"*STEP\n*FREQUENCY\n{n}\n*END STEP\n")
    subprocess.run(["ccx", name], cwd=work, capture_output=True, text=True, check=True)
    with open(os.path.join(work, name + ".dat"), encoding="utf-8") as f:
        dat = f.read()
    eig, part = dat.split("E I G E N V A L U E   O U T P U T")[1].split("P A R T I C I P A T I O N")
    num = r"\s+\d+(\s+-?[\d.]+E[+-]\d+)+\s*$"
    freqs = [float(r.split()[3]) for r in eig.splitlines() if re.match(num, r)][:n]
    P = np.array([[float(v) for v in r.split()[1:4]] for r in part.splitlines() if re.match(num, r)][:n])
    return [f for f, p in zip(freqs, np.abs(P)) if p[2] > 3 * max(p[0], p[1])]


def shell(work, L, layup, mats):
    nx = N_ALONG
    w = L / nx
    xs = np.linspace(0.0, L, 2 * nx + 1)
    nid = lambda i, j: 1 + 3 * i + j  # noqa: E731
    nodes = {nid(i, j): (x, y) for i, x in enumerate(xs) for j, y in enumerate((0.0, w / 2, w))
             if not (j == 1 and i % 2)}
    el = [[nid(i, 0), nid(i + 2, 0), nid(i + 2, 2), nid(i, 2), nid(i + 1, 0), nid(i + 2, 1), nid(i + 1, 2), nid(i, 1)]
          for i in range(0, 2 * nx, 2)]
    d = "*NODE\n" + "".join(f"{k},{x:.10g},{y:.10g},0\n" for k, (x, y) in nodes.items())
    d += "*ELEMENT, TYPE=S8R, ELSET=E\n" + "".join(f"{k + 1}," + ",".join(map(str, c)) + "\n" for k, c in enumerate(el))
    d += _materials(mats) + "*SHELL SECTION, ELSET=E, COMPOSITE, OFFSET=0\n"
    d += "".join(f"{t:.10g},,{m}\n" for t, m in layup)
    ends = [k for k, (x, _) in nodes.items() if x < 1e-12 or x > L * (1 - 1e-12)]
    d += "*BOUNDARY\n" + "".join(f"{k},3,3\n" for k in ends)
    d += f"{nid(0, 0)},1,2\n{nid(2 * nx, 0)},2,2\n"
    return _bending(work, f"shell_{L:g}", d)


def bricks(work, L, layup, mats):
    nx, nb = N_ALONG, BRICKS_PER_PLY
    H = sum(t for t, _ in layup)
    zs = [0.0]
    mat_of = []
    for t, m in layup:
        zs += [zs[-1] + t * k / (2 * nb) for k in range(1, 2 * nb + 1)]
        mat_of += [m] * nb
    zs = np.array(zs) - H / 2
    w = L / nx
    xs = np.linspace(0.0, L, 2 * nx + 1)
    nz = len(zs)
    nid = lambda i, j, k: 1 + (3 * i + j) * nz + k  # noqa: E731
    conn = {m: [] for m in mats}
    for i in range(0, 2 * nx, 2):
        for layer, m in enumerate(mat_of):
            k = 2 * layer
            c = lambda a, b, cc: nid(i + a, b, k + cc)  # noqa: E731
            conn[m].append([c(0, 0, 0), c(2, 0, 0), c(2, 2, 0), c(0, 2, 0), c(0, 0, 2), c(2, 0, 2), c(2, 2, 2), c(0, 2, 2),
                            c(1, 0, 0), c(2, 1, 0), c(1, 2, 0), c(0, 1, 0), c(1, 0, 2), c(2, 1, 2), c(1, 2, 2), c(0, 1, 2),
                            c(0, 0, 1), c(2, 0, 1), c(2, 2, 1), c(0, 2, 1)])
    used = {n for v in conn.values() for cn in v for n in cn}
    nodes = {nid(i, j, k): (x, y, z) for i, x in enumerate(xs) for j, y in enumerate((0.0, w / 2, w))
             for k, z in enumerate(zs) if nid(i, j, k) in used}
    d = "*NODE\n" + "".join(f"{k},{x:.10g},{y:.10g},{z:.10g}\n" for k, (x, y, z) in nodes.items())
    eid = 1
    for m, v in conn.items():
        if not v:
            continue
        d += f"*ELEMENT, TYPE=C3D20, ELSET=S{m}\n"
        for cn in v:
            d += f"{eid}," + ",".join(map(str, cn[:15])) + ",\n" + ",".join(map(str, cn[15:])) + "\n"
            eid += 1
    d += _materials(mats) + "".join(f"*SOLID SECTION, ELSET=S{m}, MATERIAL={m}\n" for m, v in conn.items() if v)
    k0 = int(np.argmin(np.abs(zs)))
    ends = [k for k, (x, _, z) in nodes.items() if (x < 1e-12 or x > L * (1 - 1e-12)) and abs(z - zs[k0]) < 1e-12]
    d += "*BOUNDARY\n" + "".join(f"{k},3,3\n" for k in ends)
    d += f"{ends[0]},1,2\n{ends[-1]},2,2\n"
    return _bending(work, f"bricks_{L:g}", d)


def sandwich_beam(L, layup, mats):
    """First bending frequency of the simply supported narrow strip, with and without core shear."""
    H = sum(t for t, _ in layup)
    z = -H / 2
    A = B = D = m = 0.0
    for t, name in layup:
        mat = mats[name]
        A += mat.E2 * t
        B += mat.E2 * ((z + t) ** 2 - z ** 2) / 2
        D += mat.E2 * ((z + t) ** 3 - z ** 3) / 3
        m += mat.rho * t
        z += t
    D -= B * B / A
    i = int(np.argmax([t for t, _ in layup]))     # the core
    c, core = layup[i][0], mats[layup[i][1]]
    d = c + 0.5 * (layup[i - 1][0] + layup[i + 1][0])
    S_ = core.G23 * d * d / c
    fb = (np.pi / L) ** 2 * np.sqrt(D / m) / (2 * np.pi)
    return fb, fb / np.sqrt(1 + D / S_ * (np.pi / L) ** 2)


def main(out=OUT) -> dict:
    layup, all_mats = _layup()
    mats = {name: all_mats[name] for _, name in layup}
    rows = []
    with tempfile.TemporaryDirectory() as work:
        for L in SPANS:
            fs, fb = shell(work, L, layup, mats)[0], bricks(work, L, layup, mats)[0]
            pure, sand = sandwich_beam(L, layup, mats)
            rows.append({"span": L, "shell_s8r": fs, "bricks_c3d20": fb, "beam_no_core_shear": pure,
                         "sandwich_beam": sand, "shell_over_bricks": fs / fb - 1.0})
    res = {"rows": rows}  # the layup itself is blade input data and stays private
    with open(out, "w", encoding="utf-8") as f:
        json.dump(res, f, indent=1)
    return res


if __name__ == "__main__":
    print(json.dumps(main(), indent=1))

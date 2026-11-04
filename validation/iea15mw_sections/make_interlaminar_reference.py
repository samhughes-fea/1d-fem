"""Interlaminar stresses of the 70 % span section from a CalculiX ``C3D20R`` prism.

The 2D mesh of ``examples/blade_section_2d`` (two rows per ply) is extruded into a prism
two chords long (``prism.py``) of fully integrated ``C3D20`` bricks, clamped at the root
and loaded at the tip so that a slice of bricks at mid-length carries the example's section
forces: the flapwise shear ``V_z = 200 kN`` with ``M_y = -4 MN m`` and the torque
``T = 100 kN m``. CalculiX prints the stresses at the integration points of that slice; the
middle layer of each brick's ``3 x 3 x 3`` points gives the stress at the slice's centre at
the ``3 x 3`` in-plane Gauss points of its 2D element. (Reduced-integration ``C3D20R``
bricks, which the stiffness runs use, give interlaminar shear that alternates in sign from
one thin ply element to the next here: their bricks are up to 500 times longer than a ply
row is thick.)

The output, ``reference/prism_20_interlaminar.json``, holds per brick the in-plane Gauss
point stresses in section axes (``[xx, yy, zz, yz, xz, xy]``). The interface stresses are
taken from them by the same fit as ``SolidSectionAnalysis.interface_stress`` uses
(``examples/blade_section_2d/run.py``). Needs ``ccx`` and implicit-geometry-kernel, and
a few minutes and 8 GB on four cores::

    python validation/iea15mw_sections/make_interlaminar_reference.py
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
for p in (REPO, HERE, os.path.join(REPO, "examples", "iea15mw"), os.path.join(REPO, "examples", "blade_section_2d")):
    if p not in sys.path:
        sys.path.insert(0, p)

import prism  # noqa: E402

LENGTH_CHORDS, N_ALONG = 2.0, 12
OUT = os.path.join(HERE, "reference", "prism_20_interlaminar.json")
ELEMENT = "C3D20"
#: CalculiX's integration point order in a C3D20 (xi fastest, then eta, then zeta, each at
#: -a, 0, +a): (xi, eta) of the nine points of each zeta layer, in units of a = sqrt(3/5)
IP_XI_ETA = [(i, j) for j in (-1, 0, 1) for i in (-1, 0, 1)]
N_IP = 27


def slice_forces_tip(length: float, x_c: float, forces) -> list:
    """Tip loads ``[N, V_y, V_z, T, M_y, M_z]`` that give ``forces`` (FORCES order) at ``x_c``:
    the shear's moment grows towards the root (``prism.energy_flexibility``)."""
    N, My, Mz, Vy, Vz, T = forces
    arm = length - x_c
    return [N, Vy, Vz, T, My + arm * Vz, Mz - arm * Vy]


def read_stresses(path: str, n: int) -> np.ndarray:
    """``(n N_IP, 8)`` rows ``[element, point, 6 stresses]``, stresses ``[xx, yy, zz, xy, xz, yz]`` from a ``.dat``,
    elements in the order printed."""
    rows = []
    grab = False
    with open(path, encoding="utf-8", errors="replace") as f:
        for line in f:
            s = line.split()
            if line.strip().startswith("stresses"):
                grab = True
                continue
            if grab and len(s) == 8:
                rows.append([float(v) for v in s[:8]])
    a = np.array(rows)
    if a.shape[0] != N_IP * n:
        raise RuntimeError(f"expected {N_IP * n} integration points, read {a.shape[0]}")
    return a


def run(work_dir: str, n_along: int = N_ALONG, length_chords: float = LENGTH_CHORDS) -> dict:
    import run as example                                  # examples/blade_section_2d/run.py
    import stations as S

    _, mesh, _, _, _ = example.build()
    m = S.station_model(example.INDEX)
    length = length_chords * m.chord
    j = n_along // 2
    dx = length / n_along
    x_c = (j + 0.5) * dx
    N, Vy, Vz, T, My, Mz = slice_forces_tip(length, x_c, example.section_forces())
    os.makedirs(work_dir, exist_ok=True)
    n2d = len(mesh.elements)
    slice_ids = [1 + e * n_along + j for e in range(n2d)]
    sl = ["*ELSET, ELSET=SLICE"] + [", ".join(map(str, slice_ids[a:a + 12])) + "," for a in range(0, n2d, 12)]
    book = {}

    def steps(ref, rot):
        return sl + ["*STEP", "*STATIC", "*CLOAD",
                     f"{ref}, 1, {N:.10e}", f"{ref}, 2, {Vy:.10e}", f"{ref}, 3, {Vz:.10e}",
                     f"{rot}, 1, {T:.10e}", f"{rot}, 2, {My:.10e}", f"{rot}, 3, {Mz:.10e}",
                     "*EL PRINT, ELSET=SLICE, GLOBAL=YES", "S", "*END STEP"]

    inp = os.path.join(work_dir, "prism.inp")
    # the reference ids are known only after the nodes are numbered: write once to learn them
    book = prism.write_deck(mesh, length, n_along, inp, (j, j + 1), steps=[], element=ELEMENT)
    prism.write_deck(mesh, length, n_along, inp, (j, j + 1), steps=steps(book["ref"], book["rot"]), element=ELEMENT)
    t0 = time.time()
    env = dict(os.environ, OMP_NUM_THREADS=str(os.cpu_count() or 1))
    with open(os.path.join(work_dir, "ccx.log"), "w", encoding="utf-8") as log:
        subprocess.run(["ccx", "-i", "prism"], cwd=work_dir, check=True, env=env, stdout=log,
                       stderr=subprocess.STDOUT)
    a = read_stresses(os.path.join(work_dir, "prism.dat"), n2d)
    order = {int(e): k for k, e in enumerate(slice_ids)}
    out = np.zeros((n2d, 9, 6))
    for row in range(0, a.shape[0], N_IP):
        e = order[int(a[row, 0])]
        s = a[row + 9:row + 18, 2:]                         # the middle layer: the slice's centre
        out[e] = s[:, [0, 1, 2, 5, 4, 3]]                   # -> [xx, yy, zz, yz, xz, xy]
    return {"index": example.INDEX, "length": length, "n_along": n_along, "length_chords": length_chords, "slice": j, "x_c": x_c,
            "section_forces": list(map(float, example.section_forces())), "ip_xi_eta": IP_XI_ETA,
            "ip_scale": float(np.sqrt(0.6)), "element": ELEMENT,
            "gauss_stress": [[[float(f"{v:.5g}") for v in pt] for pt in el] for el in out], "n_elements_2d": n2d, "time": time.time() - t0}


def main(n_along: int = N_ALONG, length_chords: float = LENGTH_CHORDS, out: str = OUT):
    r = run(os.path.join(HERE, "output", "prism_20_interlaminar"), n_along, length_chords)
    os.makedirs(os.path.dirname(out), exist_ok=True)
    with open(out, "w", encoding="utf-8") as f:
        json.dump(r, f)
    print("done in", round(r["time"]), "s")


if __name__ == "__main__":
    a = sys.argv[1:]
    if a:                                   # n_along, length in chords, output path (convergence runs)
        main(int(a[0]), float(a[1]), a[2])
    else:
        main()

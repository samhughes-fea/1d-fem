"""CalculiX reference for GBT on 2D sections: free vibration of a thin box cantilever.

The box (0.3 m x 0.12 m, 2 mm aluminium walls, 1.5 m long) is the 2D section mesh of
``tests/test_gbt_2d_sections.py`` extruded into ``C3D20R`` bricks (``iea15mw_sections/prism.py``),
clamped at the root (all translations) and free at the tip. Its lowest modes are local
flange modes, which the 2D section's distortional modes must reproduce. Writes
``reference.json`` here. Needs ``ccx`` on the PATH (about a minute)::

    python post_processing/validation_visualisers/gbt_2d_box/make_reference.py
"""

from __future__ import annotations

import json
import os
import re
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
for p in (REPO, os.path.join(REPO, "post_processing", "validation_visualisers", "iea15mw_sections")):
    if p not in sys.path:
        sys.path.insert(0, p)

N_ALONG = 60
N_MODES = 12


def box_mesh():
    from tests.test_gbt_2d_sections import box_mesh as mesh

    return mesh()


def main(work_dir: str = "") -> dict:
    import prism

    mesh = box_mesh()
    work_dir = work_dir or tempfile.mkdtemp(prefix="gbt_box_")
    deck = os.path.join(work_dir, "box.inp")
    prism.write_deck(mesh, 1.5, N_ALONG, deck, (0, 1))
    text = open(deck, encoding="utf-8").read().split("*RIGID BODY")[0]
    for mid, m in mesh.materials.items():
        text = text.replace(f"*MATERIAL, NAME=M{mid}\n", f"*MATERIAL, NAME=M{mid}\n*DENSITY\n{m.rho!r}\n")
    text += f"*BOUNDARY\nROOT, 1, 3, 0.0\n*STEP\n*FREQUENCY\n{N_MODES}\n*END STEP\n"
    with open(deck, "w", encoding="utf-8") as f:
        f.write(text)
    subprocess.run(["ccx", "box"], cwd=work_dir, check=True, stdout=subprocess.DEVNULL)
    dat = open(os.path.join(work_dir, "box.dat"), encoding="utf-8").read()
    block = dat.split("E I G E N V A L U E   O U T P U T")[1]
    freqs = [float(ln.split()[3]) for ln in block.splitlines()
             if re.match(r"\s+\d+\s+[-\d.E+]+\s+[-\d.E+]+\s+[-\d.E+]+", ln)][:N_MODES]
    out = {"element": "C3D20R", "n_along": N_ALONG, "frequencies_hz": freqs}
    with open(os.path.join(HERE, "reference.json"), "w", encoding="utf-8") as f:
        json.dump(out, f, indent=1)
    return out


if __name__ == "__main__":
    print(main()["frequencies_hz"])

"""CalculiX ``C3D20R`` prisms of the 2D IEA 15 MW sections (``prism.py``).

Needs ``ccx`` (CalculiX) on the PATH and about an hour per station on four cores, so the
output, ``reference/prism_<index>.json``, is committed. To regenerate it::

    python post_processing/validation_visualisers/iea15mw_sections/make_prism_reference.py [index ...]

Each prism is two chords long, 40 bricks along; the strain energy of its middle 14/40 is
used (``prism.py``). ``--chords 4`` writes ``prism_<index>_4chords.json`` from a prism four
chords long (energy of bricks 16 to 24), which shows how much of the multi-cell shear and
torsion stiffness of the two-chord prism comes from its restrained ends.
"""

from __future__ import annotations

import json
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
for p in (REPO, HERE, os.path.join(REPO, "examples", "iea15mw")):
    if p not in sys.path:
        sys.path.insert(0, p)

import prism  # noqa: E402
import stations as S  # noqa: E402

LENGTH_CHORDS, N_ALONG, STATIONS = 2.0, 40, (13, 27)


def run(index: int, work_dir: str, length_chords: float = LENGTH_CHORDS, stations=STATIONS) -> dict:
    from pre_processing.section_integration.solid.from_geometry import section_mesh_from_cut

    m = S.station_model(index)
    mesh = section_mesh_from_cut(S.geometry(index).section_at(0.5), S.laminates(index))
    t0 = time.time()
    r = prism.section_stiffness(mesh, length_chords * m.chord, N_ALONG, work_dir, stations=stations)
    return {"index": index, "K": r["stiffness"].tolist(), "n_nodes": r["n_nodes"], "n_elements": r["n_elements"],
            "time": time.time() - t0, "length_chords": length_chords, "n_along": N_ALONG, "stations": list(stations)}


def main(indices=(0, 20), long: bool = False):
    os.makedirs(os.path.join(HERE, "reference"), exist_ok=True)
    for i in indices:
        name = f"prism_{i}_4chords" if long else f"prism_{i}"
        out = (run(i, os.path.join(HERE, "output", name), 4.0, (16, 24)) if long
               else run(i, os.path.join(HERE, "output", name)))
        with open(os.path.join(HERE, "reference", f"{name}.json"), "w", encoding="utf-8") as f:
            json.dump(out, f, indent=1)
        print(i, "done in", round(out["time"]), "s")


if __name__ == "__main__":
    args = [a for a in sys.argv[1:] if a not in ("--chords", "4")]
    main(tuple(int(a) for a in args) or (0, 20), long="--chords" in sys.argv[1:])

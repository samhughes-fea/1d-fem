"""composite-laminate-model shell segments of the idealised IEA 15 MW stations.

Needs composite-laminate-model and implicit-geometry-kernel, which 1d-fem does not
depend on, so the output, ``reference/shell_<index>.json``, is committed. To regenerate::

    python post_processing/validation_visualisers/iea15mw_sections/make_shell_reference.py [index ...]

Each station is lofted prismatic over 8 chords from the same IGK geometry the 1d-fem
models cut (``examples/iea15mw/stations.py``) and meshed with curved Quad8 FSDT shells
on the walls' mid-surfaces, 32 elements along the span. ``beam_section_properties``
loads it as a cantilever through a rigid tip and measures the middle fifth. Its axes
are x chordwise and y flapwise about the pitch axis, the same as the 1d-fem section's
(y, z).
"""

from __future__ import annotations

import json
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
sys.path.insert(0, os.path.join(REPO, "examples", "iea15mw"))
import stations as S  # noqa: E402

LENGTH_CHORDS, N_SPAN = 8.0, 32


def shell_cards():
    from composite_laminate_model.laminate_materials import MaterialCard
    from composite_laminate_model.laminate_materials_3d import ElasticProperties3D, StrengthProperties3D

    big = 1e12      # strengths are not used
    return {n: MaterialCard(n, ElasticProperties3D.transversely_isotropic(p["E1"], p["E2"], p["nu12"], p["G12"], S.NU23),
                            StrengthProperties3D(big, big, big, big, big), density=p["rho"])
            for n, p in S.DATA["materials"].items()}


def shell_laminates(model):
    """Laminates per wall kind, inside first (the bottom face of an outward-normal shell)."""
    from composite_laminate_model import Laminate, Ply

    cards = shell_cards()
    return {k: Laminate([Ply(cards[n], t, 0.0) for n, t in v]) for k, v in model.laminae.items()}


def run(index: int) -> dict:
    from composite_laminate_model import BladeBuilder, FSDTModel
    from composite_laminate_model.laminate_blade import beam_section_properties

    m = S.station_model(index)
    geo = S.geometry(index, length=LENGTH_CHORDS * m.chord)
    t0 = time.time()
    blade = BladeBuilder.from_geometry(geo, n_span=N_SPAN, element_size=0.0315 * m.chord).build(shell_laminates(m), FSDTModel())
    p = beam_section_properties(blade)
    return {"index": index, "K4": p.stiffness.tolist(), "K6": p.stiffness_6x6.tolist(),
            "shear_centre": list(p.shear_centre), "tension_centre": list(p.tension_centre),
            "mass": p.mass_per_length, "time": time.time() - t0,
            "n_elements": len(blade.mesh.element_ids), "length_chords": LENGTH_CHORDS, "n_span": N_SPAN}


def main(indices=(0, 6, 20)):
    os.makedirs(os.path.join(HERE, "reference"), exist_ok=True)
    for i in indices:
        out = run(i)
        with open(os.path.join(HERE, "reference", f"shell_{i}.json"), "w", encoding="utf-8") as f:
            json.dump(out, f, indent=1)
        print(i, "done in", round(out["time"], 1), "s")


if __name__ == "__main__":
    main(tuple(int(a) for a in sys.argv[1:]) or (0, 6, 20))

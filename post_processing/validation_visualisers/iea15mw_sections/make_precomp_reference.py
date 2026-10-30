"""PreComp on the idealised IEA 15 MW sections (the layups every 1d-fem model here uses).

Needs WISDEM (``pip install wisdem``), which 1d-fem does not depend on, so the output,
``reference/precomp_idealised.json``, is committed. To regenerate it::

    python post_processing/validation_visualisers/iea15mw_sections/make_precomp_reference.py

Each section is the idealisation of ``examples/iea15mw/stations.py``: the skin all
round, the cap's carbon inside the skin between the webs, and the webs. PreComp lists
laminae from the outer surface inward and puts the cap sector boundaries on the web
positions. Its axes are x flapwise and y chordwise about the pitch axis.
"""

from __future__ import annotations

import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
sys.path.insert(0, os.path.join(REPO, "examples", "iea15mw"))
import stations as S  # noqa: E402


def precomp_inputs(index: int):
    m = S.station_model(index)
    names = list(S.DATA["materials"])
    idx = {n: k + 1 for k, n in enumerate(names)}          # 1-based for PreComp

    def sector(laminae_outer_first):
        return [(idx[n], t) for n, t in laminae_outer_first]

    skin = sector(list(reversed(m.laminae["skin"])))
    if m.web_positions:
        w0, w1 = m.web_positions
        cap = skin + sector(list(reversed(m.laminae["cap"])))
        loc, sectors = [0.0, w0, w1, 1.0], [skin, cap, skin]
    else:
        loc, sectors = [0.0, 1.0], [skin]

    def fmt(secs):
        n_lam = [len(s) for s in secs]
        flat = [x for s in secs for x in s]
        return (loc, n_lam, [1] * len(flat), [t for _, t in flat], [0.0] * len(flat), [k for k, _ in flat])

    pts = S.airfoil_points(index)                       # upper from the TE, then lower
    n_up = int(np.argmin(pts[:, 0]))
    up, lo = pts[:n_up + 1][::-1], pts[n_up:]
    # PreComp node order: LE round the upper surface to the TE and back along the lower one
    xnode = np.concatenate([up[:, 0], lo[1:, 0][::-1]])
    ynode = np.concatenate([up[:, 1], lo[1:, 1][::-1]])
    if np.hypot(*(up[-1] - lo[-1])) < 1e-12:
        xnode, ynode = xnode[:-1], ynode[:-1]
    if m.web_positions:
        web = sector(list(reversed(m.laminae["web"])))
        W = (list(m.web_positions), [len(web)] * 2, [1] * (2 * len(web)), [t for _, t in web] * 2,
             [0.0] * (2 * len(web)), [k for k, _ in web] * 2)
        nweb = 2
    else:
        nweb, W = 0, ([0], [0], [0], [0], [0], [0])
    mats = [S.DATA["materials"][n] for n in names]
    return dict(chord=m.chord, le_loc=m.pitch_axis, xnode=xnode, ynode=ynode, mats=mats,
                U=fmt(sectors), L=fmt(sectors), nweb=nweb, W=W)


def run(index: int) -> dict:
    from wisdem.precomp import properties

    p = precomp_inputs(index)
    m = p["mats"]
    out = properties.properties(
        p["chord"], 0.0, 0.0, p["le_loc"], p["xnode"], p["ynode"],
        [x["E1"] for x in m], [x["E2"] for x in m], [x["G12"] for x in m], [x["nu12"] for x in m],
        [x["rho"] for x in m], *p["U"], *p["L"], p["nweb"], *p["W"])
    (eif, eil, gj, ea, eifl, sf, sl, sft, slt, sat, x_sc, y_sc, x_tc, y_tc, mass, area, iflap, ilag, tw,
     x_cg, y_cg) = (float(np.squeeze(v)) for v in out)
    return {"EA": ea, "EI_flap": eif, "EI_edge": eil, "EI_flap_edge": eifl, "GJ": gj,
            "x_sc": x_sc, "y_sc": y_sc, "x_tc": x_tc, "y_tc": y_tc, "mass": mass, "area": area}


def main():
    res = {str(i): run(i) for i in (0, 6, 20)}
    os.makedirs(os.path.join(HERE, "reference"), exist_ok=True)
    with open(os.path.join(HERE, "reference", "precomp_idealised.json"), "w", encoding="utf-8") as f:
        json.dump(res, f, indent=1)
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()

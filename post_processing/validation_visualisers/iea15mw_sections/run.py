"""IEA 15 MW blade-section validation: 1d-fem's 2D and thin-walled sections against a
composite-laminate-model shell segment, PreComp and a CalculiX C3D20R prism.

Run from the repository root (needs implicit-geometry-kernel, ``1d-fem[geometry]``)::

    python post_processing/validation_visualisers/iea15mw_sections/run.py

Three stations of the IEA Wind 15 MW reference blade (root, maximum chord, 70 % span)
are cut by implicit-geometry-kernel with the idealisation of
``examples/iea15mw/stations.py``. Every model sees the same geometry and layups:

- **thin-walled:** ``thin_walled_from_cut`` (``hoop="free"``), computed here;
- **2D:** ``section_mesh_from_cut`` (layered ``Q8``), computed here;
- **shell:** composite-laminate-model shell segments, ``reference/shell_*.json``
  (``make_shell_reference.py``);
- **PreComp:** WISDEM's PreComp on the same idealised layups,
  ``reference/precomp_idealised.json`` (``make_precomp_reference.py``);
- **prism:** the 2D mesh extruded into CalculiX ``C3D20R`` bricks,
  ``reference/prism_*.json`` (``make_prism_reference.py``).

Writes ``summary.json``, ``report.md`` and the figures next to this file.
"""

from __future__ import annotations

import json
import os
import sys
import time
from typing import Dict, Optional

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
for p in (REPO, os.path.join(REPO, "examples", "iea15mw")):
    if p not in sys.path:
        sys.path.insert(0, p)
import stations as S  # noqa: E402

STATIONS = (0, 6, 20)
QUANTITIES = ("EA", "EI_flap", "EI_edge", "GJ")
LABELS = {"EA": "EA", "EI_flap": "EI flap", "EI_edge": "EI edge", "GJ": "GJ"}
TARGETS = {"EA": 0.01, "EI_flap": 0.01, "EI_edge": 0.01, "GJ": 0.03}
CENTRE_TARGET = 0.01          # of the chord
MODELS = ("solid_2d", "shell", "precomp", "prism")
NAMES = {"thin_walled": "1d-fem thin-walled", "solid_2d": "1d-fem 2D", "shell": "CLM shell",
         "precomp": "PreComp", "prism": "CalculiX prism"}
# reference palette slots 1-4 (light mode)
COLOURS = {"solid_2d": "#2a78d6", "shell": "#eb6834", "precomp": "#1baf7a", "prism": "#eda100"}
MATERIAL_COLOURS = {"CarbonUD": "#2a78d6", "glass_triax": "#1baf7a", "glass_biax": "#1baf7a",
                    "medium_density_foam": "#eda100", "Gelcoat": "#b8b6ad", "glass_uni": "#6250d6"}


# ---------------------------------------------------------------------------
# Comparable quantities
# ---------------------------------------------------------------------------
def summary_4x4(K4: np.ndarray, tension_centre, shear_centre, chord: float) -> Dict[str, object]:
    """EA, bending about the tension centre, GJ (shears free) from a 4 x 4 [N, M_flap, M_edge, T]."""
    K4 = np.asarray(K4, float)
    EA = K4[0, 0]
    B = K4[1:3, 1:3] - np.outer(K4[1:3, 0], K4[0, 1:3]) / EA
    return {"EA": float(EA), "EI_flap": float(B[0, 0]), "EI_edge": float(B[1, 1]), "EI_flap_edge": float(B[0, 1]),
            "GJ": float(1.0 / np.linalg.inv(K4)[3, 3]),
            "tc": [float(v) / chord for v in tension_centre],
            "sc": None if shear_centre is None else [float(v) / chord for v in shear_centre]}


def summary_1dfem(K6: np.ndarray, tension_centre, shear_centre, chord: float) -> Dict[str, object]:
    """1d-fem 6 x 6 (``FORCES`` order N, M_y, M_z, V_y, V_z, T; y chordwise, z flapwise)."""
    C = np.linalg.inv(np.asarray(K6, float))
    keep = [0, 1, 2, 5]
    out = summary_4x4(np.linalg.inv(C[np.ix_(keep, keep)]), tension_centre, shear_centre, chord)
    out["GA_edge"], out["GA_flap"] = float(K6[3][3]), float(K6[4][4])
    return out


def summary_shell(ref: dict, chord: float) -> Dict[str, object]:
    """composite-laminate-model 4 x 4 [N, M_x, M_y, M_z] (x chordwise, y flapwise)."""
    out = summary_4x4(ref["K4"], ref["tension_centre"], ref["shear_centre"], chord)
    out["GA_edge"], out["GA_flap"] = float(ref["K6"][0][0]), float(ref["K6"][1][1])
    return out


def summary_precomp(p: dict, chord: float) -> Dict[str, object]:
    """PreComp: x flapwise, y chordwise (aft positive) from the pitch axis, bending about the
    tension centre. PreComp's shear centre is not used: it sits on its tension centre."""
    return {"EA": p["EA"], "EI_flap": p["EI_flap"], "EI_edge": p["EI_edge"], "EI_flap_edge": p["EI_flap_edge"],
            "GJ": p["GJ"], "tc": [p["y_tc"] / chord, p["x_tc"] / chord], "sc": None}


def _load(name: str) -> Optional[dict]:
    path = os.path.join(HERE, "reference", name)
    if not os.path.exists(path):
        return None
    with open(path, encoding="utf-8") as f:
        return json.load(f)


# ---------------------------------------------------------------------------
# The 1d-fem models
# ---------------------------------------------------------------------------
def models(index: int, rows_per_ply: int = 1, with_2d: bool = True):
    """Thin-walled and 2D analyses of one station (2D is ``None`` where the mesher fails)."""
    from pre_processing.section_integration.solid.from_geometry import section_mesh_from_cut
    from pre_processing.section_integration.solid.section_analysis import SolidSectionAnalysis
    from pre_processing.section_integration.thin_walled.from_geometry import thin_walled_from_cut

    cut = S.geometry(index).section_at(0.5)
    lams = S.laminates(index)
    t0 = time.time()
    gs = thin_walled_from_cut(cut, lams)
    tw = gs.analysis()
    _ = tw.stiffness
    out = {"cut": cut, "thin_walled_geometry": gs, "thin_walled": tw, "times": {"thin_walled": time.time() - t0}}
    out["solid_2d"] = out["mesh"] = None
    out["solid_2d_error"] = None
    if with_2d:
        t0 = time.time()
        try:
            mesh = section_mesh_from_cut(cut, lams, rows_per_ply=rows_per_ply)
            sa = SolidSectionAnalysis(mesh)
            _ = sa.stiffness
            out["mesh"], out["solid_2d"] = mesh, sa
        except ValueError as exc:          # the layered mesher cannot fit this section
            out["solid_2d_error"] = str(exc)
        out["times"]["solid_2d"] = time.time() - t0
    return out


def carbon_matched_2d(index: int, ratio: float):
    """The 2D section with the cap material scaled by ``ratio`` (a diagnostic: equal carbon
    stiffness to the thin-walled cap, whose carbon is narrower)."""
    from pre_processing.section_integration.solid.from_geometry import section_mesh_from_cut
    from pre_processing.section_integration.solid.section_analysis import SolidSectionAnalysis
    from pre_processing.section_integration.thin_walled import Ply, PlyMaterial, WallLaminate

    m = S.station_model(index)
    c = S.ply_materials()["CarbonUD"]
    cs = PlyMaterial(E1=ratio * c.E1, E2=ratio * c.E2, G12=ratio * c.G12, nu12=c.nu12, G13=ratio * c.G13,
                     G23=ratio * c.G23, rho=c.rho)
    lams = dict(S.laminates(index))
    lams["cap"] = WallLaminate([Ply(cs, 0.0, t) for _, t in m.laminae["cap"]])
    return SolidSectionAnalysis(section_mesh_from_cut(S.geometry(index).section_at(0.5), lams))


def carbon_areas(index: int, built) -> Dict[str, float]:
    """Carbon area of the 2D mesh and of the thin-walled cap walls."""
    from pre_processing.section_integration.solid.elements import shape

    mesh, gs = built["mesh"], built["thin_walled_geometry"]
    t_c = sum(t for _, t in S.station_model(index).laminae["cap"])
    a2 = 0.0
    g, wg = np.polynomial.legendre.leggauss(3)
    for e, (conn, kind) in enumerate(zip(mesh.elements, mesh.kinds)):
        if not mesh.labels[e].startswith("spar_cap"):
            continue
        X = mesh.nodes[conn]
        a2 += sum(wx * we * np.linalg.det(shape(kind, p, q)[1] @ X) for p, wx in zip(g, wg) for q, we in zip(g, wg))
    sec = gs.section
    at = 0.0
    for w in sec.walls:
        if gs.wall_names[w.wall_id].startswith("spar_cap"):
            a, b = (np.asarray(p, float) for p in sec.endpoints(w))
            at += np.linalg.norm(b - a) * t_c
    return {"solid_2d": float(a2), "thin_walled": float(at)}


# ---------------------------------------------------------------------------
# Comparison
# ---------------------------------------------------------------------------
def compare(stations=STATIONS, rows_per_ply: int = 1) -> dict:
    precomp = _load("precomp_idealised.json") or {}
    rows = {}
    built_all = {}
    for i in stations:
        m = S.station_model(i)
        b = models(i, rows_per_ply=rows_per_ply)
        built_all[i] = b
        tw = b["thin_walled"]
        res = {"name": m.name, "chord": m.chord, "times": b["times"],
               "thin_walled": summary_1dfem(tw.stiffness, tw.tension_centre, tw.shear_centre, m.chord)}
        if b["solid_2d"] is not None:
            sa = b["solid_2d"]
            res["solid_2d"] = summary_1dfem(sa.stiffness, sa.tension_centre, sa.shear_centre, m.chord)
            res["solid_2d"]["n_elements"] = len(b["mesh"].elements)
        else:
            res["solid_2d_error"] = b["solid_2d_error"]
        shell = _load(f"shell_{i}.json")
        if shell:
            res["shell"] = summary_shell(shell, m.chord)
        if str(i) in precomp:
            res["precomp"] = summary_precomp(precomp[str(i)], m.chord)
        prism = _load(f"prism_{i}.json")
        if prism:
            res["prism"] = summary_1dfem(prism["K"], (np.nan, np.nan), None, m.chord)
            res["prism"]["tc"] = None
            res["prism_vs_2d_6x6"] = prism_6x6_check(prism, b["solid_2d"])
        res["ratios"] = {k: {q: res[k][q] / res["thin_walled"][q] - 1.0 for q in QUANTITIES}
                         for k in MODELS if k in res}
        if "prism" in res and "solid_2d" in res:
            res["prism_vs_2d"] = {q: res["prism"][q] / res["solid_2d"][q] - 1.0 for q in QUANTITIES}
        res["centre_offsets"] = {}
        for k in MODELS:
            if k not in res:
                continue
            for c in ("tc", "sc"):
                a, r = res[k].get(c), res["thin_walled"][c]
                if a is not None and r is not None:
                    res["centre_offsets"][f"{k}_{c}"] = float(np.hypot(a[0] - r[0], a[1] - r[1]))
        if b["solid_2d"] is not None and m.web_positions:
            areas = carbon_areas(i, b)
            ratio = areas["thin_walled"] / areas["solid_2d"]
            sm = carbon_matched_2d(i, ratio)
            cm = summary_1dfem(sm.stiffness, sm.tension_centre, sm.shear_centre, m.chord)
            res["carbon"] = {"areas": areas, "ratio": ratio,
                             "matched_ratios": {q: cm[q] / res["thin_walled"][q] - 1.0 for q in QUANTITIES}}
        rows[i] = res
    return {"stations": rows, "built": built_all}


def prism_6x6_check(prism: dict, sa) -> Optional[dict]:
    """Prism 6 x 6 against the 2D 6 x 6: each term relative to the geometric mean of its diagonals."""
    if sa is None:
        return None
    Kp, K2 = np.asarray(prism["K"]), np.asarray(sa.stiffness)
    d = np.sqrt(np.outer(np.diag(K2), np.diag(K2)))
    rel = (Kp - K2) / d
    return {"max_normalised_difference": float(np.abs(rel).max()), "normalised_difference": rel.tolist(),
            "diagonal_ratio": (np.diag(Kp) / np.diag(K2) - 1.0).tolist()}


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------
def _style(ax):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.spines["left"].set_color("#8b8a83")
    ax.spines["bottom"].set_color("#8b8a83")
    ax.tick_params(colors="#52514e", labelsize=8)


def figure_sections(result: dict, path: str):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.collections import PolyCollection

    st = list(result["stations"])
    fig, axes = plt.subplots(len(st), 1, figsize=(9, 2.6 * len(st)))
    for ax, i in zip(np.atleast_1d(axes), st):
        b = result["built"][i]
        m = S.station_model(i)
        if b["mesh"] is not None:
            mesh = b["mesh"]
            mats = {mid: name for mid, name in _material_names(mesh).items()}
            polys, cols = [], []
            for e, conn in enumerate(mesh.elements):
                polys.append(mesh.nodes[conn[:4]])
                cols.append(MATERIAL_COLOURS.get(mats.get(mesh.material_ids[e], ""), "#b8b6ad"))
            ax.add_collection(PolyCollection(polys, facecolors=cols, edgecolors="#fcfcfb", linewidths=0.15))
            ax.autoscale_view()
            note = f"2D mesh, {len(mesh.elements)} Q8 elements"
        else:
            note = "thin-walled mid-lines (the 2D mesher cannot fit this section)"
        gs = b["thin_walled_geometry"]
        for w in gs.section.walls:
            a, c = (np.asarray(p, float) for p in gs.section.endpoints(w))
            ax.plot([a[0], c[0]], [a[1], c[1]], color="#0b0b0b", lw=0.6)
        r = result["stations"][i]
        for key, mk in (("tc", "+"), ("sc", "x")):
            v = r["thin_walled"][key]
            ax.plot(v[0] * m.chord, v[1] * m.chord, mk, color="#e34948", ms=8, mew=1.5)
        ax.set_aspect("equal")
        ax.set_title(f"{m.name}, r = {S.station(i)['r']:.1f} m, chord {m.chord:.2f} m: {note}",
                     fontsize=9, color="#0b0b0b", loc="left")
        _style(ax)
    from matplotlib.patches import Patch
    from matplotlib.lines import Line2D
    handles = [Patch(color=MATERIAL_COLOURS[k], label=lab) for k, lab in
               (("CarbonUD", "carbon UD"), ("glass_triax", "glass triax / biax"), ("medium_density_foam", "foam"),
                ("Gelcoat", "gelcoat"))]
    handles += [Line2D([], [], color="#e34948", marker="+", ls="", label="tension centre"),
                Line2D([], [], color="#e34948", marker="x", ls="", label="shear centre"),
                Line2D([], [], color="#0b0b0b", lw=0.6, label="thin-walled mid-lines")]
    fig.legend(handles=handles, loc="lower center", ncol=4, fontsize=8, frameon=False)
    fig.tight_layout(rect=(0, 0.06, 1, 1))
    fig.savefig(path, dpi=150)
    plt.close(fig)


def _material_names(mesh) -> Dict[int, str]:
    """Material id of the 2D mesh to the IEA material name (matched by E1 and E2)."""
    out = {}
    for mid, mat in mesh.materials.items():
        e1 = getattr(mat, "E1", None)
        for name, p in S.DATA["materials"].items():
            if e1 is not None and abs(e1 - p["E1"]) < 1e-6 * p["E1"] and abs(mat.E2 - p["E2"]) < 1e-6 * p["E2"]:
                out[mid] = name
    return out


def figure_ratios(result: dict, path: str):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    st = list(result["stations"])
    fig, axes = plt.subplots(1, len(st), figsize=(10, 3.6), sharey=False)
    for ax, i in zip(np.atleast_1d(axes), st):
        r = result["stations"][i]
        x = np.arange(len(QUANTITIES))
        present = [k for k in MODELS if k in r["ratios"]]
        w = 0.8 / max(1, len(present))
        for j, k in enumerate(present):
            v = [100 * r["ratios"][k][q] for q in QUANTITIES]
            xs = x - 0.4 + w * (j + 0.5)
            ax.vlines(xs, 0, v, color=COLOURS[k], lw=2)
            ax.plot(xs, v, "o", color=COLOURS[k], ms=5, label=NAMES[k], mec="#fcfcfb", mew=1)
        for q_i, q in enumerate(QUANTITIES):
            t = 100 * TARGETS[q]
            ax.fill_between([q_i - 0.45, q_i + 0.45], -t, t, color="#e9e8e3", zorder=0, lw=0)
        ax.axhline(0, color="#8b8a83", lw=0.8)
        ax.set_xticks(x, [LABELS[q] for q in QUANTITIES])
        ax.set_title(f"{r['name']}", fontsize=9, loc="left", color="#0b0b0b")
        ax.set_ylabel("difference from thin-walled, %" if i == st[0] else "", fontsize=8, color="#52514e")
        _style(ax)
    h, lab = np.atleast_1d(axes)[0].get_legend_handles_labels()
    for ax in np.atleast_1d(axes)[1:]:
        for hh, ll in zip(*ax.get_legend_handles_labels()):
            if ll not in lab:
                h.append(hh)
                lab.append(ll)
    fig.legend(h, lab, loc="lower center", ncol=4, fontsize=8, frameon=False)
    fig.text(0.01, 0.01, "Grey bands: the plan's targets (1 % for EA and EI, 3 % for GJ).", fontsize=7, color="#52514e")
    fig.tight_layout(rect=(0, 0.1, 1, 1))
    fig.savefig(path, dpi=150)
    plt.close(fig)


# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------
def _pct(v):
    return "" if v is None else f"{100 * v:+.2f} %"


def report_md(result: dict) -> str:
    lines = ["# IEA 15 MW section validation", "",
             "Generated by `run.py`. Differences are relative to the 1d-fem thin-walled section "
             "(`hoop=\"free\"`). Bending is about the tension centre with the shear forces free; centre "
             "offsets are distances in the section plane as a fraction of the chord.", ""]
    for i, r in result["stations"].items():
        lines += [f"## {r['name']} (station {i}, chord {r['chord']:.2f} m)", "",
                  "| Model | " + " | ".join(LABELS[q] for q in QUANTITIES) + " | TC offset | SC offset |",
                  "|---|" + "---|" * (len(QUANTITIES) + 2)]
        tw = r["thin_walled"]
        lines.append("| thin-walled (value) | " + " | ".join(f"{tw[q]:.4g}" for q in QUANTITIES) + " | | |")
        for k in MODELS:
            if k not in r["ratios"]:
                continue
            off = r["centre_offsets"]
            tc, sc = off.get(f"{k}_tc"), off.get(f"{k}_sc")
            lines.append(f"| {NAMES[k]} | " + " | ".join(_pct(r["ratios"][k][q]) for q in QUANTITIES) +
                         f" | {'' if tc is None else f'{100 * tc:.2f} %'} | {'' if sc is None else f'{100 * sc:.2f} %'} |")
        if "solid_2d_error" in r:
            lines += ["", f"The 2D mesher could not mesh this section: `{r['solid_2d_error']}`."]
        if "carbon" in r:
            c = r["carbon"]
            lines += ["", f"Carbon area: 2D mesh {c['areas']['solid_2d']:.4f} m², thin-walled caps "
                          f"{c['areas']['thin_walled']:.4f} m² (2D is {100 * (1 / c['ratio'] - 1):.1f} % wider). With the "
                          "2D cap stiffness scaled to the same carbon: " +
                      ", ".join(f"{LABELS[q]} {_pct(v)}" for q, v in c["matched_ratios"].items()) + "."]
        if "prism_vs_2d" in r:
            lines += ["", "CalculiX prism against the 2D section (same mesh): " +
                      ", ".join(f"{LABELS[q]} {_pct(v)}" for q, v in r["prism_vs_2d"].items()) +
                      f"; largest 6 x 6 term difference {100 * r['prism_vs_2d_6x6']['max_normalised_difference']:.2f} % "
                      "of the geometric mean of its diagonals."]
        lines.append("")
    return "\n".join(lines)


def _jsonable(result: dict) -> dict:
    return {"stations": {str(i): r for i, r in result["stations"].items()}}


def main(stations=STATIONS, plot: bool = True, out_dir: str = HERE) -> dict:
    result = compare(stations)
    os.makedirs(out_dir, exist_ok=True)
    summary = _jsonable(result)
    with open(os.path.join(out_dir, "summary.json"), "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=1)
    with open(os.path.join(out_dir, "report.md"), "w", encoding="utf-8") as f:
        f.write(report_md(result))
    if plot:
        figure_sections(result, os.path.join(out_dir, "sections.png"))
        figure_ratios(result, os.path.join(out_dir, "stiffness_differences.png"))
    return summary


if __name__ == "__main__":
    s = main()
    print(open(os.path.join(HERE, "report.md"), encoding="utf-8").read())

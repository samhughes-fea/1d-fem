# examples/thin_walled/three_cell_blade/run.py
"""Example 2: three-cell blade section cantilever, from section analysis to ply stresses.

Run from the repository root::

    python examples/thin_walled/three_cell_blade/run.py

Writes the job inputs to ``jobs/`` and the 1d-fem results to ``results/`` (both
not tracked), and the figures and ``summary.json`` next to this file.
"""

from __future__ import annotations

import json
import os
import sys
from typing import Dict, List, Tuple

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))

from common import (  # noqa: E402
    CantileverSpec,
    Material,
    SectionSpec,
    bredt_multicell,
    by_ply,
    build_analysis,
    contour,
    plot_deflection,
    plot_interlaminar,
    plot_ply_stresses,
    plot_section,
    plot_section_fields,
    ply_materials,
    run,
)

GFRP = Material(E1=41e9, E2=10e9, G12=4e9, nu12=0.28, G13=4e9, G23=3.5e9, rho=1900.0)
CFRP = Material(E1=135e9, E2=10e9, G12=5e9, nu12=0.3, G13=5e9, G23=3.8e9, rho=1550.0)
L, N_EL = 10.0, 20        # span, elements
P, T = 10e3, 2e3          # tip flapwise force (z) and torque, N and N m
X_REF = 0.3               # reference axis at 30 % chord
CHORD = 1.0
WEBS = (0.25, 0.5)        # web chord positions; the spar caps run between them
# walls per skin segment: leading edge to front web, the cap, rear web to trailing edge
N_LE, N_CAP, N_TE = 7, 4, 8


def naca4(m: float, p: float, t: float):
    """Upper (``side = 1``) and lower (``side = -1``) surface ``(x, z)`` of a NACA 4-digit
    airfoil at chord parameter ``xc``, with the closed trailing edge (last coefficient -0.1036)."""

    def surface(xc: float, side: int) -> Tuple[float, float]:
        yt = 5 * t * (0.2969 * np.sqrt(xc) - 0.1260 * xc - 0.3516 * xc ** 2 + 0.2843 * xc ** 3 - 0.1036 * xc ** 4)
        if xc < p:
            yc, dyc = m / p ** 2 * (2 * p * xc - xc ** 2), 2 * m / p ** 2 * (p - xc)
        else:
            yc, dyc = m / (1 - p) ** 2 * (1 - 2 * p + 2 * p * xc - xc ** 2), 2 * m / (1 - p) ** 2 * (p - xc)
        th = np.arctan(dyc)
        return float(xc - side * yt * np.sin(th)), float(yc + side * yt * np.cos(th))

    return surface


def _surface_points(surface, side: int) -> List[Tuple[float, float]]:
    """Points of one skin from the leading edge to the trailing edge, with points at the webs."""
    from scipy.optimize import brentq

    def at_x(x: float) -> Tuple[float, float]:
        if x <= 0.0:
            return 0.0, 0.0
        if x >= 1.0:
            return 1.0, 0.0
        xc = brentq(lambda q: surface(q, side)[0] - x, 1e-9, 1.0 - 1e-9)
        return surface(xc, side)

    # cosine spacing toward the leading edge, even spacing on the cap and aft panel
    le = WEBS[0] * (1 - np.cos(np.linspace(0, np.pi / 2, N_LE + 1)))
    cap = np.linspace(WEBS[0], WEBS[1], N_CAP + 1)
    te = np.linspace(WEBS[1], 1.0, N_TE + 1)
    xs = np.concatenate([le, cap[1:], te[1:]])
    return [at_x(x) for x in xs]


GROUPS = ("LE upper skin", "upper cap", "TE upper skin", "TE lower skin", "lower cap", "LE lower skin",
          "front web", "rear web")


def section_layout():
    """Points, walls and the wall groups of the NACA 4418 three-cell section."""
    surf = naca4(0.04, 0.4, 0.18)
    up = _surface_points(surf, 1)
    lo = _surface_points(surf, -1)
    pts: Dict[int, Tuple[float, float]] = {}
    pid = 1
    up_id, lo_id = [], []
    for k, (x, z) in enumerate(up):
        pts[pid] = (x * CHORD - X_REF, z * CHORD)
        up_id.append(pid)
        pid += 1
    lo_id.append(up_id[0])
    for x, z in lo[1:-1]:
        pts[pid] = (x * CHORD - X_REF, z * CHORD)
        lo_id.append(pid)
        pid += 1
    lo_id.append(up_id[-1])
    i_w1, i_w2 = N_LE, N_LE + N_CAP       # index of the web points along each skin
    walls, groups = [], {g: [] for g in GROUPS}
    wid = 1

    def add(a, b, lam, g):
        nonlocal wid
        walls.append((wid, a, b, lam))
        groups[g].append(wid)
        wid += 1

    # contour order: leading edge, along the upper skin to the trailing edge, back along the lower skin
    for i in range(len(up_id) - 1):
        g = "LE upper skin" if i < i_w1 else ("upper cap" if i < i_w2 else "TE upper skin")
        add(up_id[i], up_id[i + 1], 2 if g == "upper cap" else 1, g)
    for i in range(len(lo_id) - 1, 0, -1):
        j = i - 1
        g = "LE lower skin" if j < i_w1 else ("lower cap" if j < i_w2 else "TE lower skin")
        add(lo_id[i], lo_id[i - 1], 2 if g == "lower cap" else 1, g)
    add(up_id[i_w1], lo_id[i_w1], 3, "front web")
    add(up_id[i_w2], lo_id[i_w2], 3, "rear web")
    return pts, walls, groups


def section(n_sub: int = 2) -> SectionSpec:
    skin = [(1, a, 0.5e-3) for a in (0, 45, -45, -45, 45, 0)]          # 3 mm triax
    cap = [(1, 45, 0.5e-3), (1, -45, 0.5e-3)] + [(2, 0, 0.5e-3)] * 16 + [(1, -45, 0.5e-3), (1, 45, 0.5e-3)]  # 10 mm
    web = [(1, a, 0.5e-3) for a in (45, -45, -45, 45)] * 2              # 4 mm biax
    pts, walls, _ = section_layout()
    return SectionSpec(points=pts, walls=walls, materials={1: GFRP, 2: CFRP},
                       laminates={1: skin, 2: cap, 3: web}, n_sub=n_sub)


GROUP_WALLS = section_layout()[2]
#: one wall per group for the tables: the middle wall of the group
WALL_NAMES = {ws[len(ws) // 2]: g for g, ws in GROUP_WALLS.items()}


def cases(n_el: int = N_EL, n_sub: int = 2):
    sec = section(n_sub)
    load = [0.0, 0.0, P, T, 0.0, 0.0]
    return {
        "Euler-Bernoulli": CantileverSpec(L, n_el, "LinearEulerBernoulliBeamElement3D", sec, load),
        "Timoshenko": CantileverSpec(L, n_el, "LinearTimoshenkoBeamElement3D", sec, load),
        "Benscoter": CantileverSpec(L, n_el, "LinearBenscoterBeamElement3D", sec, load, warping=1),
        "GBT": CantileverSpec(L, n_el, "LinearGBTBeamElement3D", sec, load, modes=2, shear_modes=4),
    }


def exact_tip(analysis, length: float, load):
    """Tip ``(v, w, twist)`` of the Saint-Venant beam by unit loads: int F_u^T C F dx."""
    C = analysis.flexibility
    xs, ws = np.polynomial.legendre.leggauss(8)
    x = 0.5 * length * (xs + 1)
    w = 0.5 * length * ws

    def forces(xx, Fy, Fz, Mx):  # [N, M_y, M_z, V_y, V_z, T] of tip loads at the reference axis
        return np.array([0.0, Fz * (xx - length), -Fy * (xx - length), Fy, Fz, Mx])

    F = [forces(xi, 0.0, load[2], load[3]) for xi in x]
    unit = [(1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)]
    return [float(sum(wi * forces(xi, *u) @ C @ Fi for xi, wi, Fi in zip(x, w, F))) for u in unit]


def _mid(an, wall_id: int) -> float:
    return 0.5 * an.section.length(an.section.wall(wall_id))


def main(n_el: int = N_EL, plot: bool = True, work_dir: str = HERE) -> dict:
    from onedfem.pre_processing.section_integration.thin_walled.stress_recovery import SaintVenantRecovery

    an = build_analysis(section())
    results = {name: run(spec, work_dir, name.replace(" ", "_").replace(",", "")) for name, spec in cases(n_el).items()}
    v_ex, w_ex, tw_ex = exact_tip(an, L, [0, 0, P, T])
    summary = {
        "section": {"airfoil": "NACA 4418, closed trailing edge", "chord": CHORD, "n_walls": len(an.section.walls),
                    "groups": GROUP_WALLS},
        "section_stiffness_6x6": an.stiffness.tolist(),
        "shear_centre": [float(c) for c in an.shear_centre],
        "exact_saint_venant_tip": {"v": v_ex, "w": w_ex, "twist": tw_ex},
        "tip": {k: {"v": float(r.U[-1, 1]), "w": float(r.U[-1, 2]), "twist": float(r.U[-1, 3])}
                for k, r in results.items()},
    }
    # torque shear flow: 1d-fem against multi-cell Bredt-Batho, in every wall
    q_b = bredt_multicell(an, T)
    rec = SaintVenantRecovery(an)
    q_fem = {w.wall_id: float(rec.at(w.wall_id, _mid(an, w.wall_id), [0, 0, 0, 0, 0, T]).resultants[2])
             for w in an.section.walls}
    summary["torque_shear_flow"] = {WALL_NAMES[w]: {"bredt": q_b[w], "1d-fem": q_fem[w]} for w in WALL_NAMES}
    summary["torque_shear_flow_largest_difference"] = max(abs(q_fem[w] / q_b[w] - 1.0) for w in q_b)
    # mid-span: upper cap centre and front web centre, exact section forces of the cantilever
    e_mid = n_el // 2
    x_mid = (e_mid + 0.5) * L / n_el
    F_mid = [0, P * (x_mid - L), 0, 0, P, T]
    cap_w = GROUP_WALLS["upper cap"][len(GROUP_WALLS["upper cap"]) // 2]
    web_w = GROUP_WALLS["front web"][0]
    cap_s, web_s = 0.0, _mid(an, web_w)    # the cap point is the start of its middle wall: mid-cap
    cap = rec.at(cap_w, cap_s, F_mid, points_per_ply=5)
    web = rec.at(web_w, web_s, F_mid, points_per_ply=5)
    cap_dense = rec.at(cap_w, cap_s, F_mid, points_per_ply=21)
    web_dense = rec.at(web_w, web_s, F_mid, points_per_ply=21)
    summary["mid_span_upper_cap"] = {
        "wall": cap_w, "N_xx": float(cap.resultants[0]), "N_xs": float(cap.resultants[2]),
        "max_sigma_11": float(max(abs(q.stress_ply[0]) for q in cap.plies)),
        "max_tau_xn": float(max(abs(q.tau_xn) for q in cap.plies)),
    }
    summary["mid_span_front_web"] = {
        "wall": web_w, "N_xs": float(web.resultants[2]), "Q_x": float(web.transverse_shear[0]),
        "max_tau_12": float(max(abs(q.stress_ply[2]) for q in web.plies)),
    }
    rows = [r for r in results["Timoshenko"].wall_rows
            if int(r["element_id"]) == e_mid and r["wall_id"] == str(cap_w) and abs(float(r["s"]) - cap_s) < 1e-9]
    summary["csv_vs_direct_N_xx"] = [float(rows[0]["N_xx"]), float(cap.resultants[0])] if rows else None
    # root element: N_xx around the whole contour, GBT against Timoshenko
    root = {k: [(int(r["wall_id"]), float(r["s"]), float(r["N_xx"])) for r in contour(results[k].wall_rows, 0)]
            for k in ("Timoshenko", "GBT")}
    summary["root_N_xx"] = root
    caps = set(GROUP_WALLS["upper cap"])
    peak = {k: max(abs(n) for w, s, n in v if w in caps) for k, v in root.items()}
    summary["root_upper_cap_peak_N_xx_gbt_over_timoshenko"] = peak["GBT"] / peak["Timoshenko"]
    with open(os.path.join(work_dir, "summary.json"), "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)
    if plot:
        _plots(an, rec, results, {k: c.element_type for k, c in cases(n_el).items()}, cap_dense, web_dense, root,
               q_b, n_el, (w_ex, tw_ex), cap_w, web_w)
    return summary


def _plots(an, rec, results, element_types, cap, web, root, q_b, n_el, exact, cap_w, web_w):
    """Figures. Fields around the section are the section's own element fields, sampled inside
    each wall element (``thin_walled_drawing``); curves along the beam are each element's
    shape-function field (``beam_curve``); through-thickness curves are drawn ply by ply; the
    root contour values from the job's CSV are shown as the discrete points they are."""
    from onedfem.post_processing.graphical_visualisers import figure_style as fs
    from onedfem.post_processing.graphical_visualisers.thin_walled_drawing import unrolled_field

    fs.use()
    import matplotlib.pyplot as plt

    mats = ply_materials(section(), {1: "GFRP", 2: "CFRP"})

    # 1. the section at true thickness, ply by ply, and a close-up of the front web joint
    p_web = an.section.endpoints(an.section.wall(web_w))[0]
    plot_section(an, mats, os.path.join(HERE, "section.png"), "NACA 4418 three-cell section, walls at true thickness",
                 f"{len(an.section.walls)} walls, {an.n_sub} cubic elements each (ticks: element ends); "
                 "skins 3 mm, caps 10 mm, webs 4 mm",
                 (p_web[0] - 0.045, p_web[1] - 0.03, 0.09, 0.045), "Front web at the upper cap")

    # 2. torque shear flow and the root flap state, drawn on the section
    F_root = [0.0, -P * L, 0.0, 0.0, P, 0.0]
    F_T = [0, 0, 0, 0, 0, T]
    offsets = {"LE upper skin": (-0.03, 0.045), "upper cap": (0.0, -0.03), "TE upper skin": (0.0, 0.045),
               "LE lower skin": (-0.03, -0.04), "lower cap": (0.0, 0.03), "TE lower skin": (0.0, -0.04),
               "front web": (-0.07, 0.0), "rear web": (0.07, 0.0)}
    torque_labels = []
    for g, off in offsets.items():
        w = GROUP_WALLS[g][len(GROUP_WALLS[g]) // 2]
        q = rec.at(w, _mid(an, w), F_T).resultants[2]
        torque_labels.append((w, _mid(an, w), off, f"{q / 1e3:.2f} (Bredt {q_b[w] / 1e3:.2f})"))
    plot_section_fields(an, rec, mats, [
        ("Torque shear flow N_xs under T = 2 kN·m, kN/m", F_T, 2, 1e3, "kN/m"),
        ("Axial force N_xx at the root under the 10 kN flap force", F_root, 0, 1e3, "kN/m"),
        ("Shear flow N_xs at the root under the 10 kN flap force", F_root, 2, 1e3, "kN/m"),
    ], os.path.join(HERE, "section_fields.png"), 0.035, (10, 9.0), ylim=(-0.12, 0.18), labels={0: torque_labels})

    # 3. deflection and twist along the beam
    plot_deflection(results, element_types, L, n_el, exact, os.path.join(HERE, "deflection.png"),
                    "Euler-Bernoulli misses the flap force's offset from the shear centre")

    # 4. shear lag around the root contour; interlaminar shear at mid-span
    fig, ax = plt.subplots(1, 2, figsize=(10.5, 3.9), gridspec_kw={"width_ratios": [1.6, 1]})
    order = [w.wall_id for w in an.section.walls if w.wall_id not in GROUP_WALLS["front web"] + GROUP_WALLS["rear web"]]
    starts = {w: s0 for w, s0 in zip(order, np.cumsum([0.0] + [an.section.length(an.section.wall(w)) for w in order]))}
    for i, (k, pts) in enumerate(root.items()):
        c = fs.colour(k, i)
        xs = [starts[w] + s for w, s, n in pts if w in starts]
        ns = [n / 1e3 for w, s, n in pts if w in starts]
        (fs.solves if k == "GBT" else fs.projected)(ax[0], xs, ns, color=c, label=f"{k} (wall_stresses.csv)",
                                                    ms=3.6)
    x0 = 0.5 * L / n_el     # the CSV rows are at the root element's mid-point
    F0 = [0, P * (x0 - L), 0, 0, P, T]
    unrolled_field(ax[0], an, lambda w, s: rec.at(w, s, F0).resultants[0] / 1e3, order,
                   color=fs.MUTED, label=f"Saint-Venant section field at x = {x0:.2f} m", mark_gauss=False,
                   groups={g.replace(" skin", ""): GROUP_WALLS[g] for g in GROUPS[:6]})
    ax[0].set(xlabel="distance around the contour from the leading edge (m)", ylabel="N_xx (kN/m)")
    fs.panel_title(ax[0], "Root element: N_xx around the skins", "GBT, with the root warping fixed, shows the "
                                                                "shear lag at the cap edges")
    ax[0].legend(loc="lower left")
    plot_interlaminar(ax[1], [(cap, "upper cap"), (web, "front web")])
    fig.tight_layout()
    fs.save(fig, os.path.join(HERE, "cap_stresses.png"))

    # 5. ply-axis stresses through the upper cap, ply bands shaded by material
    plot_ply_stresses(an.laminates[an.section.wall(cap_w).laminate_id], cap, os.path.join(HERE, "ply_stresses.png"),
                      "Ply-axis stresses at mid-span, upper cap centre, each ply drawn on its own",
                      {41: "GFRP", 135: "CFRP"})
    plt.close("all")


if __name__ == "__main__":
    s = main()
    print(json.dumps({k: s[k] for k in ("shear_centre", "tip", "exact_saint_venant_tip", "torque_shear_flow",
                                        "torque_shear_flow_largest_difference", "csv_vs_direct_N_xx",
                                        "root_upper_cap_peak_N_xx_gbt_over_timoshenko")}, indent=2))

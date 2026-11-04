# examples/thin_walled/box_spar/run.py
"""Example 1: composite box spar cantilever, from section analysis to ply stresses.

Run from the repository root::

    python examples/thin_walled/box_spar/run.py

Writes the job inputs to ``jobs/``, the 1d-fem results to ``results/`` (not
tracked) and the figures and ``summary.json`` next to this file.
"""

from __future__ import annotations

import json
import os
import sys

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

CFRP = Material(E1=135e9, E2=10e9, G12=5e9, nu12=0.3, G13=5e9, G23=3.8e9, rho=1550.0)
B, H = 0.6, 0.3          # box width (y) and height (z), m
L, N_EL = 6.0, 24        # span, elements
P, T = 20e3, 5e3         # tip flapwise force (z) and torque, N and N m
T_PLY = 0.5e-3


def section(n_sub: int = 4) -> SectionSpec:
    cap = [(1, a, T_PLY) for a in (45, -45, 0, 0, 0, 0, 0, 0, 0, 0, -45, 45)]   # 6 mm
    web = [(1, a, T_PLY) for a in (45, -45, -45, 45)]                          # 2 mm
    return SectionSpec(
        points={1: (-B / 2, -H / 2), 2: (B / 2, -H / 2), 3: (B / 2, H / 2), 4: (-B / 2, H / 2)},
        walls=[(1, 1, 2, 1), (2, 2, 3, 2), (3, 3, 4, 1), (4, 4, 1, 2)],  # bottom cap, rear web, top cap, front web
        materials={1: CFRP},
        laminates={1: cap, 2: web},
        n_sub=n_sub,
    )


def cases(n_el: int = N_EL, n_sub: int = 4):
    sec = section(n_sub)
    load = [0.0, 0.0, P, T, 0.0, 0.0]
    gbt_extras = None
    out = {
        "Euler-Bernoulli": CantileverSpec(L, n_el, "LinearEulerBernoulliBeamElement3D", sec, load),
        "Timoshenko": CantileverSpec(L, n_el, "LinearTimoshenkoBeamElement3D", sec, load),
        "Benscoter": CantileverSpec(L, n_el, "LinearBenscoterBeamElement3D", sec, load, warping=1),
        "GBT": CantileverSpec(L, n_el, "LinearGBTBeamElement3D", sec, load, modes=2, shear_modes=6),
    }
    from onedfem.pre_processing.section_integration.thin_walled.cross_section_modes import CrossSectionModes

    n_sv = CrossSectionModes(build_analysis(sec)).n_sv_warping
    gbt_extras = tuple(range(1, n_sv + 7))  # all warping amplitudes free at the root
    out["GBT, free warping"] = CantileverSpec(L, n_el, "LinearGBTBeamElement3D", sec, load, modes=2, shear_modes=6,
                                              free_root_extras=gbt_extras)
    return out


def exact_tip(analysis, length: float, load):
    """Tip displacement and rotation of the Saint-Venant beam by unit loads: int F_u^T C F dx."""
    C = analysis.flexibility
    xs, ws = np.polynomial.legendre.leggauss(8)
    x = 0.5 * length * (xs + 1)
    w = 0.5 * length * ws

    def forces(xx, Fz, Mx):
        return np.array([0.0, Fz * (xx - length), 0.0, 0.0, Fz, Mx])  # [N, M_y, M_z, V_y, V_z, T]

    w_tip = sum(wi * forces(xi, 1.0, 0.0) @ C @ forces(xi, load[2], load[3]) for xi, wi in zip(x, w))
    twist = sum(wi * forces(xi, 0.0, 1.0) @ C @ forces(xi, load[2], load[3]) for xi, wi in zip(x, w))
    return float(w_tip), float(twist)


def main(n_el: int = N_EL, plot: bool = True, work_dir: str = HERE) -> dict:
    from onedfem.pre_processing.section_integration.thin_walled.stress_recovery import SaintVenantRecovery

    an = build_analysis(section())
    results = {name: run(spec, work_dir, name.replace(" ", "_").replace(",", "")) for name, spec in cases(n_el).items()}
    w_ex, tw_ex = exact_tip(an, L, [0, 0, P, T])
    summary = {
        "section_stiffness_6x6": an.stiffness.tolist(),
        "shear_centre": list(an.shear_centre),
        "exact_saint_venant_tip": {"w": w_ex, "twist": tw_ex},
        "tip": {k: {"w": float(r.U[-1, 2]), "twist": float(r.U[-1, 3])} for k, r in results.items()},
    }
    # torque shear flow: 1d-fem against Bredt-Batho
    q_b = bredt_multicell(an, T)
    rec = SaintVenantRecovery(an)
    summary["torque_shear_flow"] = {
        str(w): {"bredt": q_b[w], "1d-fem": float(rec.at(w, 0.5 * an.section.length(an.section.wall(w)),
                                                        [0, 0, 0, 0, 0, T]).resultants[2])}
        for w in (1, 2)
    }
    # mid-span top cap centre: through-thickness stresses, exact section forces of the cantilever
    e_mid = n_el // 2
    x_mid = (e_mid + 0.5) * L / n_el
    F_mid = [0, P * (x_mid - L), 0, 0, P, T]
    cap = rec.at(3, B / 2, F_mid, points_per_ply=5)
    web = rec.at(2, H / 2, F_mid, points_per_ply=5)
    # dense in each ply for the through-thickness figures, drawn ply by ply
    cap_dense = rec.at(3, B / 2, F_mid, points_per_ply=21)
    web_dense = rec.at(2, H / 2, F_mid, points_per_ply=21)
    summary["mid_span_top_cap"] = {
        "N_xx": float(cap.resultants[0]), "N_xs": float(cap.resultants[2]),
        "max_sigma_11": float(max(abs(q.stress_ply[0]) for q in cap.plies)),
        "max_tau_xn": float(max(abs(q.tau_xn) for q in cap.plies)),
    }
    # the CSV of the Timoshenko job agrees with the direct recovery at the ply faces
    rows = [r for r in results["Timoshenko"].wall_rows
            if int(r["element_id"]) == e_mid and r["wall_id"] == "3" and abs(float(r["s"]) - B / 2) < 1e-9]
    summary["csv_vs_direct_N_xx"] = [float(rows[0]["N_xx"]), float(cap.resultants[0])]
    # shear lag at the root: N_xx across the top cap, GBT (plane root) against Timoshenko
    root = {k: [(float(r["s"]), float(r["N_xx"])) for r in contour(results[k].wall_rows, 0) if r["wall_id"] == "3"]
            for k in ("Timoshenko", "GBT")}
    summary["root_top_cap_N_xx"] = root
    with open(os.path.join(work_dir, "summary.json"), "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)
    if plot:
        _plots(an, rec, results, {k: c.element_type for k, c in cases(n_el).items()}, cap_dense, web_dense, root,
               n_el, (w_ex, tw_ex))
    return summary


def _plots(an, rec, results, element_types, cap, web, root, n_el, exact):
    """Figures. Fields around the section are the section's own element fields, sampled inside
    each wall element; curves along the beam are each element's shape-function field
    (``beam_curve``); through-thickness curves are drawn ply by ply; the root contour values
    from the job's CSV are shown as the discrete points they are."""
    from onedfem.post_processing.graphical_visualisers import figure_style as fs
    from onedfem.post_processing.graphical_visualisers.thin_walled_drawing import unrolled_field

    fs.use()
    import matplotlib.pyplot as plt

    sec = section()
    mats = ply_materials(sec, {1: "CFRP"})
    plot_section(an, mats, os.path.join(HERE, "section.png"), "CFRP box spar, walls at true thickness",
                 f"{an.n_sub} cubic elements per wall (ticks: element ends); caps 6 mm, webs 2 mm",
                 (B / 2 - 0.03, H / 2 - 0.025, 0.045, 0.035), "Rear web at the top cap", figsize=(10, 4.6),
                 width_ratios=(1.7, 1))
    F_T = [0, 0, 0, 0, 0, T]
    F_root = [0.0, -P * L, 0.0, 0.0, P, 0.0]
    labels = [(w, an.section.length(an.section.wall(w)) / 2, off,
               f"{rec.at(w, an.section.length(an.section.wall(w)) / 2, F_T).resultants[2] / 1e3:.2f} "
               f"(Bredt {q / 1e3:.2f})")
              for (w, off), q in zip(((1, (0, -0.06)), (2, (0.0, 0.1)), (3, (0, 0.06)), (4, (0.0, -0.1))),
                                     bredt_multicell(an, T).values())]
    plot_section_fields(an, rec, mats, [
        ("Torque shear flow N_xs, T = 5 kN·m", F_T, 2, 1e3, "kN/m"),
        ("Root N_xx, 20 kN flap force", F_root, 0, 1e3, "kN/m"),
        ("Root N_xs, 20 kN flap force", F_root, 2, 1e3, "kN/m"),
    ], os.path.join(HERE, "section_fields.png"), 0.07, (13, 4.6), ylim=(-0.27, 0.27), xlim=(-0.45, 0.45),
        labels={0: labels}, across=True)
    plot_deflection(results, element_types, L, n_el, exact, os.path.join(HERE, "deflection.png"),
                    "root warping restraint stiffens torsion by 0.2 %")

    fig, ax = plt.subplots(1, 2, figsize=(10.5, 3.9), gridspec_kw={"width_ratios": [1.3, 1]})
    x0 = 0.5 * L / n_el
    for i, (k, pts) in enumerate(root.items()):
        s_, n = np.array(pts).T
        (fs.solves if k == "GBT" else fs.projected)(ax[0], s_, n / 1e3, color=fs.colour(k, i),
                                                    label=f"{k} (wall_stresses.csv)", ms=4.5)
    unrolled_field(ax[0], an, lambda w, s_: rec.at(w, s_, [0, P * (x0 - L), 0, 0, P, T]).resultants[0] / 1e3, [3],
                   color=fs.MUTED, label=f"Saint-Venant section field at x = {x0:.3f} m", shade=False)
    ax[0].set(xlabel="s along the top cap, rear to front (m)", ylabel="N_xx (kN/m)")
    fs.panel_title(ax[0], "Root element: shear lag in the top cap",
                   "GBT with the root warping fixed loads the corners 10 % more")
    ax[0].legend(loc="lower center")
    plot_interlaminar(ax[1], [(cap, "top cap centre"), (web, "rear web centre")])
    fig.tight_layout()
    fs.save(fig, os.path.join(HERE, "cap_stresses.png"))
    plot_ply_stresses(an.laminates[1], cap, os.path.join(HERE, "ply_stresses.png"),
                      "Ply-axis stresses at mid-span, top cap centre (compression side), each ply drawn on its own",
                      {135: "CFRP"})
    plt.close("all")


if __name__ == "__main__":
    s = main()
    print(json.dumps({k: s[k] for k in ("tip", "exact_saint_venant_tip", "torque_shear_flow")}, indent=2))

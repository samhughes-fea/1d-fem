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

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))

from common import (  # noqa: E402
    CantileverSpec,
    Material,
    SectionSpec,
    bredt_multicell,
    build_analysis,
    contour,
    run,
)

GFRP = Material(E1=41e9, E2=10e9, G12=4e9, nu12=0.28, G13=4e9, G23=3.5e9, rho=1900.0)
CFRP = Material(E1=135e9, E2=10e9, G12=5e9, nu12=0.3, G13=5e9, G23=3.8e9, rho=1550.0)
L, N_EL = 10.0, 20        # span, elements
P, T = 10e3, 2e3          # tip flapwise force (z) and torque, N and N m
X_REF = 0.3               # reference axis at 30 % chord

# airfoil-like outline, chord 1 m: leading edge, upper and lower skin points, trailing edge
_PTS = {1: (0.0, 0.0), 2: (0.1, 0.06), 3: (0.25, 0.085), 4: (0.5, 0.08), 5: (1.0, 0.0),
        6: (0.5, -0.06), 7: (0.25, -0.065), 8: (0.1, -0.05)}
WALL_NAMES = {1: "LE upper skin", 2: "upper skin", 3: "upper cap", 4: "TE upper skin",
              5: "TE lower skin", 6: "lower cap", 7: "lower skin", 8: "LE lower skin",
              9: "front web", 10: "rear web"}


def section(n_sub: int = 4) -> SectionSpec:
    skin = [(1, a, 0.5e-3) for a in (0, 45, -45, -45, 45, 0)]          # 3 mm triax
    cap = [(1, 45, 0.5e-3), (1, -45, 0.5e-3)] + [(2, 0, 0.5e-3)] * 16 + [(1, -45, 0.5e-3), (1, 45, 0.5e-3)]  # 10 mm
    web = [(1, a, 0.5e-3) for a in (45, -45, -45, 45)] * 2              # 4 mm biax
    return SectionSpec(
        points={p: (y - X_REF, z) for p, (y, z) in _PTS.items()},
        walls=[(1, 1, 2, 1), (2, 2, 3, 1), (3, 3, 4, 2), (4, 4, 5, 1),
               (5, 5, 6, 1), (6, 6, 7, 2), (7, 7, 8, 1), (8, 8, 1, 1),
               (9, 3, 7, 3), (10, 4, 6, 3)],
        materials={1: GFRP, 2: CFRP},
        laminates={1: skin, 2: cap, 3: web},
        n_sub=n_sub,
    )


def cases(n_el: int = N_EL, n_sub: int = 4):
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


def main(n_el: int = N_EL, plot: bool = True, work_dir: str = HERE) -> dict:
    from pre_processing.section_integration.thin_walled.stress_recovery import SaintVenantRecovery

    an = build_analysis(section())
    results = {name: run(spec, work_dir, name.replace(" ", "_").replace(",", "")) for name, spec in cases(n_el).items()}
    v_ex, w_ex, tw_ex = exact_tip(an, L, [0, 0, P, T])
    summary = {
        "section_stiffness_6x6": an.stiffness.tolist(),
        "shear_centre": [float(c) for c in an.shear_centre],
        "exact_saint_venant_tip": {"v": v_ex, "w": w_ex, "twist": tw_ex},
        "tip": {k: {"v": float(r.U[-1, 1]), "w": float(r.U[-1, 2]), "twist": float(r.U[-1, 3])}
                for k, r in results.items()},
    }
    # torque shear flow in every wall: 1d-fem against multi-cell Bredt-Batho
    q_b = bredt_multicell(an, T)
    rec = SaintVenantRecovery(an)
    summary["torque_shear_flow"] = {
        WALL_NAMES[w]: {"bredt": q_b[w], "1d-fem": float(rec.at(w, 0.5 * an.section.length(an.section.wall(w)),
                                                                [0, 0, 0, 0, 0, T]).resultants[2])}
        for w in WALL_NAMES
    }
    # mid-span: upper cap centre and front web centre, exact section forces of the cantilever
    e_mid = n_el // 2
    x_mid = (e_mid + 0.5) * L / n_el
    F_mid = [0, P * (x_mid - L), 0, 0, P, T]
    cap_s = 0.5 * an.section.length(an.section.wall(3))
    web_s = 0.5 * an.section.length(an.section.wall(9))
    cap = rec.at(3, cap_s, F_mid, points_per_ply=5)
    web = rec.at(9, web_s, F_mid, points_per_ply=5)
    summary["mid_span_upper_cap"] = {
        "N_xx": float(cap.resultants[0]), "N_xs": float(cap.resultants[2]),
        "max_sigma_11": float(max(abs(q.stress_ply[0]) for q in cap.plies)),
        "max_tau_xn": float(max(abs(q.tau_xn) for q in cap.plies)),
    }
    summary["mid_span_front_web"] = {
        "N_xs": float(web.resultants[2]), "Q_x": float(web.transverse_shear[0]),
        "max_tau_12": float(max(abs(q.stress_ply[2]) for q in web.plies)),
    }
    rows = [r for r in results["Timoshenko"].wall_rows
            if int(r["element_id"]) == e_mid and r["wall_id"] == "3" and abs(float(r["s"]) - cap_s) < 1e-9]
    summary["csv_vs_direct_N_xx"] = [float(rows[0]["N_xx"]), float(cap.resultants[0])] if rows else None
    # root element: N_xx around the whole contour, GBT against Timoshenko
    root = {k: [(int(r["wall_id"]), float(r["s"]), float(r["N_xx"])) for r in contour(results[k].wall_rows, 0)]
            for k in ("Timoshenko", "GBT")}
    summary["root_N_xx"] = root
    with open(os.path.join(work_dir, "summary.json"), "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)
    if plot:
        _plots(an, results, cap, web, root, q_b, n_el)
    return summary


def _plots(an, results, cap, web, root, q_b, n_el):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    sec = an.section
    # section, shear centre and torque shear flow
    fig, ax = plt.subplots(figsize=(9, 2.6))
    qmax = max(abs(v) for v in q_b.values())
    for w in sec.walls:
        (y0, z0), (y1, z1) = sec.points[w.start], sec.points[w.end]
        lw = 1.5 + 6.0 * abs(q_b[w.wall_id]) / qmax
        ax.plot([y0, y1], [z0, z1], color="C3" if w.laminate_id == 2 else "C0", lw=lw, solid_capstyle="round")
        ax.annotate(f"{q_b[w.wall_id] / 1e3:.1f}", (0.5 * (y0 + y1), 0.5 * (z0 + z1)), fontsize=7,
                    ha="center", va="center", bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.8))
    ax.plot(*an.shear_centre, "k+", ms=10, label="shear centre")
    ax.plot(0, 0, "kx", ms=7, label="reference axis")
    ax.set(aspect="equal", xlabel="y (m)", ylabel="z (m)",
           title="Torque shear flow (kN/m, line width); caps in red")
    ax.legend(fontsize=7, loc="upper right")
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "section_shear_flow.png"), dpi=150)

    x = np.linspace(0, L, n_el + 1)
    fig, ax = plt.subplots(1, 2, figsize=(10, 3.8))
    for k, r in results.items():
        ax[0].plot(x, r.U[:, 2] * 1e3, label=k)
        ax[1].plot(x, np.rad2deg(r.U[:, 3]), label=k)
    ax[0].set(xlabel="x (m)", ylabel="flapwise deflection w (mm)")
    ax[1].set(xlabel="x (m)", ylabel="twist (deg)")
    ax[0].legend(fontsize=8)
    for a in ax:
        a.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "deflection.png"), dpi=150)

    fig, ax = plt.subplots(1, 2, figsize=(10, 3.8))
    for k, pts in root.items():
        upper = [(s, n) for w, s, n in pts if w == 3]
        lower = [(s, n) for w, s, n in pts if w == 6]
        for lab, data, ls in (("upper cap", upper, "-"), ("lower cap", lower, "--")):
            s, n = np.array(data).T
            ax[0].plot(s, n / 1e3, "o" + ls, label=f"{k}, {lab}", color="C0" if k == "Timoshenko" else "C1")
    ax[0].set(xlabel="s along the cap (m)", ylabel="N_xx (kN/m)", title="Root element: shear lag in the caps")
    ax[0].legend(fontsize=7)
    for q, lab in ((cap, "upper cap"), (web, "front web")):
        ax[1].plot([p.tau_xn / 1e3 for p in q.plies], [p.n * 1e3 for p in q.plies], label=f"tau_xn, {lab}")
        ax[1].plot([p.tau_sn / 1e3 for p in q.plies], [p.n * 1e3 for p in q.plies], ls="--", label=f"tau_sn, {lab}")
    ax[1].set(xlabel="interlaminar shear (kPa)", ylabel="n (mm)", title="Mid-span, through the thickness")
    ax[1].legend(fontsize=7)
    for a in ax:
        a.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "cap_stresses.png"), dpi=150)

    fig, ax = plt.subplots(1, 3, figsize=(11, 3.8), sharey=True)
    n = [p.n * 1e3 for p in cap.plies]
    for i, lab in enumerate(("sigma_11", "sigma_22", "tau_12")):
        ax[i].plot([p.stress_ply[i] / 1e6 for p in cap.plies], n)
        ax[i].set(xlabel=f"{lab} (MPa)")
        ax[i].grid(alpha=0.3)
    ax[0].set_ylabel("n through the upper cap (mm)")
    fig.suptitle("Ply-axis stresses at mid-span, upper cap centre")
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "ply_stresses.png"), dpi=150)
    plt.close("all")


if __name__ == "__main__":
    s = main()
    print(json.dumps({k: s[k] for k in ("shear_centre", "tip", "exact_saint_venant_tip", "torque_shear_flow",
                                        "csv_vs_direct_N_xx")}, indent=2))

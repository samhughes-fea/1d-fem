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
    build_analysis,
    contour,
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
    from pre_processing.section_integration.thin_walled.cross_section_modes import CrossSectionModes

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
    from pre_processing.section_integration.thin_walled.stress_recovery import SaintVenantRecovery

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
        _plots(results, cap, web, root, n_el)
    return summary


def _plots(results, cap, web, root, n_el):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    x = np.linspace(0, L, n_el + 1)
    fig, ax = plt.subplots(1, 2, figsize=(10, 3.8))
    for k, r in results.items():
        ax[0].plot(x, r.U[:, 2] * 1e3, label=k, ls="--" if "free" in k else "-")
        ax[1].plot(x, np.rad2deg(r.U[:, 3]), label=k, ls="--" if "free" in k else "-")
    ax[0].set(xlabel="x (m)", ylabel="flapwise deflection w (mm)")
    ax[1].set(xlabel="x (m)", ylabel="twist (deg)")
    ax[0].legend(fontsize=8)
    for a in ax:
        a.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "deflection.png"), dpi=150)

    fig, ax = plt.subplots(1, 2, figsize=(10, 3.8))
    for k, pts in root.items():
        s, n = np.array(pts).T
        ax[0].plot(s, n / 1e3, "o-", label=k)
    ax[0].set(xlabel="s along the top cap (m)", ylabel="N_xx (kN/m)", title="Root element: shear lag in the cap")
    ax[0].legend()
    for q, lab in ((cap, "top cap centre"), (web, "rear web centre")):
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
    ax[0].set_ylabel("n through the top cap (mm)")
    fig.suptitle("Ply-axis stresses at mid-span, top cap centre")
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "ply_stresses.png"), dpi=150)
    plt.close("all")


if __name__ == "__main__":
    s = main()
    print(json.dumps({k: s[k] for k in ("tip", "exact_saint_venant_tip", "torque_shear_flow")}, indent=2))

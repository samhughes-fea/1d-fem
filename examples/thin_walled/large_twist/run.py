# examples/thin_walled/large_twist/run.py
"""Example 4: large twist of an I-beam with restrained warping (Wagner stiffening).

Run from the repository root::

    python examples/thin_walled/large_twist/run.py

A steel I-beam cantilever, root clamped with its warping restrained, carries a tip
torque. Twisting stretches the outer fibres (``eps = r^2 phi'^2 / 2``) and the axial
stress of that stretch stiffens the beam (the Wagner effect). The torque is applied in
increments through the nonlinear static job runner with
``CorotationalBenscoterBeamElement3D``, and the tip twist is checked against the
nonlinear Vlasov equation ``GJ phi' - E Gamma phi''' + 1/2 E I_n phi'^3 = T``, which is
solved by collocation. A linear Benscoter job shows what the linear theory predicts.

Writes the job inputs to ``jobs/``, the 1d-fem results to ``results/`` (not tracked) and
``summary.json`` and ``twist.png`` next to this file.
"""

from __future__ import annotations

import csv
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))

from common import CantileverSpec, Material, SectionSpec, build_analysis, results_dir, write_job, write_supports  # noqa: E402

E, NU = 2.0e11, 0.3
G = E / (2 * (1 + NU))
STEEL = Material(E1=E, E2=E, G12=G, nu12=NU, G13=G, G23=G, rho=7850.0)
H, B, T_F, T_W = 0.3, 0.15, 6e-3, 4e-3
L, N_EL = 2.0, 16
TORQUES = (1000.0, 2000.0, 4000.0, 6000.0)


def section() -> SectionSpec:
    return SectionSpec(
        points={1: (-B / 2, H / 2), 2: (0.0, H / 2), 3: (B / 2, H / 2),
                4: (-B / 2, -H / 2), 5: (0.0, -H / 2), 6: (B / 2, -H / 2)},
        walls=[(1, 1, 2, 1), (2, 2, 3, 1), (3, 4, 5, 1), (4, 5, 6, 1), (5, 5, 2, 2)],
        materials={1: STEEL},
        laminates={1: [(1, 0.0, T_F)], 2: [(1, 0.0, T_W)]},
        n_sub=2,
    )


def vlasov_nonlinear(T: float, along: bool = False):
    """Tip twist of the nonlinear Vlasov equation (root: phi = phi' = 0, tip: phi'' = 0);
    with ``along``, the collocation solution ``phi(x)`` as a callable."""
    from numpy.polynomial.legendre import leggauss
    from scipy.integrate import solve_bvp

    an = build_analysis(section())
    GJ = 1.0 / an.flexibility[5, 5]
    A = T_F * 2 * B + T_W * H
    I_y = 2 * B * T_F * (H / 2) ** 2 + T_W * H**3 / 12
    I_z = 2 * T_F * B**3 / 12
    EG = E * T_F * B**3 * H**2 / 24
    x, w = leggauss(10)
    Irr = sum(T_F * 0.5 * B * np.sum(w * ((0.5 * B * x) ** 2 + z0**2) ** 2) for z0 in (H / 2, -H / 2))
    Irr += T_W * 0.5 * H * np.sum(w * (0.5 * H * x) ** 4)
    I_n = Irr - (I_y + I_z) ** 2 / A

    def rhs(s, y):
        return np.vstack([y[1], y[2], (GJ * y[1] + 0.5 * E * I_n * y[1] ** 3 - T) / EG])

    s = np.linspace(0, L, 200)
    y0 = np.zeros((3, s.size))
    y0[1] = T / GJ * (1 - np.exp(-5 * s))
    sol = solve_bvp(rhs, lambda a, b: np.array([a[0], a[1], b[2]]), s, y0, tol=1e-10, max_nodes=100000)
    if along:
        return lambda xx: sol.sol(xx)[0]
    return float(sol.sol(L)[0])


def _run_job(work_dir, name, etype, T, nonlinear):
    import logging

    from onedfem.workflow_orchestrator.run_job import process_job

    spec = CantileverSpec(L, N_EL, etype, section(), [0.0, 0.0, 0.0, T, 0.0, 0.0], warping=1)
    job = os.path.join(work_dir, "jobs", name)
    write_job(job, spec)
    write_supports(job, spec)
    if nonlinear:
        with open(os.path.join(job, "simulation_settings.txt"), "w", encoding="utf-8") as f:
            f.write("[Simulation]\n[Type]\nstatic_nonlinear\n\n[Newton]\ntolerance = 1e-4\n"
                    "relative_tolerance = 1e-8\nrelative_reference = first_residual\nmax_iterations = 40\n"
                    "tolerance_delta_u = 1e-10\n\n[Solver]\ntype = direct\ntolerance = 1e-12\n\n"
                    "[Nonlinear]\nnum_increments = 8\nline_search = false\n")
    res = os.path.join(work_dir, "results", name)
    logging.getLogger().setLevel(logging.WARNING)
    process_job(job, results_dir(res), {}, {}, force_serial=True, max_processes_per_job=1)
    if nonlinear:  # the runner stops at a failed increment with only a warning: check it reached the full load
        with open(os.path.join(res, "logs", "newton_history.csv"), newline="") as fh:
            last = list(csv.DictReader(fh))[-1]
        if float(last["load_factor"]) < 1.0 - 1e-12 or last["residual_ok"] != "1":
            raise RuntimeError(f"{name}: Newton did not reach the full load (last row {last})")
    with open(os.path.join(res, "primary_results", "global", "U_global.csv"), newline="") as fh:
        U = np.array([float(r[1]) for r in list(csv.reader(fh))[1:]]).reshape(N_EL + 1, -1)
    _FIELDS[name] = U
    return float(U[-1, 3])


_FIELDS: dict = {}


def main(n_el: int = N_EL, plot: bool = True, work_dir: str = HERE, torques=TORQUES) -> dict:
    global N_EL
    N_EL = n_el
    rows = []
    for T in torques:
        rows.append({
            "torque": T,
            "corotational": _run_job(work_dir, f"cr_{T:.0f}", "CorotationalBenscoterBeamElement3D", T, True),
            "linear": _run_job(work_dir, f"lin_{T:.0f}", "LinearBenscoterBeamElement3D", T, False),
            "vlasov_nonlinear": vlasov_nonlinear(T),
        })
    summary = {"rows": rows, "max_error": max(abs(r["corotational"] / r["vlasov_nonlinear"] - 1) for r in rows)}
    with open(os.path.join(work_dir, "summary.json"), "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)
    if plot:
        _plots(rows, torques, work_dir)
    return summary


def _plots(rows, torques, out_dir):
    """Torque against tip twist (one marker per solve), twist along the beam through each
    element's interpolation, and the tip section turned by the computed twist."""
    from onedfem.post_processing.graphical_visualisers import figure_style as fs
    from onedfem.post_processing.graphical_visualisers.beam_field_sampling import beam_curve
    from onedfem.post_processing.graphical_visualisers.thin_walled_drawing import draw_section

    from common import ply_materials

    fs.use()
    import matplotlib.pyplot as plt

    fig = plt.figure(figsize=(12, 3.9))
    gs = fig.add_gridspec(1, 3, width_ratios=[1, 1.25, 0.8])
    ax = fig.add_subplot(gs[0])
    Ts = np.linspace(0, max(torques), 40)
    fs.analytic(ax, [0] + [vlasov_nonlinear(t) for t in Ts[1:]], Ts / 1e3, label="nonlinear Vlasov (collocation)")
    fs.solves(ax, [r["corotational"] for r in rows], [r["torque"] / 1e3 for r in rows],
              color=fs.colour("corotational"), label="corotational Benscoter, one job per torque")
    fs.solves(ax, [r["linear"] for r in rows], [r["torque"] / 1e3 for r in rows], color=fs.colour("linear"),
              marker="s", label="linear Benscoter")
    ax.set(xlabel="tip twist (rad)", ylabel="tip torque (kN·m)")
    fs.panel_title(ax, "Wagner stiffening", "27 % less twist than linear at 6 kN·m")
    ax.legend(loc="lower right")

    ax = fig.add_subplot(gs[1])
    x = np.linspace(0, L, N_EL + 1)
    xa = np.linspace(0, L, 101)
    for i, r in enumerate(rows):
        c = fs.colour("", i)
        U = _FIELDS.get(f"cr_{r['torque']:.0f}")
        if U is None:
            continue
        xs, v = beam_curve(x, U, 3, "LinearBenscoterBeamElement3D")
        fs.field(ax, xs, v, color=c, label=f"{r['torque'] / 1e3:.0f} kN·m")
        fs.nodal(ax, x, U[:, 3], color=c)
        fs.analytic(ax, xa, vlasov_nonlinear(r["torque"], along=True)(xa), label="Vlasov" if i == 0 else None,
                    lw=0.8)
    ax.set(xlabel="x (m)", ylabel="twist φ (rad)")
    fs.panel_title(ax, "Twist along the beam, corotational", f"{N_EL} elements, linear interpolation in each; "
                                                            "warping restrained at the root")
    ax.legend(loc="upper left")
    fs.marks_legend(ax, ("field", "nodal", "analytic"), loc="lower right")

    ax = fig.add_subplot(gs[2])
    an = build_analysis(section())
    mats = ply_materials(section(), {1: "steel"})
    draw_section(ax, an, mats, legend=False, nodes=False)
    phi = rows[-1]["corotational"]
    R = np.array([[np.cos(phi), -np.sin(phi)], [np.sin(phi), np.cos(phi)]])
    for w in an.section.walls:
        a_, b_ = an.section.endpoints(w)
        p = np.array([R @ a_, R @ b_])
        ax.plot(p[:, 0], p[:, 1], color=fs.colour("corotational"), lw=1.8, zorder=5)
    ax.set(xlabel="y (m)", ylabel="z (m)", xlim=(-0.2, 0.2), ylim=(-0.2, 0.2))
    fs.panel_title(ax, "Tip section", f"turned by the computed {phi:.3f} rad at {rows[-1]['torque'] / 1e3:.0f} kN·m")
    fig.tight_layout()
    fs.save(fig, os.path.join(out_dir, "twist.png"))


if __name__ == "__main__":
    s = main()
    for r in s["rows"]:
        print(f"T = {r['torque']:6.0f} N m: corotational {r['corotational']:.4f} rad, nonlinear Vlasov "
              f"{r['vlasov_nonlinear']:.4f} rad, linear {r['linear']:.4f} rad")
    print(f"largest difference from the nonlinear Vlasov solution: {100 * s['max_error']:.2f} %")

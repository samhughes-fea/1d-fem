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


def vlasov_nonlinear(T: float) -> float:
    """Tip twist of the nonlinear Vlasov equation (root: phi = phi' = 0, tip: phi'' = 0)."""
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
    return float(sol.sol(L)[0])


def _run_job(work_dir, name, etype, T, nonlinear):
    import logging

    from workflow_orchestrator.run_job import process_job

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
    return float(U[-1, 3])


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
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        Ts = np.linspace(0, max(torques), 40)
        fig, ax = plt.subplots(figsize=(5.5, 3.8))
        ax.plot([0] + [vlasov_nonlinear(t) for t in Ts[1:]], Ts / 1e3, "k-", label="nonlinear Vlasov")
        ax.plot([r["corotational"] for r in rows], [r["torque"] / 1e3 for r in rows], "o",
                label="1d-fem corotational Benscoter")
        ax.plot([r["linear"] for r in rows], [r["torque"] / 1e3 for r in rows], "s", mfc="none",
                label="1d-fem linear Benscoter")
        ax.set_xlabel("tip twist (rad)")
        ax.set_ylabel("tip torque (kN m)")
        ax.set_title("I-beam, warping restrained at the root")
        ax.legend()
        fig.tight_layout()
        fig.savefig(os.path.join(work_dir, "twist.png"), dpi=130)
        plt.close(fig)
    return summary


if __name__ == "__main__":
    s = main()
    for r in s["rows"]:
        print(f"T = {r['torque']:6.0f} N m: corotational {r['corotational']:.4f} rad, nonlinear Vlasov "
              f"{r['vlasov_nonlinear']:.4f} rad, linear {r['linear']:.4f} rad")
    print(f"largest difference from the nonlinear Vlasov solution: {100 * s['max_error']:.2f} %")

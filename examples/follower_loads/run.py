"""Follower loads: a cantilever under a tip force that stays perpendicular to the beam.

Run from the repository root::

    python examples/follower_loads/run.py

A ``follower_load.txt`` (rows ``x y z F_x F_y F_z M_x M_y M_z`` under ``[Follower load]``)
puts loads on nodes that turn with the node. A nonlinear static job of corotational
Timoshenko elements carries a tip force that stays normal to the deformed tip, up to
``P L^2 / EI = 6`` in load steps. The reference is the planar elastica with a follower
force, solved by shooting with continuation in the load, and the dead-load
solution of the same beam is shown alongside.

Writes the job inputs to ``jobs/`` and results to ``results/`` (not tracked), and the figure
and ``summary.json`` next to this file.
"""

from __future__ import annotations

import csv
import json
import logging
import os
import shutil
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
if REPO not in sys.path:
    sys.path.insert(0, REPO)

L, E, A, I = 1.0, 1e9, 1e-2, 1e-7      # slender: shear and stretch are negligible
EI = E * I
N_EL = 40
ELEMENT = "CorotationalTimoshenkoBeamElement3D"
STEPS = (1.0, 2.0, 3.0, 4.5, 6.0)     # P L^2 / EI, one job each


def elastica(k: float, follower: bool, n: int = 401):
    """Planar cantilever under a tip force ``k EI / L^2``, initially along ``+y``.

    ``theta(s)`` is the slope; ``EI theta' = M``, ``M' = -(t x f)_z`` with ``t = (cos, sin)``.
    The force is ``(0, P)`` (dead) or ``P (-sin theta_L, cos theta_L)`` (follower). Solved by
    shooting from the tip (``theta(L) = theta_L``, ``M(L) = 0``) for ``theta(0) = 0``, with
    continuation in the load. Returns ``s, x, y, theta``.
    """
    from scipy.integrate import solve_ivp
    from scipy.optimize import brentq

    s = np.linspace(0.0, L, n)

    def shoot(tl, P, dense=False):
        fx, fy = (-P * np.sin(tl), P * np.cos(tl)) if follower else (0.0, P)

        def rhs(_, q):
            th, m = q
            return [m / EI, -(np.cos(th) * fy - np.sin(th) * fx)]

        sol = solve_ivp(rhs, (L, 0.0), [tl, 0.0], rtol=1e-12, atol=1e-14, dense_output=dense)
        return sol if dense else sol.y[0, -1]

    tl = 0.0
    for kk in np.linspace(0.0, k, int(np.ceil(k / 0.1)) + 1)[1:]:
        P = kk * EI / L ** 2
        lo, hi = tl - 1e-9, tl + 0.5
        while shoot(lo, P) * shoot(hi, P) > 0:
            hi += 0.5
        tl = brentq(lambda t: shoot(t, P), lo, hi, xtol=1e-14)
    sol = shoot(tl, k * EI / L ** 2, dense=True)
    th = sol.sol(s)[0]
    ds = np.diff(s)
    x = np.concatenate([[0.0], np.cumsum(0.5 * ds * (np.cos(th[1:]) + np.cos(th[:-1])))])
    y = np.concatenate([[0.0], np.cumsum(0.5 * ds * (np.sin(th[1:]) + np.sin(th[:-1])))])
    return s, x, y, th


def write_job(job_dir: str, k_max: float, follower: bool) -> str:
    if os.path.isdir(job_dir):
        shutil.rmtree(job_dir)
    os.makedirs(job_dir)

    def w(name, text):
        with open(os.path.join(job_dir, name), "w", encoding="utf-8") as f:
            f.write(text)

    n = N_EL
    P = k_max * EI / L ** 2
    w("grid.txt", "[Grid]\n[node_id] [x] [y] [z]\n" + "".join(f"{i} {i * L / n:.12g} 0 0\n" for i in range(n + 1)))
    w("element.txt", "[Element]\n[element_id] [node1] [node2] [element_type] [axial_order] [bending_y_order] "
      "[bending_z_order] [shear_y_order] [shear_z_order] [torsion_order] [load_order]\n"
      + "".join(f"{e} {e} {e + 1} {ELEMENT} 2 2 2 1 1 2 2\n" for e in range(n)))
    w("material.txt", "[Material]\n[element_id] [E] [G] [nu] [rho]\n"
      + "".join(f"{e} {E:g} {E / 2.6:.10g} 0.3 1000\n" for e in range(n)))
    w("section.txt", "[Section]\n[element_id] [A] [I_x] [I_y] [I_z] [J_t]\n"
      + "".join(f"{e} {A:g} 0 {I:g} {I:g} {2 * I:g}\n" for e in range(n)))
    row = f"{L:g} 0 0 0 {P:.12g} 0 0 0 0\n"
    if follower:
        w("point_load.txt", "[Point load]\n[x] [y] [z] [F_x] [F_y] [F_z] [M_x] [M_y] [M_z]\n")
        w("follower_load.txt", "[Follower load]\n" + row)
    else:
        w("point_load.txt", "[Point load]\n[x] [y] [z] [F_x] [F_y] [F_z] [M_x] [M_y] [M_z]\n" + row)
    w("prescribed_displacement.txt", "[Prescribed Displacement]\n[id] [node_id] [dof] [value] [type]\n"
      + "".join(f"{i} 0 {d} 0.0 displacement\n" for i, d in enumerate(("UX", "UY", "UZ", "RX", "RY", "RZ"))))
    # load steps of P L^2 / EI = 0.5
    factors = ", ".join(f"{f:.12g}" for f in np.linspace(0.0, 1.0, int(np.ceil(k_max / 0.5)) + 1)[1:])
    w("simulation_settings.txt", "[Simulation]\n[Type]\nstatic_nonlinear\n\n[Newton]\ntolerance = 1e-6\n"
      "relative_tolerance = 1e-10\nrelative_reference = first_residual\nmax_iterations = 40\n"
      "tolerance_delta_u = 1e-12\n\n[Solver]\ntype = direct\ntolerance = 1e-12\n\n"
      f"[Nonlinear]\nload_factors = {factors}\nline_search = false\n")
    return job_dir


def run_job(job_dir: str, work: str) -> np.ndarray:
    """Run the job; return ``U`` (n_nodes, 6) at the last load step."""
    from workflow_orchestrator.run_job import process_job

    res = os.path.join(work, "results", os.path.basename(job_dir))
    if os.path.isdir(res):
        shutil.rmtree(res)
    for sub in ("element_stiffness_matrices", "element_force_vectors", "primary_results", "secondary_results",
                "logs", "maps", "diagnostics"):
        os.makedirs(os.path.join(res, sub))
    process_job(job_dir, res, {}, {}, force_serial=True, max_processes_per_job=1)
    with open(os.path.join(res, "primary_results", "global", "U_global.csv"), newline="", encoding="utf-8") as f:
        U = np.array([float(r[1]) for r in list(csv.reader(f))[1:]])
    return U.reshape(N_EL + 1, -1)


def main(plot: bool = True, work_dir: str = HERE, out_dir: str = HERE) -> dict:
    logging.disable(logging.WARNING)
    try:
        out = {"steps": list(STEPS), "follower": [], "dead": []}
        for follower in (True, False):
            key = "follower" if follower else "dead"
            for k in STEPS:
                U = run_job(write_job(os.path.join(work_dir, "jobs", f"{key}_{k:g}"), k, follower), work_dir)
                _, xr, yr, thr = elastica(k, follower)
                tip = U[-1]
                out[key].append({
                    "k": k, "u_x": float(tip[0]), "u_y": float(tip[1]), "theta": float(tip[5]),
                    "ref_u_x": float(xr[-1] - L), "ref_u_y": float(yr[-1]), "ref_theta": float(thr[-1]),
                    "nodes": (np.array([np.linspace(0, L, N_EL + 1) + U[:, 0], U[:, 1]]).T).tolist()})
    finally:
        logging.disable(logging.NOTSET)
    for key in ("follower", "dead"):
        for r in out[key]:
            r["error_tip_over_L"] = float(np.hypot(r["u_x"] - r["ref_u_x"], r["u_y"] - r["ref_u_y"]) / L)
            r["error_theta"] = float(r["theta"] - r["ref_theta"])
    with open(os.path.join(out_dir, "summary.json"), "w", encoding="utf-8") as f:
        json.dump({k: ([{kk: vv for kk, vv in r.items() if kk != "nodes"} for r in v] if isinstance(v, list) and v
                       and isinstance(v[0], dict) else v) for k, v in out.items()}, f, indent=1)
    if plot:
        figures(out, out_dir)
    return out


def figures(out: dict, out_dir: str) -> None:
    """Nodes of the deformed beam are triangles (explicit nodal values); the elastica
    reference is a line from its own ODE solution."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(10, 4.2))
    for ax, key, title in zip(axes, ("follower", "dead"), ("follower tip force (stays normal to the tip)",
                                                          "dead tip force (stays vertical)")):
        for i, r in enumerate(out[key]):
            _, xr, yr, _ = elastica(r["k"], key == "follower")
            ax.plot(xr, yr, "-", color="#c3c2b7", lw=1.0, label="elastica (ODE solution)" if i == 0 else None)
            nd = np.array(r["nodes"])
            ax.plot(nd[:, 0], nd[:, 1], "^", ms=3, color="#2a78d6",
                    label=f"1d-fem nodes, {N_EL} corotational Timoshenko elements" if i == 0 else None)
        ax.set_aspect("equal")
        ax.set_title(f"P L² / EI = {', '.join(f'{k:g}' for k in out['steps'])}: {title}", fontsize=8, loc="left")
        ax.set_xlabel("x / L", fontsize=8)
        ax.set_ylabel("y / L", fontsize=8)
        ax.tick_params(labelsize=7)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
    axes[1].legend(fontsize=7, frameon=False, loc="upper left")
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, "follower_cantilever.png"), dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    o = main()
    for key in ("follower", "dead"):
        print(key, [(r["k"], round(r["error_tip_over_L"], 6), round(r["error_theta"], 6)) for r in o[key]])

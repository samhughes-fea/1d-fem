# examples/thin_walled/channel_buckling/run.py
"""Example 3: flexural-torsional buckling and P-delta of a channel column.

Run from the repository root::

    python examples/thin_walled/channel_buckling/run.py

A steel channel with fork supports at both ends. Its shear centre and centroid are
both off the web, where the beam axis (the section origin) is placed, so every
second-order coupling of the section is exercised:

1. Linear buckling with Benscoter elements under compression at the centroid: the
   flexural-torsional load against Timoshenko and Gere's closed form.
2. P-delta with Timoshenko elements under compression at the web, which is eccentric to
   the tension centre: the midspan deflection against the secant formula, for a range of
   loads, with the linear solution for comparison.

Writes the job inputs to ``jobs/``, the 1d-fem results to ``results/`` (not tracked) and
``summary.json`` and ``p_delta.png`` next to this file.
"""

from __future__ import annotations

import csv
import glob
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))

from common import CantileverSpec, Material, SectionSpec, build_analysis, results_dir, write_job  # noqa: E402

E, NU = 2.0e11, 0.3
G = E / (2 * (1 + NU))
STEEL = Material(E1=E, E2=E, G12=G, nu12=NU, G13=G, G23=G, rho=7850.0)
H, B, T_W = 0.2, 0.08, 4e-3     # web depth (z), flange width (y), wall thickness, m
L_BUCKLE, L_PDELTA, N_EL = 4.0, 6.0, 32


def section() -> SectionSpec:
    """Web on y = 0 (the beam axis), flanges towards +y."""
    return SectionSpec(
        points={1: (B, -H / 2), 2: (0.0, -H / 2), 3: (0.0, H / 2), 4: (B, H / 2)},
        walls=[(1, 1, 2, 1), (2, 2, 3, 1), (3, 3, 4, 1)],
        materials={1: STEEL},
        laminates={1: [(1, 0.0, T_W)]},
        n_sub=2,
    )


def closed_form_ft(length: float) -> dict:
    """Flexural-torsional buckling of the channel from mid-line formulas (TG sec. 5.4)."""
    t, h, b = T_W, H, B
    A = t * (h + 2 * b)
    y_c = b * b / (h + 2 * b)
    y_s = -3 * b * b / (h + 6 * b)
    I_y = t * h**3 / 12 + 2 * b * t * (h / 2) ** 2
    I_z = 2 * t * b**3 / 3 - A * y_c**2
    J = t**3 * (h + 2 * b) / 3
    Gamma = t * b**3 * h**2 / 12 * (3 * b + 2 * h) / (6 * b + h)
    y0 = y_s - y_c
    r02 = (I_y + I_z) / A + y0**2
    P_w = np.pi**2 * E * I_y / length**2
    P_v = np.pi**2 * E * I_z / length**2
    P_phi = (G * J + np.pi**2 * E * Gamma / length**2) / r02
    P_ft = float(min(np.roots([1 - y0**2 / r02, -(P_w + P_phi), P_w * P_phi])))
    return dict(A=A, y_c=y_c, y_s=y_s, I_y=I_y, I_z=I_z, J=J, Gamma=Gamma, P_v=P_v, P_w=P_w, P_phi=P_phi,
                P_ft=P_ft)


def _write_fork_job(job_dir, etype, length, n_el, warping, loads, settings):
    write_job(job_dir, CantileverSpec(length, n_el, etype, section(), [0.0] * 6, warping=warping))
    with open(os.path.join(job_dir, "point_load.txt"), "w", encoding="utf-8") as f:
        f.write("[Point load]\n[x] [y] [z] [F_x] [F_y] [F_z] [M_x] [M_y] [M_z]\n")
        for x, vals in loads:
            f.write(f"{x:.9g} 0 0 " + " ".join(f"{v:.9g}" for v in vals) + "\n")
    rows = [(0, d) for d in ("UX", "UY", "UZ", "RX")] + [(n_el, d) for d in ("UY", "UZ", "RX")]
    with open(os.path.join(job_dir, "prescribed_displacement.txt"), "w", encoding="utf-8") as f:
        f.write("[Prescribed Displacement]\n[id] [node_id] [dof] [value] [type] [comment]\n")
        f.writelines(f"{i} {n} {d} 0.0 displacement\n" for i, (n, d) in enumerate(rows))
    with open(os.path.join(job_dir, "simulation_settings.txt"), "w", encoding="utf-8") as f:
        f.write(settings)


def _run(job_dir, res_dir):
    import logging

    from workflow_orchestrator.run_job import process_job

    logging.getLogger().setLevel(logging.WARNING)
    process_job(job_dir, results_dir(res_dir), {}, {}, force_serial=True, max_processes_per_job=1)


def _compression_at(y, P):
    """End loads of a compression ``P`` acting at ``(y, 0)`` of both end sections."""
    return [(0.0, [0.0, 0.0, 0.0, 0.0, 0.0, -y * P]), (None, [-P, 0.0, 0.0, 0.0, 0.0, y * P])]


def buckling(n_el: int, work_dir: str) -> dict:
    ref = closed_form_ft(L_BUCKLE)
    an = build_analysis(section())
    y_t = an.tension_centre[0]
    loads = [(L_BUCKLE if x is None else x, v) for x, v in _compression_at(y_t, 1.0)]
    settings = ("[Simulation]\n[Type]\nbuckling\n\n[Solver]\ntype = direct\ntolerance = 1e-12\n\n"
                "[Buckling]\nenabled = true\nnum_modes = 3\nbuckling_prestress = linear_static\n"
                "buckling_load_factor = 1.0\n")
    job = os.path.join(work_dir, "jobs", "buckling")
    _write_fork_job(job, "LinearBenscoterBeamElement3D", L_BUCKLE, n_el, 1, loads, settings)
    res = os.path.join(work_dir, "results", "buckling")
    _run(job, res)
    f = glob.glob(os.path.join(res, "primary_results", "modal_results", "*_buckling_load_factors.txt"))[0]
    lam = np.atleast_1d(np.loadtxt(f))
    modes = np.loadtxt(f.replace("load_factors", "mode_shapes")).reshape(-1, lam.size)
    stride = modes.shape[0] // (n_el + 1)
    # flexure along the symmetry axis is uncoupled; the flexural-torsional mode twists
    twists = [k for k in np.argsort(lam) if np.abs(modes[3::stride, k]).max() > 1e-6 * np.abs(modes[:, k]).max()]
    lam_ft = float(lam[twists[0]])
    lam = np.sort(lam)
    return {
        "shear_centre_y": an.shear_centre[0], "tension_centre_y": y_t,
        "closed_form": {k: ref[k] for k in ("P_v", "P_w", "P_phi", "P_ft")},
        "load_factors": lam.tolist(),
        "P_ft": lam_ft, "P_ft_error": lam_ft / ref["P_ft"] - 1.0,
    }


def p_delta(n_el: int, work_dir: str, ratios=(0.2, 0.4, 0.6, 0.8)) -> dict:
    an = build_analysis(section())
    EIz = 1.0 / an.flexibility[2, 2]
    P_v = np.pi**2 * EIz / L_PDELTA**2
    e = -an.tension_centre[0]                     # load at the web: eccentric to the tension centre
    rows = []
    for r in ratios:
        P = r * P_v
        loads = [(L_PDELTA if x is None else x, v) for x, v in _compression_at(0.0, P)]
        out = {}
        for so in ("none", "p_delta"):
            settings = ("[Simulation]\n[Type]\nstatic\n\n[Solver]\ntype = direct\ntolerance = 1e-12\n"
                        f"max_iterations = 500\n\n[Static]\nsecond_order = {so}\np_delta_tolerance = 1e-10\n")
            name = f"p_delta_{r:.1f}_{so}"
            job = os.path.join(work_dir, "jobs", name)
            _write_fork_job(job, "LinearTimoshenkoBeamElement3D", L_PDELTA, n_el, 0, loads, settings)
            res = os.path.join(work_dir, "results", name)
            _run(job, res)
            with open(os.path.join(res, "primary_results", "global", "U_global.csv"), newline="") as fh:
                U = np.array([float(row[1]) for row in list(csv.reader(fh))[1:]]).reshape(n_el + 1, -1)
            out[so] = float(U[n_el // 2, 1] - U[0, 1])
        k = np.sqrt(P / EIz)
        rows.append({"P_over_Pv": r, "linear": out["none"], "p_delta": out["p_delta"],
                     "secant": float(-e * (1 / np.cos(k * L_PDELTA / 2) - 1))})
    return {"eccentricity": e, "P_v": P_v, "rows": rows,
            "max_error": max(abs(rw["p_delta"] / rw["secant"] - 1) for rw in rows)}


def main(n_el: int = N_EL, plot: bool = True, work_dir: str = HERE) -> dict:
    summary = {"buckling": buckling(n_el, work_dir), "p_delta": p_delta(n_el, work_dir)}
    with open(os.path.join(work_dir, "summary.json"), "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)
    if plot:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        rows = summary["p_delta"]["rows"]
        r = [rw["P_over_Pv"] for rw in rows]
        fig, ax = plt.subplots(figsize=(5.5, 3.8))
        ax.plot(r, [1e3 * rw["secant"] for rw in rows], "k-", label="secant formula")
        ax.plot(r, [1e3 * rw["p_delta"] for rw in rows], "o", label="1d-fem P-delta")
        ax.plot(r, [1e3 * rw["linear"] for rw in rows], "s", mfc="none", label="1d-fem linear")
        ax.set_xlabel("P / P_v")
        ax.set_ylabel("midspan deflection v (mm)")
        ax.set_title("Channel compressed at the web (eccentric to the tension centre)")
        ax.legend()
        fig.tight_layout()
        fig.savefig(os.path.join(work_dir, "p_delta.png"), dpi=130)
        plt.close(fig)
    return summary


if __name__ == "__main__":
    s = main()
    b, p = s["buckling"], s["p_delta"]
    print(f"flexural-torsional buckling: {b['P_ft']:.1f} N vs closed form {b['closed_form']['P_ft']:.1f} N "
          f"({100 * b['P_ft_error']:+.2f} %)")
    print(f"P-delta: largest error against the secant formula {100 * p['max_error']:.2f} %")

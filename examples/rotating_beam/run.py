"""Rotating beams: centrifugal stiffening, spin softening and Coriolis terms, against exact results.

Run from the repository root::

    python examples/rotating_beam/run.py

Each case is an ordinary 1d-fem eigen job with a ``[Rotation]`` block in
``simulation_settings.txt``; ``omega_sweep`` makes the job write a Campbell table
(``primary_results/rotation/campbell.csv``, rotating-frame frequencies at each speed).

1. **Flap and lag of a uniform cantilever rotating about its root** against Hodges and
   Rutkowski (1981). Flap: their tabulated first-mode frequencies. Lag: the exact relation
   ``omega_lag^2 = omega_flap^2 - Omega^2`` for the same bending stiffness.
2. **Whirl of a shaft spinning about its own axis.** In the rotating frame the Coriolis
   terms split each bending pair into ``omega_n - Omega`` and ``omega_n + Omega``.

Writes the job inputs to ``jobs/`` and results to ``results/`` (not tracked), and the
figures and ``summary.json`` next to this file.
"""

from __future__ import annotations

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

L, E, RHO, A = 1.0, 1e9, 1000.0, 1e-3
N_EL = 40
ELEMENT = "LinearEulerBernoulliBeamElement3D"
# Hodges & Rutkowski (1981): first flap frequency of a uniform cantilever rotating about its
# root, nondimensionalised by sqrt(EI / (m L^4)), at nondimensional speeds Omega_bar
HODGES_FLAP = {0: 3.5160, 1: 3.6816, 2: 4.1373, 3: 4.7973, 4: 5.5850, 5: 6.4495, 6: 7.3604, 7: 8.2996,
               8: 9.2568, 9: 10.2257, 10: 11.2023, 11: 12.1843, 12: 13.1702}


def write_job(job_dir: str, Iy: float, Iz: float, rotation: str, n: int = N_EL, num_modes: int = 4) -> str:
    """A cantilever along ``x`` clamped at ``x = 0``, with a ``[Rotation]`` block."""
    if os.path.isdir(job_dir):
        shutil.rmtree(job_dir)
    os.makedirs(job_dir)

    def w(name, text):
        with open(os.path.join(job_dir, name), "w", encoding="utf-8") as f:
            f.write(text)

    w("grid.txt", "[Grid]\n[node_id] [x] [y] [z]\n" + "".join(f"{i} {i * L / n:.12g} 0 0\n" for i in range(n + 1)))
    w("element.txt", "[Element]\n[element_id] [node1] [node2] [element_type] [axial_order] [bending_y_order] "
      "[bending_z_order] [shear_y_order] [shear_z_order] [torsion_order] [load_order]\n"
      + "".join(f"{e} {e} {e + 1} {ELEMENT} 3 4 4 0 0 3 2\n" for e in range(n)))
    w("material.txt", "[Material]\n[element_id] [E] [G] [nu] [rho]\n"
      + "".join(f"{e} {E:g} {E / 2.6:.10g} 0.3 {RHO:g}\n" for e in range(n)))
    w("section.txt", "[Section]\n[element_id] [A] [I_x] [I_y] [I_z] [J_t]\n"
      + "".join(f"{e} {A:g} 0 {Iy:g} {Iz:g} {Iy + Iz:g}\n" for e in range(n)))
    w("point_load.txt", "[Point load]\n[x] [y] [z] [F_x] [F_y] [F_z] [M_x] [M_y] [M_z]\n")
    w("prescribed_displacement.txt", "[Prescribed Displacement]\n[id] [node_id] [dof] [value] [type]\n"
      + "".join(f"{i} 0 {d} 0.0 displacement\n" for i, d in enumerate(("UX", "UY", "UZ", "RX", "RY", "RZ"))))
    w("simulation_settings.txt", f"[Simulation]\n[Type]\neigen\n\n[Eigen]\nnum_modes = {num_modes}\n\n"
      f"[Rotation]\n{rotation}")
    return job_dir


def run(job_dir: str, work: str) -> np.ndarray:
    """Run an eigen job; return its Campbell table (rows ``[omega, f_1, ...]``)."""
    from onedfem.workflow_orchestrator.run_job import process_job

    res = os.path.join(work, "results", os.path.basename(job_dir))
    if os.path.isdir(res):
        shutil.rmtree(res)
    for sub in ("element_stiffness_matrices", "element_force_vectors", "primary_results", "secondary_results",
                "logs", "maps", "diagnostics"):
        os.makedirs(os.path.join(res, sub))
    process_job(job_dir, res, {}, {}, force_serial=True, max_processes_per_job=1)
    return np.loadtxt(os.path.join(res, "primary_results", "rotation", "campbell.csv"), delimiter=",", skiprows=1)


def hodges(work: str) -> dict:
    """Flap (rotor axis ``z``, flap is ``w``) and lag (``v``) against Hodges and Rutkowski.

    The reference is an inextensible beam. Here the centrifugal stretch adds
    ``Omega_bar^2 I / (A L^2)`` to the tension, so the section is slender
    (``I / (A L^2) = 1e-6``) to make that negligible.
    """
    Iy = 1e-9
    s = np.sqrt(E * Iy / (RHO * A * L ** 4))            # rad/s per unit nondimensional frequency
    speeds = sorted(HODGES_FLAP)
    sweep = ", ".join(f"{b * s:.12g}" for b in speeds)
    rot = f"axis = 0, 0, 1\norigin = 0, 0, 0\ngyroscopic = false\nomega_sweep = {sweep}\n"
    # flap: lag 10 times stiffer so the first mode is flap; lag: the same with y and z swapped
    flap = run(write_job(os.path.join(work, "jobs", "flap"), Iy, 100 * Iy, rot), work)
    lag = run(write_job(os.path.join(work, "jobs", "lag"), 100 * Iy, Iy, rot), work)
    f_flap = flap[:, 1] * 2 * np.pi / s
    f_lag = lag[:, 1] * 2 * np.pi / s
    ref_flap = np.array([HODGES_FLAP[b] for b in speeds])
    ref_lag = np.sqrt(ref_flap ** 2 - np.array(speeds, dtype=float) ** 2)
    return {"omega_bar": speeds, "flap": f_flap.tolist(), "lag": f_lag.tolist(), "hodges_flap": ref_flap.tolist(),
            "lag_exact": ref_lag.tolist(), "flap_difference": (f_flap / ref_flap - 1).tolist(),
            "lag_difference": (f_lag / ref_lag - 1).tolist()}


def whirl(work: str) -> dict:
    """A round shaft spinning about its own axis: forward and backward whirl."""
    I = 1e-8                                           # slender, little rotary inertia
    w1 = 1.875104 ** 2 * np.sqrt(E * I / (RHO * A * L ** 4))
    fracs = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8]
    sweep = ", ".join(f"{f * w1:.12g}" for f in fracs)
    rows = run(write_job(os.path.join(work, "jobs", "whirl"), I, I,
                         f"axis = 1, 0, 0\norigin = 0, 0, 0\ngyroscopic = true\nomega_sweep = {sweep}\n"), work)
    f = rows[:, 1:3] * 2 * np.pi / w1
    # at Omega = 0 the job solves the still beam: both bending planes at omega_1
    return {"omega_over_omega1": fracs, "backward": f[:, 0].tolist(), "forward": f[:, 1].tolist(),
            "backward_difference": (f[1:, 0] - (1 - np.array(fracs[1:]))).tolist(),
            "forward_difference": (f[1:, 1] - (1 + np.array(fracs[1:]))).tolist()}


def figures(summary: dict, out_dir: str) -> None:
    """Campbell diagrams: one marker per solved speed, references as hollow markers or exact
    lines, and the difference from the reference at every speed."""
    from onedfem.post_processing.graphical_visualisers import figure_style as fs

    fs.use()
    import matplotlib.pyplot as plt

    h, wh = summary["hodges"], summary["whirl"]
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.0), gridspec_kw={"width_ratios": [1.2, 1.2, 1]})
    ax = axes[0]
    b = np.array(h["omega_bar"], dtype=float)
    fs.reference(ax, b, h["hodges_flap"], label="flap, Hodges & Rutkowski (1981)", ms=7)
    fs.solves(ax, b, h["flap"], color=fs.colour("Timoshenko"), label="flap, 1d-fem")
    ax.plot(b, h["lag_exact"], "^", ls="none", mfc="white", mec=fs.INK, ms=7, mew=0.9, label="lag, √(flap² − Ω̄²)")
    fs.solves(ax, b, h["lag"], color=fs.colour("GBT"), label="lag, 1d-fem")
    ax.set(xlabel="rotation speed Ω̄ = Ω √(m L⁴ / EI)", ylabel="first frequency ω̄ (rotating frame)")
    fs.panel_title(ax, "Cantilever rotating about its root", f"{N_EL} Euler-Bernoulli elements, one eigen job per speed")
    ax.legend(loc="upper left")
    ax = axes[1]
    r = np.array(wh["omega_over_omega1"])
    rr = np.linspace(0, r.max(), 2)
    fs.analytic(ax, rr, 1 - rr, label="exact: ω₁ − Ω and ω₁ + Ω")
    fs.analytic(ax, rr, 1 + rr, label=None)
    fs.solves(ax, r, wh["backward"], color=fs.colour("GBT"), label="backward whirl, 1d-fem")
    fs.solves(ax, r, wh["forward"], color=fs.colour("Timoshenko"), label="forward whirl, 1d-fem")
    ax.set(xlabel="spin Ω / ω₁", ylabel="frequency / ω₁ (rotating frame)")
    fs.panel_title(ax, "Shaft spinning about its own axis", "Coriolis (gyroscopic) terms")
    ax.legend(loc="upper left")
    ax = axes[2]
    bd = b[-len(h["flap_difference"]):]
    fs.solves(ax, bd / b.max(), 100 * np.array(h["flap_difference"]), color=fs.colour("Timoshenko"),
              label="flap, against Hodges")
    fs.solves(ax, bd / b.max(), 100 * np.array(h["lag_difference"]), color=fs.colour("GBT"), label="lag, against exact")
    rd = r[-len(wh["backward_difference"]):]
    fs.solves(ax, rd / r.max(), 100 * np.array(wh["backward_difference"]), color=fs.colour("Benscoter"), marker="s",
              label="whirl, backward (of ω₁)")
    fs.solves(ax, rd / r.max(), 100 * np.array(wh["forward_difference"]), color=fs.colour("2D sections"), marker="s",
              label="whirl, forward (of ω₁)")
    ax.axhline(0, color=fs.MUTED, lw=0.6)
    ax.set(xlabel="speed / largest speed solved", ylabel="difference (%)")
    fs.panel_title(ax, "Difference from the reference", "flap and lag within 0.02 %, whirl within 0.03 % of ω₁")
    ax.legend(loc="lower left")
    fig.tight_layout()
    fs.save(fig, os.path.join(out_dir, "campbell.png"))


def main(plot: bool = True, work_dir: str = HERE, out_dir: str = HERE) -> dict:
    logging.disable(logging.WARNING)
    try:
        summary = {"hodges": hodges(work_dir), "whirl": whirl(work_dir)}
    finally:
        logging.disable(logging.NOTSET)
    with open(os.path.join(out_dir, "summary.json"), "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=1)
    if plot:
        figures(summary, out_dir)
    return summary


if __name__ == "__main__":
    s = main()
    print("flap", np.round(np.array(s["hodges"]["flap_difference"]) * 100, 3))
    print("lag", np.round(np.array(s["hodges"]["lag_difference"]) * 100, 3))
    print("whirl", np.round(s["whirl"]["backward_difference"], 5), np.round(s["whirl"]["forward_difference"], 5))

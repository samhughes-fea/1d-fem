# examples/thin_walled/layup_sensitivity/run.py
"""Design sensitivities of a composite box spar to its ply angles and thicknesses.

Run from the repository root::

    python examples/thin_walled/layup_sensitivity/run.py

A 6 m CFRP box spar (0.6 x 0.3 m) with an off-axis ply pair in each cap. The bottom cap's
pair is at ``+theta`` and the top cap's at ``-theta`` in wall axes, which is the same fibre
direction in the global frame on both caps, so the spar twists when it bends (bend-twist
coupling, as in an aeroelastically tailored blade).

1. **Sensitivities.** With ``layup_sensitivity.txt`` in a static and an eigen job, 1d-fem
   writes the derivatives of the tip flap deflection, the tip twist, the compliance and the
   first four frequencies with respect to every chosen ply angle and thickness
   (``layup_sensitivity.csv``).
2. **Check.** Each derivative against central differences of two re-solved jobs.
3. **Design.** Newton iterations on ``theta`` (the two caps linked, ``theta_top = -theta``)
   to reach a target twist per unit flap deflection, each step a job with its sensitivities.

Writes ``jobs/`` and ``results/`` (not tracked), and ``summary.json``, ``sensitivities.png``
and ``design.png`` next to this file.
"""

from __future__ import annotations

import csv
import json
import logging
import os
import sys
from typing import Dict, List, Optional

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
for p in (REPO, os.path.dirname(HERE)):
    if p not in sys.path:
        sys.path.insert(0, p)

from common import CantileverSpec, Material, SectionSpec, results_dir, write_job, write_supports  # noqa: E402

from onedfem.pre_processing.parsing.layup_sensitivity_parser import Response, format_layup_sensitivity_file  # noqa: E402
from anisotropic_section_model.thin_walled.layup_sensitivity import LayupParameter  # noqa: E402

CFRP = Material(E1=135e9, E2=10e9, G12=5e9, nu12=0.3, G13=5e9, G23=3.8e9, rho=1550.0)
B, H, T_PLY = 0.6, 0.3, 0.5e-3
L, N_EL = 6.0, 12
P_FLAP = 20e3                                  # tip force along z, N
THETA0 = 20.0                                  # off-axis pair, degrees
TARGET = -0.18                                 # target tip twist per tip flap deflection, deg per 100 mm
THETA_START = 5.0                              # first design iterate, degrees
MAX_STEP = 10.0                                # largest Newton step, degrees
SWEEP = tuple(range(0, 45, 5))                 # angles of the jobs plotted behind the design loop
BOTTOM, WEB, TOP = 1, 2, 3                     # laminate ids
CAP_PLIES = ("45", "-45", "theta", "theta", "0", "0", "0", "0", "-45", "45")
WEB_PLIES = (45.0, -45.0, -45.0, 45.0)
# the chosen parameters: the off-axis pair and a 0-degree ply of each cap, and the webs
PARAMS = [LayupParameter(BOTTOM, 2, "angle"), LayupParameter(BOTTOM, 3, "angle"),
          LayupParameter(TOP, 2, "angle"), LayupParameter(TOP, 3, "angle"),
          LayupParameter(WEB, 0, "angle"), LayupParameter(WEB, 1, "angle"),
          LayupParameter(BOTTOM, 4, "thickness"), LayupParameter(TOP, 4, "thickness"),
          LayupParameter(WEB, 0, "thickness")]
STATIC = [Response("displacement", N_EL, "UZ"), Response("displacement", N_EL, "RX"), Response("compliance")]
EIGEN = [Response("frequency", mode=k) for k in (1, 2, 3, 4)]


def cap(theta: float, sign: float, t0: float = T_PLY) -> List[tuple]:
    out = []
    for k, a in enumerate(CAP_PLIES):
        ang = sign * theta if a == "theta" else float(a)
        out.append((1, ang, t0))
    return out


def section(theta: float = THETA0, changes: Optional[Dict[LayupParameter, float]] = None) -> SectionSpec:
    lams = {BOTTOM: cap(theta, +1.0), TOP: cap(theta, -1.0), WEB: [(1, a, T_PLY) for a in WEB_PLIES]}
    for prm, dv in (changes or {}).items():
        mid, ang, t = lams[prm.laminate_id][prm.ply]
        lams[prm.laminate_id][prm.ply] = (mid, ang + dv, t) if prm.kind == "angle" else (mid, ang, t + dv)
    return SectionSpec(
        points={1: (-B / 2, -H / 2), 2: (B / 2, -H / 2), 3: (B / 2, H / 2), 4: (-B / 2, H / 2)},
        walls=[(1, 1, 2, BOTTOM), (2, 2, 3, WEB), (3, 3, 4, TOP), (4, 4, 1, WEB)],
        materials={1: CFRP}, laminates=lams, n_sub=2)


def run_job(work: str, name: str, sec: SectionSpec, kind: str, params=PARAMS) -> Dict[str, Dict[str, dict]]:
    """Solve one job; ``{response: {parameter: row}}`` from ``layup_sensitivity.csv``."""
    from onedfem.workflow_orchestrator.run_job import process_job

    spec = CantileverSpec(L, N_EL, "LinearTimoshenkoBeamElement3D", sec, [0.0, 0.0, P_FLAP, 0.0, 0.0, 0.0])
    job = os.path.join(work, "jobs", name)
    write_job(job, spec)
    write_supports(job, spec)
    if kind == "eigen":
        with open(os.path.join(job, "simulation_settings.txt"), "w", encoding="utf-8") as f:
            f.write("[Simulation]\n[Type]\neigen\n\n[Solver]\ntype = direct\n\n[Eigen]\nnum_modes = 6\n")
    with open(os.path.join(job, "layup_sensitivity.txt"), "w", encoding="utf-8") as f:
        f.write(format_layup_sensitivity_file(params, STATIC if kind == "static" else EIGEN))
    res = results_dir(os.path.join(work, "results", name))
    process_job(job, res, {}, {}, force_serial=True, max_processes_per_job=1)
    out: Dict[str, Dict[str, dict]] = {}
    with open(os.path.join(res, "secondary_results", "sensitivity", "layup_sensitivity.csv"), newline="",
              encoding="utf-8") as f:
        for r in csv.DictReader(f):
            out.setdefault(r["response"], {})[r["parameter"]] = r
    return out


def value(rows, response: str) -> float:
    return float(next(iter(rows[response].values()))["response_value"])


def check(work: str, base: Dict[str, Dict[str, dict]], kind: str) -> List[dict]:
    """Every derivative against central differences of two re-solved jobs."""
    out = []
    for prm in PARAMS:
        h = 1e-3 if prm.kind == "angle" else 1e-8
        rp = run_job(work, f"{kind}_{prm.label}_p", section(changes={prm: +h}), kind, [prm])
        rm = run_job(work, f"{kind}_{prm.label}_m", section(changes={prm: -h}), kind, [prm])
        for resp, by_param in base.items():
            d = float(by_param[prm.label]["derivative"])
            fd = (value(rp, resp) - value(rm, resp)) / (2 * h)
            out.append({"parameter": prm.label, "response": resp, "adjoint": d, "central_difference": fd,
                        "relative_error": abs(d - fd) / max(abs(fd), 1e-300)})
    return out


def ratio_and_gradient(rows) -> tuple:
    """Twist per flap deflection (deg per 100 mm) and its derivative per degree of theta (linked caps)."""
    uz, rx = value(rows, f"UZ@{N_EL}"), value(rows, f"RX@{N_EL}")
    g = {}
    for resp in (f"UZ@{N_EL}", f"RX@{N_EL}"):
        d = rows[resp]
        # theta enters both plies of the pair, +theta on the bottom cap and -theta on the top cap
        g[resp] = sum(float(d[LayupParameter(BOTTOM, k, "angle").label]["derivative"])
                      - float(d[LayupParameter(TOP, k, "angle").label]["derivative"]) for k in (2, 3))
    c = np.degrees(1.0) * 0.1                   # rad per m -> deg per 100 mm
    r = c * rx / uz
    dr = c * (g[f"RX@{N_EL}"] * uz - rx * g[f"UZ@{N_EL}"]) / uz**2
    return r, dr


def design(work: str, theta0: float = THETA_START, iterations: int = 6) -> List[dict]:
    hist = []
    theta = theta0
    for it in range(iterations):
        rows = run_job(work, f"design_{it}", section(theta), "static")
        r, dr = ratio_and_gradient(rows)
        hist.append({"iteration": it, "theta": theta, "twist_per_flap": r, "gradient": dr})
        if abs(r - TARGET) < 1e-6 * abs(TARGET):
            break
        theta = theta + float(np.clip((TARGET - r) / dr, -MAX_STEP, MAX_STEP))
    return hist


def main(plot: bool = True, work_dir: str = HERE, out_dir: str = HERE, verify: bool = True) -> dict:
    logging.disable(logging.WARNING)
    try:
        st = run_job(work_dir, "static", section(), "static")
        ei = run_job(work_dir, "eigen", section(), "eigen")
        checks = (check(work_dir, st, "static") + check(work_dir, ei, "eigen")) if verify else []
        hist = design(work_dir)
        sweep = []
        for th in SWEEP:
            r, dr = ratio_and_gradient(run_job(work_dir, f"sweep_{th:g}", section(th), "static"))
            sweep.append({"theta": th, "twist_per_flap": r, "gradient": dr})
    finally:
        logging.disable(logging.NOTSET)
    table = {}
    for rows in (st, ei):
        for resp, by in rows.items():
            table[resp] = {p: {"value": float(r["response_value"]), "derivative": float(r["derivative"]),
                               "relative_derivative": float(r["relative_derivative"]) if r["relative_derivative"]
                               else None} for p, r in by.items()}
    result = {"theta0": THETA0, "target_twist_per_flap": TARGET, "theta_start": THETA_START, "sensitivities": table, "check": checks,
              "max_relative_error": max((c["relative_error"] for c in checks), default=None), "design": hist, "sweep": sweep}
    with open(os.path.join(out_dir, "summary.json"), "w", encoding="utf-8") as f:
        json.dump(result, f, indent=1)
    if plot:
        figures(result, out_dir)
    return result


def figures(result: dict, out_dir: str) -> None:
    """Sensitivity bars, the adjoint-against-central-difference check (one point per derivative),
    the design loop (one marker per job) and the section with its coupled layup."""
    from onedfem.post_processing.graphical_visualisers import figure_style as fs

    from common import build_analysis, plot_section, ply_materials

    fs.use()
    import matplotlib.pyplot as plt

    sens = result["sensitivities"]
    angles = [p.label for p in PARAMS if p.kind == "angle"]
    thick = [p.label for p in PARAMS if p.kind == "thickness"]
    names = {f"UZ@{N_EL}": "tip flap deflection", f"RX@{N_EL}": "tip twist", "compliance": "compliance",
             "f1": "f1", "f2": "f2", "f3": "f3", "f4": "f4"}
    nice = {"L1": "bottom cap", "L2": "web", "L3": "top cap"}

    def lab(s):
        lid, ply, _ = s.split(".")
        return f"{nice[lid]} ply {ply[1:]}"

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    resp = [f"UZ@{N_EL}", f"RX@{N_EL}", "f1", "f2", "f3", "f4"]
    for ax, plist, title, key in ((axes[0], angles, "Per degree of ply angle", "derivative"),
                                  (axes[1], thick, "Logarithmic thickness derivative", "relative_derivative")):
        x = np.arange(len(plist))
        w = 0.8 / len(resp)
        for i, r in enumerate(resp):
            vals = []
            for p in plist:
                d = sens[r][p]
                vals.append(100 * d["derivative"] / d["value"] if key == "derivative" else d["relative_derivative"])
            ax.bar(x + (i - len(resp) / 2 + 0.5) * w, vals, w, color=fs.colour("", i), label=names[r], lw=0)
        ax.axhline(0, color=fs.MUTED, lw=0.6)
        ax.set_xticks(x)
        ax.set_xticklabels([lab(p) for p in plist], rotation=20, ha="right")
        ax.grid(axis="x", visible=False)
        fs.panel_title(ax, title, "% of the response per degree" if key == "derivative" else "(t / J) dJ/dt")
    axes[0].legend(ncol=2, loc="lower left")
    fig.suptitle("Box spar: adjoint sensitivities of the tip response and the first four frequencies to the layup",
                 x=0.01, ha="left")
    fig.tight_layout()
    fs.save(fig, os.path.join(out_dir, "sensitivities.png"))

    fig, ax = plt.subplots(1, 2, figsize=(10.5, 3.9))
    chk = result["check"]
    a_ = np.array([abs(c["adjoint"]) for c in chk])
    f_ = np.array([abs(c["central_difference"]) for c in chk])
    e_ = np.array([max(c["relative_error"], 1e-16) for c in chk])
    lo, hi = min(a_.min(), f_.min()) / 3, max(a_.max(), f_.max()) * 3
    ax[0].plot([lo, hi], [lo, hi], color=fs.MUTED, lw=0.8, zorder=0)
    fs.solves(ax[0], f_, a_, color=fs.colour("Timoshenko"), label=f"{len(chk)} derivatives")
    ax[0].set(xscale="log", yscale="log", xlabel="|central difference|", ylabel="|adjoint|")
    fs.panel_title(ax[0], "Adjoint against central differences", "each derivative from one job's "
                                                                  "layup_sensitivity.csv")
    ax[0].legend(loc="upper left")
    resp_of = sorted({c["response"] for c in chk})
    for i, r in enumerate(resp_of):
        sel = [k for k, c in enumerate(chk) if c["response"] == r]
        fs.solves(ax[1], [i] * len(sel), e_[sel], color=fs.colour("", i), ms=4.5)
    ax[1].set_xticks(range(len(resp_of)))
    ax[1].set_xticklabels([names.get(r, r) for r in resp_of], rotation=20, ha="right")
    ax[1].set(yscale="log", ylabel="relative difference")
    ax[1].axhline(result["max_relative_error"], color=fs.MUTED, lw=0.6, ls="--")
    fs.panel_title(ax[1], "Relative difference", f"largest {result['max_relative_error']:.1e}, the step error of "
                                                 "the differences")
    fig.tight_layout()
    fs.save(fig, os.path.join(out_dir, "check.png"))

    hist, sweep = result["design"], result["sweep"]
    fig, ax = plt.subplots(figsize=(6.4, 4.0))
    fs.solves(ax, [h["theta"] for h in sweep], [h["twist_per_flap"] for h in sweep], color=fs.GRID, ms=6,
              label="1d-fem job (sweep at 5° steps)")
    for h in hist:
        th = np.array([h["theta"] - 3.0, h["theta"] + 3.0])
        ax.plot(th, h["twist_per_flap"] + h["gradient"] * (th - h["theta"]), color=fs.colour("Timoshenko"), lw=0.8)
    fs.solves(ax, [h["theta"] for h in hist], [h["twist_per_flap"] for h in hist], color=fs.colour("Timoshenko"),
              label="Newton iterate, with its adjoint tangent")
    shown = []
    for h in hist:                              # one label for iterates that sit on top of each other
        if [k for k, t in shown if abs(t - h["theta"]) < 0.3]:
            continue
        same = [g["iteration"] for g in hist if abs(g["theta"] - h["theta"]) < 0.3]
        shown.append((h["iteration"], h["theta"]))
        ax.annotate(", ".join(map(str, same)), (h["theta"], h["twist_per_flap"]), textcoords="offset points",
                    xytext=(5, 4), fontsize=7, color=fs.colour("Timoshenko"))
    ax.axhline(result["target_twist_per_flap"], color=fs.colour("GBT"), lw=1, label="target")
    ax.set(xlabel="off-axis ply angle θ (bottom cap +θ, top cap −θ), deg",
           ylabel="tip twist per flap deflection (deg / 100 mm)")
    fs.panel_title(ax, "Design loop on θ", "each step is one job; gradient from the job's own sensitivities")
    ax.legend(loc="upper right")
    fig.tight_layout()
    fs.save(fig, os.path.join(out_dir, "design.png"))

    sec = section()
    an = build_analysis(sec)
    plot_section(an, ply_materials(sec, {k: "CFRP" for k in sec.materials}), os.path.join(out_dir, "section.png"),
                 "Box spar with bend-twist coupled caps",
                 f"caps [45, −45, θ, θ, 0₄, −45, 45] at θ = {THETA0:g}° (+θ bottom, −θ top in wall axes)",
                 (B / 2 - 0.03, H / 2 - 0.025, 0.045, 0.035), "Rear web at the top cap", figsize=(10, 4.6),
                 width_ratios=(1.7, 1))


if __name__ == "__main__":
    out = main()
    print(f"largest adjoint / central-difference difference: {out['max_relative_error']:.1e}")
    for h in out["design"]:
        print(f"iteration {h['iteration']}: theta {h['theta']:.4f} deg, twist per flap {h['twist_per_flap']:.5f}")

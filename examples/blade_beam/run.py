"""A blade beam from an implicit-geometry-kernel loft, against a shell model of the same blade.

Run from the repository root (needs implicit-geometry-kernel, ``1d-fem[geometry]``)::

    python examples/blade_beam/run.py

The outboard 20.2 m of the IEA Wind 15 MW blade (r = 64.5 m to 84.7 m; chord 4.07 m to
3.27 m, twist 1.5 deg to -0.4 deg, ``examples/iea15mw/segment.py``) is lofted once by
implicit-geometry-kernel. ``blade_beam_job`` cuts a section at every element mid-point and
writes the 1d-fem job: the reference axis, the element frames, the twist and pretwist all
come from the loft. Three beam models are compared with a composite-laminate-model shell
blade meshed from the same loft (``data/shell_reference.json``, ``make_shell_reference.py``):

1. **Static:** flap and edge tip forces and a tip torque, for Timoshenko and Benscoter
   elements on thin-walled sections and Timoshenko elements on 2D sections, tapered and
   prismatic.
2. **Modal:** the first frequencies of the cantilever, Timoshenko and Benscoter.
3. **Nonlinear:** corotational Timoshenko elements under a growing flapwise tip force,
   against CalculiX ``NLGEOM`` shells.

Writes ``jobs/`` and ``results/`` (not tracked), and the figures and ``summary.json`` next
to this file.
"""

from __future__ import annotations

import csv
import glob
import json
import logging
import os
import shutil
import sys
import time
from typing import Dict, Optional

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
for p in (REPO, os.path.join(REPO, "examples", "iea15mw")):
    if p not in sys.path:
        sys.path.insert(0, p)
import segment as G  # noqa: E402
import stations as S  # noqa: E402

TW, BENSCOTER, TS_2D = "timoshenko_thin_walled", "benscoter_thin_walled", "timoshenko_2d"
MODELS = {TW: ("thin_walled", "LinearTimoshenkoBeamElement3D"),
          BENSCOTER: ("thin_walled", "LinearBenscoterBeamElement3D"),
          TS_2D: ("2d", "LinearTimoshenkoBeamElement3D")}
NAMES = {TW: "Timoshenko, thin-walled", BENSCOTER: "Benscoter, thin-walled", TS_2D: "Timoshenko, 2D sections",
         "shell": "CLM shell"}
COLOURS = {TW: "#2a78d6", BENSCOTER: "#eb6834", TS_2D: "#1baf7a", "shell": "#0b0b0b"}
# blade frame: x chordwise, y flapwise, z spanwise; loads [F_x, F_y, F_z, M_x, M_y, M_z]
STATIC = {"flap": ([0, 1e4, 0, 0, 0, 0], 1), "edge": ([1e4, 0, 0, 0, 0, 0], 0), "torque": ([0, 0, 0, 0, 0, 1e4], 5)}
NONLINEAR_FRACTIONS = (0.1, 0.25, 0.45)
EIGEN = "[Simulation]\n[Type]\neigen\n\n[Eigen]\nenabled = true\nnum_modes = 6\n"


def reference() -> dict:
    with open(os.path.join(HERE, "data", "shell_reference.json"), encoding="utf-8") as f:
        return json.load(f)


def _job(name, work, model, n, prismatic=False, **kw):
    from pre_processing.mesh_library.blade_beam_job import blade_beam_job

    kind, et = MODELS[model]
    job = blade_beam_job(os.path.join(work, "jobs", name), G.geometry(prismatic), S.laminates(20), n,
                         kind=kind, element_type=kw.pop("element_type", et), **kw)
    # stresses at one station only (the root element's mid-point): the stress files of a
    # blade job with hundreds of walls per section are large and slow to write
    with open(os.path.join(job.job_dir, "stress_stations.txt"), "w", encoding="utf-8") as f:
        f.write(f"[Stress_Station]\n[x]\n{0.5 * job.nodes[1, 0]:.10g}\n")
    return job


def _run(job, work) -> str:
    from examples.thin_walled.common import results_dir
    from workflow_orchestrator.run_job import process_job

    res = results_dir(os.path.join(work, "results", os.path.basename(job.job_dir)))
    process_job(job.job_dir, res, {}, {}, force_serial=True, max_processes_per_job=1)
    return res


def _tip(res, job) -> np.ndarray:
    """Tip translations and rotations in the blade frame."""
    from pre_processing.mesh_library.blade_beam_job import job_to_blade

    with open(os.path.join(res, "primary_results", "global", "U_global.csv"), encoding="utf-8") as f:
        U = np.array([float(r[1]) for r in list(csv.reader(f))[1:]])
    U = U.reshape(job.nodes.shape[0], -1)
    return job_to_blade(U[job.node_at(G.LENGTH)][:6])


def _modes(res, job):
    """Frequencies (Hz) and the blade-frame tip motion [u_x, u_y, u_z, twist] of each mode."""
    from pre_processing.mesh_library.blade_beam_job import job_to_blade

    f = sorted(glob.glob(os.path.join(res, "**", "*_frequencies.txt"), recursive=True))[0]
    freq = np.atleast_1d(np.loadtxt(f))
    Phi = np.loadtxt(f.replace("_frequencies.txt", "_mode_shapes.txt"))
    Phi = Phi.reshape(Phi.shape[0], -1)
    nd = Phi.shape[0] // job.nodes.shape[0]
    tips = []
    for k in range(Phi.shape[1]):
        u = Phi[:, k].reshape(job.nodes.shape[0], nd)[job.node_at(G.LENGTH)][:6]
        b = job_to_blade(u)
        tips.append([b[0], b[1], b[2], b[5]])
    return freq.tolist(), tips


def classify(tip, chord: float) -> str:
    """flap, edge or torsion: the largest of the tip's flap and edge motion and its twist
    times half a chord."""
    ux, uy, _, tw = tip
    k = int(np.argmax([abs(uy), abs(ux), abs(tw) * 0.5 * chord]))
    return ("flap", "edge", "torsion")[k]


def run_static(work, n) -> dict:
    out = {}
    for model in MODELS:
        out[model] = {}
        for case, (load, _) in STATIC.items():
            from pre_processing.mesh_library.blade_beam_job import blade_to_job

            job = _job(f"{model}_{case}", work, model, n, loads={G.LENGTH: blade_to_job(load)})
            out[model][case] = _tip(_run(job, work), job).tolist()
    out["untapered"] = {}
    for case, (load, _) in STATIC.items():
        from pre_processing.mesh_library.blade_beam_job import blade_to_job

        job = _job(f"untapered_{case}", work, TW, n, taper=False, loads={G.LENGTH: blade_to_job(load)})
        out["untapered"][case] = _tip(_run(job, work), job).tolist()
    out["prismatic"] = {}
    for case, (load, _) in STATIC.items():
        from pre_processing.mesh_library.blade_beam_job import blade_to_job

        job = _job(f"prismatic_{case}", work, TW, n, prismatic=True, loads={G.LENGTH: blade_to_job(load)})
        out["prismatic"][case] = _tip(_run(job, work), job).tolist()
    return out


def run_modal(work, n) -> dict:
    out = {}
    for model in (TW, BENSCOTER):
        job = _job(f"{model}_modal", work, model, n, settings=EIGEN)
        freq, tips = _modes(_run(job, work), job)
        out[model] = {"frequencies_hz": freq, "tips": tips, "kinds": [classify(t, G.TIP["chord"]) for t in tips]}
    return out


def run_nonlinear(work, n, P: float) -> dict:
    from pre_processing.mesh_library.blade_beam_job import blade_to_job

    rows = []
    for frac in NONLINEAR_FRACTIONS:
        load = {G.LENGTH: blade_to_job([0, frac * P, 0, 0, 0, 0])}
        job = _job(f"nonlinear_{frac:.2f}", work, TW, n, loads=load, element_type="CorotationalTimoshenkoBeamElement3D")
        rows.append([frac, *_tip(_run(job, work), job).tolist()])
    return {"rows": rows, "P": P}


N_COMPARED = 4   # beam modes compared; the shell's higher modes include panel modes


def compare(static, modal, nonlinear, ref) -> dict:
    s_ref = ref["static"]["tip"]
    out = {"static": {}, "static_prismatic": {}, "modal": {}, "nonlinear": {}}
    for model in MODELS:
        out["static"][model] = {case: static[model][case][k] / s_ref[case][k] - 1.0 for case, (_, k) in STATIC.items()}
    out["static_untapered"] = {TW: {case: static["untapered"][case][k] / s_ref[case][k] - 1.0
                                    for case, (_, k) in STATIC.items()}}
    sp = ref["static_prismatic"]["tip"]
    out["static_prismatic"][TW] = {case: static["prismatic"][case][k] / sp[case][k] - 1.0 for case, (_, k) in STATIC.items()}
    # shell modes, classified the same way; beam modes matched in order within each kind
    shell_kinds = [classify(t[1:], G.TIP["chord"]) for t in ref["modal"]["tip_series"]]
    sf = ref["modal"]["frequencies_hz"]
    for model, m in modal.items():
        pairs = []
        used = set()
        for j, (f, kind) in enumerate(zip(m["frequencies_hz"][:N_COMPARED], m["kinds"])):
            cand = [i for i, k in enumerate(shell_kinds) if k == kind and i not in used]
            if not cand:
                continue
            i = cand[0]
            used.add(i)
            pairs.append({"beam_mode": j + 1, "shell_mode": i + 1, "kind": kind, "beam_hz": f, "shell_hz": sf[i],
                          "difference": f / sf[i] - 1.0})
        out["modal"][model] = pairs
    out["modal"]["shell_kinds"] = shell_kinds
    shell = {round(r[0], 4): r for r in ref["nonlinear"]["tip_series"]}
    lin = static[TW]["flap"][1] / 1e4      # beam flap compliance
    lin_s = s_ref["flap"][1] / 1e4
    for r in nonlinear["rows"]:
        frac = r[0]
        sh = shell.get(round(frac, 4))
        P = frac * nonlinear["P"]
        row = {"load": P, "beam_flap": r[2], "beam_shortening": r[3], "beam_linear": lin * P,
               "beam_nonlinear_ratio": r[2] / (lin * P)}
        if sh is not None:
            row.update({"shell_flap": sh[2], "shell_shortening": sh[3], "shell_linear": lin_s * P,
                        "shell_nonlinear_ratio": sh[2] / (lin_s * P), "flap_difference": r[2] / sh[2] - 1.0})
        out["nonlinear"][f"{frac:.2f}"] = row
    return out


# ---------------------------------------------------------------------------
def _style(ax):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.tick_params(labelsize=7, colors="#52514e")


def figures(result, ref, out_dir):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    # 1. the loft and its cuts
    geo = G.geometry()
    fig = plt.figure(figsize=(9, 3.6))
    ax = fig.add_axes((-0.25, -0.45, 1.5, 1.8), projection="3d")
    for z in np.linspace(0, G.LENGTH, 9):
        cut = geo.section_at(float(z))
        for w in cut.walls:
            xy = np.asarray(w.xy)
            ax.plot(np.full(len(xy), z), xy[:, 0], xy[:, 1], color="#2a78d6" if w.kind != "skin" else "#8b8a83", lw=0.6)
    ax.set_box_aspect((G.LENGTH, 4.5, 1.2), zoom=1.0)
    ax.set_xlabel("span z, m", fontsize=7)
    ax.set_yticks([])
    ax.set_zticks([])
    ax.view_init(elev=22, azim=-62)
    fig.text(0.01, 0.95, "IEA 15 MW outboard segment lofted by implicit-geometry-kernel (caps and webs blue)",
             fontsize=9)
    fig.savefig(os.path.join(out_dir, "loft.png"), dpi=150)
    plt.close(fig)

    # 2. static differences
    cmp_ = result["comparison"]
    fig, ax = plt.subplots(figsize=(7.5, 3.4))
    cases = list(STATIC)
    x = np.arange(len(cases))
    series = [(m, cmp_["static"][m], NAMES[m]) for m in MODELS]
    series.insert(0, ("untapered", cmp_["static_untapered"][TW], "Timoshenko, thin-walled, no taper coupling"))
    series.append(("prismatic", cmp_["static_prismatic"][TW], "Timoshenko, thin-walled, prismatic segment"))
    w = 0.8 / len(series)
    cols = dict(COLOURS, prismatic="#6250d6", untapered="#a3a29b")
    for j, (key, vals, lab) in enumerate(series):
        xs = x - 0.4 + w * (j + 0.5)
        v = [100 * vals[c] for c in cases]
        ax.vlines(xs, 0, v, color=cols[key], lw=2)
        ax.plot(xs, v, "o", color=cols[key], ms=5, mec="#fcfcfb", label=lab)
    ax.fill_between([-0.5, len(cases) - 0.5], -3, 3, color="#e9e8e3", zorder=0, lw=0)
    ax.axhline(0, color="#8b8a83", lw=0.8)
    ax.set_xticks(x, ["flap tip deflection", "edge tip deflection", "tip twist"])
    ax.set_ylabel("beam against shell, %", fontsize=8)
    ax.set_ylim(-6, 8)
    ax.legend(fontsize=7, frameon=False, loc="upper center", ncol=2)
    ax.set_title("Static tip response (grey: 3 % target)", fontsize=9, loc="left")
    _style(ax)
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, "static.png"), dpi=150)
    plt.close(fig)

    # 3. frequencies
    fig, ax = plt.subplots(figsize=(7.5, 3.4))
    sf = ref["modal"]["frequencies_hz"]
    kinds = cmp_["modal"]["shell_kinds"]
    for i, (f, k) in enumerate(zip(sf, kinds)):
        ax.plot([0.6, 1.4], [f, f], color="#0b0b0b", lw=1.2)
        ax.text(1.45, f, f"{k}", fontsize=6, va="center", color="#52514e")
    for j, model in enumerate((TW, BENSCOTER)):
        m = result["modal"][model]
        for f, k in zip(m["frequencies_hz"], m["kinds"]):
            ax.plot([1.9 + j, 2.7 + j], [f, f], color=COLOURS[model], lw=1.2)
            ax.text(2.75 + j, f, k, fontsize=6, va="center", color="#52514e")
    ax.set_xticks([1, 2.3, 3.3], ["CLM shell\n(CalculiX)", "Timoshenko", "Benscoter"])
    ax.set_ylabel("frequency, Hz", fontsize=8)
    ax.set_xlim(0.4, 4.2)
    ax.set_title("Cantilever frequencies", fontsize=9, loc="left")
    _style(ax)
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, "modes.png"), dpi=150)
    plt.close(fig)

    # 4. nonlinear tip path
    fig, axes = plt.subplots(1, 2, figsize=(9, 3.4))
    sh = np.array([r[:4] for r in ref["nonlinear"]["tip_series"]], dtype=float)
    P = ref["nonlinear"]["P"]
    lin = result["static"][TW]["flap"][1] / 1e4
    lin_s = ref["static"]["tip"]["flap"][1] / 1e4
    nl = np.array(result["nonlinear"]["rows"], dtype=float)
    ax = axes[0]
    ax.plot(sh[:, 0] * P / 1e6, sh[:, 2], "s", color=COLOURS["shell"], ms=4, mfc="none",
            label="CLM shell, CalculiX NLGEOM increments")
    ax.plot(nl[:, 0] * P / 1e6, nl[:, 2], "o", color=COLOURS[TW], ms=5, label="corotational Timoshenko load steps")
    f = np.linspace(0, 0.5, 2)
    ax.plot(f * P / 1e6, lin * f * P, "--", color=COLOURS[TW], lw=0.8, label="linear Timoshenko (exactly linear in load)")
    ax.plot(f * P / 1e6, lin_s * f * P, "--", color=COLOURS["shell"], lw=0.8, label="linear shell (exactly linear in load)")
    ax.set_xlabel("flapwise tip force, MN", fontsize=8)
    ax.set_ylabel("flap tip deflection, m", fontsize=8)
    ax.legend(fontsize=7, frameon=False)
    _style(ax)
    ax = axes[1]
    ax.plot(sh[:, 2], -sh[:, 3], "s", color=COLOURS["shell"], ms=4, mfc="none", label="shell increments")
    ax.plot(nl[:, 2], -nl[:, 3], "o", color=COLOURS[TW], ms=5, label="corotational beam load steps")
    ax.set_xlabel("flap tip deflection, m", fontsize=8)
    ax.set_ylabel("tip shortening, m", fontsize=8)
    _style(ax)
    fig.suptitle("Large flapwise deflection", fontsize=9, x=0.01, ha="left")
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, "nonlinear.png"), dpi=150)
    plt.close(fig)


def main(n_elements: int = 20, plot: bool = True, work_dir: Optional[str] = None, parts=("static", "modal", "nonlinear"),
         out_dir: str = HERE) -> dict:
    """Run the example (see the module docstring); logging below warnings is off while it runs."""
    logging.disable(logging.WARNING)
    try:
        return _main(n_elements=n_elements, plot=plot, work_dir=work_dir, parts=parts, out_dir=out_dir)
    finally:
        logging.disable(logging.NOTSET)


def _main(n_elements: int = 20, plot: bool = True, work_dir: Optional[str] = None, parts=("static", "modal", "nonlinear"),
          out_dir: str = HERE) -> dict:
    work = work_dir or HERE
    for d in ("jobs", "results"):
        shutil.rmtree(os.path.join(work, d), ignore_errors=True)
    ref = reference()
    t0 = time.time()
    result: Dict[str, object] = {"n_elements": n_elements, "length": G.LENGTH}
    result["static"] = run_static(work, n_elements) if "static" in parts else None
    result["modal"] = run_modal(work, n_elements) if "modal" in parts else None
    result["nonlinear"] = run_nonlinear(work, n_elements, ref["nonlinear"]["P"]) if "nonlinear" in parts else None
    if all(result[p] is not None for p in ("static", "modal", "nonlinear")):
        result["comparison"] = compare(result["static"], result["modal"], result["nonlinear"], ref)
        if plot:
            figures(result, ref, out_dir)
    result["time"] = time.time() - t0
    with open(os.path.join(out_dir, "summary.json"), "w", encoding="utf-8") as f:
        json.dump(result, f, indent=1)
    return result


if __name__ == "__main__":
    r = main()
    print(json.dumps(r.get("comparison"), indent=1))

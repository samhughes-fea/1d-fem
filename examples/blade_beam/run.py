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
2. **Modal:** the first frequencies of the cantilever, Timoshenko and Benscoter on
   thin-walled sections, Timoshenko on 2D sections, and GBT on 2D sections with four
   distortional modes, whose inertia brings the torsion modes down to the shell's. A
   second GBT beam has its sandwich cores made stiff in transverse shear, as the shell's
   single brick through the wall makes them (``sandwich_strip_check.py``).
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
# torsion modes with the sections' own distortion: GBT on the 2D sections, against the same
# beam without the modes, and GBT with the sandwich cores stiff in transverse shear (the
# shell's kinematics through the wall, see sandwich_strip_check.py)
GBT_2D, GBT_2D_STIFF = "gbt_2d", "gbt_2d_stiff_core"
GBT_MODES = 4
#: transverse shear moduli of the core material are multiplied by this in GBT_2D_STIFF
STIFF_CORE = 100.0
CORE = "medium_density_foam"
MODAL_MODELS = {TW: {}, BENSCOTER: {}, TS_2D: {},
                GBT_2D: dict(model=TS_2D, element_type="LinearGBTBeamElement3D", gbt_modes=GBT_MODES),
                GBT_2D_STIFF: dict(model=TS_2D, element_type="LinearGBTBeamElement3D", gbt_modes=GBT_MODES,
                                   stiff_core=True)}
NAMES.update({GBT_2D: f"GBT, 2D sections, {GBT_MODES} distortional modes",
              GBT_2D_STIFF: f"GBT, 2D sections, {GBT_MODES} modes, cores stiff in transverse shear"})
COLOURS = {TW: "#2a78d6", BENSCOTER: "#eb6834", TS_2D: "#1baf7a", "shell": "#0b0b0b",
           GBT_2D: "#6250d6", GBT_2D_STIFF: "#a3a29b"}
# blade frame: x chordwise, y flapwise, z spanwise; loads [F_x, F_y, F_z, M_x, M_y, M_z]
STATIC = {"flap": ([0, 1e4, 0, 0, 0, 0], 1), "edge": ([1e4, 0, 0, 0, 0, 0], 0), "torque": ([0, 0, 0, 0, 0, 1e4], 5)}
NONLINEAR_FRACTIONS = (0.1, 0.25, 0.45)
EIGEN = "[Simulation]\n[Type]\neigen\n\n[Eigen]\nenabled = true\nnum_modes = 6\n"


def reference() -> dict:
    with open(os.path.join(HERE, "data", "shell_reference.json"), encoding="utf-8") as f:
        return json.load(f)


def _job(name, work, model, n, prismatic=False, **kw):
    from onedfem.pre_processing.mesh_library.blade_beam_job import blade_beam_job

    kind, et = MODELS[model]
    laminates = _stiff_core_laminates() if kw.pop("stiff_core", False) else S.laminates(20)
    job = blade_beam_job(os.path.join(work, "jobs", name), G.geometry(prismatic), laminates, n,
                         kind=kind, element_type=kw.pop("element_type", et), **kw)
    # stresses at one station only (the root element's mid-point): the stress files of a
    # blade job with hundreds of walls per section are large and slow to write
    with open(os.path.join(job.job_dir, "stress_stations.txt"), "w", encoding="utf-8") as f:
        f.write(f"[Stress_Station]\n[x]\n{0.5 * job.nodes[1, 0]:.10g}\n")
    return job


def _stiff_core_laminates() -> dict:
    """The station laminates with the core's transverse shear moduli times ``STIFF_CORE``."""
    import dataclasses

    from onedfem.pre_processing.section_integration.thin_walled import Ply, WallLaminate

    out = {}
    for k, lam in S.laminates(20).items():
        plies = []
        for p in lam.plies:
            m = p.material
            if p.material == S.ply_materials()[CORE]:
                m = dataclasses.replace(m, G13=STIFF_CORE * m.G13, G23=STIFF_CORE * m.G23)
            plies.append(dataclasses.replace(p, material=m))
        out[k] = WallLaminate(plies)
    return out


def _run(job, work) -> str:
    from examples.thin_walled.common import results_dir
    from onedfem.workflow_orchestrator.run_job import process_job

    res = results_dir(os.path.join(work, "results", os.path.basename(job.job_dir)))
    process_job(job.job_dir, res, {}, {}, force_serial=True, max_processes_per_job=1)
    return res


def _tip(res, job) -> np.ndarray:
    """Tip translations and rotations in the blade frame."""
    from onedfem.pre_processing.mesh_library.blade_beam_job import job_to_blade

    with open(os.path.join(res, "primary_results", "global", "U_global.csv"), encoding="utf-8") as f:
        U = np.array([float(r[1]) for r in list(csv.reader(f))[1:]])
    U = U.reshape(job.nodes.shape[0], -1)
    return job_to_blade(U[job.node_at(G.LENGTH)][:6])


def _modes(res, job):
    """Frequencies (Hz) and the blade-frame tip motion [u_x, u_y, u_z, twist] of each mode."""
    from onedfem.pre_processing.mesh_library.blade_beam_job import job_to_blade

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
            from onedfem.pre_processing.mesh_library.blade_beam_job import blade_to_job

            job = _job(f"{model}_{case}", work, model, n, loads={G.LENGTH: blade_to_job(load)})
            out[model][case] = _tip(_run(job, work), job).tolist()
    out["untapered"] = {}
    for case, (load, _) in STATIC.items():
        from onedfem.pre_processing.mesh_library.blade_beam_job import blade_to_job

        job = _job(f"untapered_{case}", work, TW, n, taper=False, loads={G.LENGTH: blade_to_job(load)})
        out["untapered"][case] = _tip(_run(job, work), job).tolist()
    out["prismatic"] = {}
    for case, (load, _) in STATIC.items():
        from onedfem.pre_processing.mesh_library.blade_beam_job import blade_to_job

        job = _job(f"prismatic_{case}", work, TW, n, prismatic=True, loads={G.LENGTH: blade_to_job(load)})
        out["prismatic"][case] = _tip(_run(job, work), job).tolist()
    return out


def run_modal(work, n) -> dict:
    out = {}
    for name, kw in MODAL_MODELS.items():
        kw = dict(kw)
        model = kw.pop("model", name)
        job = _job(f"{name}_modal", work, model, n, settings=EIGEN, **kw)
        freq, tips = _modes(_run(job, work), job)
        out[name] = {"frequencies_hz": freq, "tips": tips, "kinds": [classify(t, G.TIP["chord"]) for t in tips]}
    return out


def run_nonlinear(work, n, P: float) -> dict:
    from onedfem.pre_processing.mesh_library.blade_beam_job import blade_to_job

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
    """Loft, static differences, frequencies and the large-deflection path. Load steps,
    shell increments and frequencies are discrete values and are drawn as markers."""
    from onedfem.post_processing.graphical_visualisers import figure_style as fs

    fs.use()
    import matplotlib.pyplot as plt

    # 1. the loft: planform with the 20 element cuts, and cuts overlaid
    geo = G.geometry()
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(12, 3.6), gridspec_kw={"width_ratios": [1.9, 1]})
    zs = np.linspace(0, G.LENGTH, 81)
    ext = []
    for z in zs:
        pts = np.vstack([np.asarray(w.xy) for w in geo.section_at(float(z)).walls])
        ext.append((pts[:, 0].min(), pts[:, 0].max()))
    ext = np.array(ext)
    ax.fill_between(zs, ext[:, 0], ext[:, 1], color=fs.GRID, alpha=0.6, lw=0)
    ax.plot(zs, ext[:, 0], color=fs.INK, lw=0.9, label="leading edge")
    ax.plot(zs, ext[:, 1], color=fs.MUTED, lw=0.9, label="trailing edge")
    n_el = int(result.get("n_elements", 20))
    for k in range(n_el):
        z = (k + 0.5) * G.LENGTH / n_el
        for w in geo.section_at(float(z)).walls:
            if w.kind != "skin":
                xy = np.asarray(w.xy)
                ax.plot([z, z], [xy[:, 0].min(), xy[:, 0].max()], color=COLOURS[TW], lw=1.2,
                        solid_capstyle="butt")
    ax.plot([], [], color=COLOURS[TW], lw=1.2, label="caps and webs at each element's mid-point cut")
    ax.set(xlabel="span from the segment root (m)", ylabel="chordwise y (m)")
    fs.panel_title(ax, "IEA 15 MW outboard segment, r = 64.5 to 84.7 m",
                   f"lofted once by implicit-geometry-kernel; {n_el} beam elements, one section cut each")
    ax.legend(loc="upper right", ncol=3, fontsize=7)
    for z, c in ((0.0, fs.colour("", 0)), (0.5 * G.LENGTH, fs.colour("", 2)), (G.LENGTH, fs.colour("", 1))):
        cut = geo.section_at(float(z))
        for i, w in enumerate(cut.walls):
            xy = np.asarray(w.xy)
            ax2.plot(xy[:, 0], xy[:, 1], color=c, lw=0.9 if w.kind == "skin" else 1.6,
                     label=f"{z:.1f} m" if i == 0 else None)
    ax2.set_aspect("equal")
    ax2.set(xlabel="y (m)", ylabel="z (m)")
    fs.panel_title(ax2, "Root, mid and tip cuts", "wall mid-lines; caps and webs thicker")
    ax2.legend(loc="lower center", ncol=3, fontsize=7)
    fig.tight_layout()
    fs.save(fig, os.path.join(out_dir, "loft.png"))

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
    ax.grid(axis="x", visible=False)
    fs.panel_title(ax, "Static tip response, beam against shell", "tip loads 10 kN flap, 10 kN edge, 10 kN·m torque; "
                                                                "grey band: the 3 % target")
    fig.tight_layout()
    fs.save(fig, os.path.join(out_dir, "static.png"))

    # 3. frequencies
    fig, (ax, axd) = plt.subplots(1, 2, figsize=(13, 4.0), gridspec_kw={"width_ratios": [1.6, 1]})
    sf = ref["modal"]["frequencies_hz"]
    kinds = cmp_["modal"]["shell_kinds"]
    for i, (f, k) in enumerate(zip(sf, kinds)):
        ax.plot([0.7, 1.25], [f, f], color="#0b0b0b", lw=1.2)
        ax.text(1.3, f, f"{k}", fontsize=6, va="center", color="#52514e")
    columns = list(MODAL_MODELS)
    for j, model in enumerate(columns):
        m = result["modal"][model]
        for f, k in zip(m["frequencies_hz"], m["kinds"]):
            ax.plot([1.9 + j, 2.45 + j], [f, f], color=COLOURS[model], lw=1.2)
            ax.text(2.5 + j, f, k, fontsize=6, va="center", color="#52514e")
    ax.set_xticks([0.975] + [2.175 + j for j in range(len(columns))],
                  ["CLM shell\n(CalculiX)", "Timoshenko", "Benscoter", "Timoshenko 2D",
                   f"GBT 2D, {GBT_MODES} modes", "GBT 2D, cores stiff\nin transverse shear"])
    ax.tick_params(axis="x", labelsize=6)
    ax.set_ylabel("frequency, Hz", fontsize=8)
    ax.set_xlim(0.4, 2.0 + len(columns) + 0.2)
    ax.grid(axis="x", visible=False)
    fs.panel_title(ax, "Cantilever frequencies", "each line one computed mode, labelled by its kind")
    labels = {"timoshenko_thin_walled": "Timoshenko", "benscoter_thin_walled": "Benscoter",
              "timoshenko_2d": "Timoshenko 2D", "gbt_2d": f"GBT 2D, {GBT_MODES} modes",
              "gbt_2d_stiff_core": "GBT 2D, stiff cores"}
    markers = "osD^v"
    nm = len(columns)
    for j, model in enumerate(columns):
        rows_ = cmp_["modal"][model]
        xs = np.array([r["shell_mode"] for r in rows_], dtype=float) + (j - (nm - 1) / 2) * 0.12
        fs.solves(axd, xs, [100 * r["difference"] for r in rows_], color=COLOURS[model], marker=markers[j % 5],
                  label=labels.get(model, model))
    axd.axhspan(-3, 3, color=fs.GRID, alpha=0.5, lw=0, zorder=0)
    axd.axhline(0, color=fs.MUTED, lw=0.6)
    names_ = [f"{i + 1} {k}" for i, k in enumerate(kinds[:4])]
    axd.set_xticks(range(1, 5), names_)
    axd.grid(axis="x", visible=False)
    axd.set(ylabel="beam against shell (%)")
    fs.panel_title(axd, "First four modes, difference from the shell", "plane-section beams put torsion 7-11 % high")
    axd.legend(loc="upper left", fontsize=6.5)
    fig.tight_layout()
    fs.save(fig, os.path.join(out_dir, "modes.png"))

    # 4. nonlinear tip path
    fig, axes = plt.subplots(1, 2, figsize=(9, 3.4))
    sh = np.array([r[:4] for r in ref["nonlinear"]["tip_series"]], dtype=float)
    P = ref["nonlinear"]["P"]
    lin = result["static"][TW]["flap"][1] / 1e4
    lin_s = ref["static"]["tip"]["flap"][1] / 1e4
    nl = np.array(result["nonlinear"]["rows"], dtype=float)
    ax = axes[0]
    fs.reference(ax, sh[:, 0] * P / 1e6, sh[:, 2], label="CLM shell, CalculiX NLGEOM increments", ms=4)
    fs.solves(ax, nl[:, 0] * P / 1e6, nl[:, 2], color=COLOURS[TW], label="corotational Timoshenko load steps")
    f = np.linspace(0, 0.5, 2)
    ax.plot(f * P / 1e6, lin * f * P, "--", color=COLOURS[TW], lw=0.8, label="linear Timoshenko (exactly linear in load)")
    ax.plot(f * P / 1e6, lin_s * f * P, "--", color=COLOURS["shell"], lw=0.8, label="linear shell (exactly linear in load)")
    ax.set(xlabel="flapwise tip force (MN)", ylabel="flap tip deflection (m)")
    fs.panel_title(ax, "Tip deflection", "within 0.6 % of the shell along the path")
    ax.legend(loc="upper left")
    ax = axes[1]
    fs.reference(ax, sh[:, 2], -sh[:, 3], label="shell increments", ms=4)
    fs.solves(ax, nl[:, 2], -nl[:, 3], color=COLOURS[TW], label="corotational beam load steps")
    ax.set(xlabel="flap tip deflection (m)", ylabel="tip shortening (m)")
    fs.panel_title(ax, "Tip shortening")
    ax.legend(loc="upper left")
    fig.suptitle("Large flapwise deflection, up to 7 % of the span", x=0.01, ha="left")
    fig.tight_layout()
    fs.save(fig, os.path.join(out_dir, "nonlinear.png"))


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


def replot(out_dir: str = HERE) -> None:
    """Redraw the figures from the committed ``summary.json`` and shell reference, without solving."""
    with open(os.path.join(out_dir, "summary.json"), encoding="utf-8") as f:
        result = json.load(f)
    figures(result, reference(), out_dir)


if __name__ == "__main__":
    if "--replot" in sys.argv:
        replot()
        sys.exit(0)
    r = main()
    print(json.dumps(r.get("comparison"), indent=1))

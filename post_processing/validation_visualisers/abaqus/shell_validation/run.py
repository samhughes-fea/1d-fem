# post_processing/validation_visualisers/abaqus/shell_validation/run.py
"""Generate, run and compare the shell validation cases.

Run from the repository root::

    # Abaqus: writes the inputs, runs each job, extracts the ODB and compares
    python -m post_processing.validation_visualisers.abaqus.shell_validation.run --solver abaqus

    # only write the Abaqus inputs (to run them elsewhere), or only compare existing results
    python -m post_processing.validation_visualisers.abaqus.shell_validation.run --solver abaqus --generate-only
    python -m post_processing.validation_visualisers.abaqus.shell_validation.run --solver abaqus --compare-only

    # the same models in CalculiX (quadratic S8R shells), as a cross-check
    python -m post_processing.validation_visualisers.abaqus.shell_validation.run --solver calculix --fine 0.5

The job files and raw results go to
``post_processing/validation_visualisers/abaqus_results/shell_validation/<solver>/<case>/``
(not tracked). The comparison writes ``summary_<solver>.json``, ``report_<solver>.md`` and
one plot per nonlinear case next to this file.

Abaqus is found as in ``abaqus/config.py`` (``ABAQUS_CAE_ROOT`` or ``abaqus`` on the
PATH). The 1d-fem results are cached in ``beam_results.json``; ``--rerun-beam`` runs the
beam jobs again.
"""

from __future__ import annotations

import argparse
import dataclasses
import json
import os
import subprocess
import sys
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from post_processing.validation_visualisers.abaqus.shell_validation.cases import CASES, build_case  # noqa: E402
from post_processing.validation_visualisers.abaqus.shell_validation.results import (  # noqa: E402
    read_abaqus_csv,
    read_calculix_dat,
    station_motions,
)
from post_processing.validation_visualisers.abaqus.shell_validation.shell_model import write_inp  # noqa: E402

RESULTS = ROOT / "post_processing" / "validation_visualisers" / "abaqus_results" / "shell_validation"
BEAM_CACHE = HERE / "beam_results.json"
QUANTITIES = ("u", "v", "w", "phi")

# what each nonlinear case compares: (station x/L, quantity)
COMPARED = {
    "ibeam_ltb_path": [("0.5", "v"), ("0.5", "phi")],
    "ibeam_large_twist": [("1", "phi")],
    "composite_box": [("1", "u"), ("1", "w"), ("1", "phi")],
}
GLOBAL_RESIDUAL = 0.3      # a shell mode is global when the rigid section fit leaves less than this
TWIST_SHARE = 0.2          # ... and torsional when the twist carries more than this share of it


# ------------------------------------------------------------------ solvers
def _abaqus_cmd() -> str:
    from post_processing.validation_visualisers.abaqus.config import ABAQUS_CAE_CMD

    return ABAQUS_CAE_CMD


def run_abaqus(inp: Path, cpus: int) -> Path:
    job = inp.stem
    cmd = [_abaqus_cmd(), f"job={job}", f"input={inp.name}", f"cpus={cpus}", "interactive", "ask_delete=OFF"]
    subprocess.run(cmd, cwd=inp.parent, check=True, shell=(os.name == "nt"))
    odb = inp.with_suffix(".odb")
    csv_path = inp.with_suffix(".csv")
    subprocess.run([_abaqus_cmd(), "python", str(HERE / "extract_odb.py"), odb.name, csv_path.name],
                   cwd=inp.parent, check=True, shell=(os.name == "nt"))
    return csv_path


def run_calculix(inp: Path, ccx: str = "ccx") -> Path:
    with open(inp.with_suffix(".out"), "w") as log:
        subprocess.run([ccx, inp.stem], cwd=inp.parent, check=True, stdout=log, stderr=subprocess.STDOUT)
    return inp.with_suffix(".dat")


# ------------------------------------------------------------------ comparison
def _twist_share(m, radius: float) -> float:
    t = abs(m.phi) * radius
    return t / max(t, abs(m.v), abs(m.w), 1e-300)


def compare_buckling(case, mesh, frames, beam: Dict) -> Dict:
    """The lowest flexural and the lowest torsional global modes of each model, where the
    beam's two lowest modes include one of that kind."""
    from post_processing.validation_visualisers.abaqus.shell_validation.onedfem import _radius

    r = _radius(case.model.section)
    interior = [s for s in case.model.stations if 0.2 <= s <= 0.8]
    motions = station_motions(mesh, frames, interior, linear=True)
    shell = []
    for f, value in enumerate(frames.values):
        res = max(motions[k][f].residual for k in motions)
        mid = motions[min(motions, key=lambda k: abs(float(k) - 0.5))][f]
        shell.append({"factor": value, "residual": res, "twist_share": _twist_share(mid, r),
                      "global": res < GLOBAL_RESIDUAL})
    rows = []
    for kind in ("flexural", "torsional"):
        def is_kind(share):
            return (share > TWIST_SHARE) == (kind == "torsional")

        b = [lam for lam, sh in zip(beam["load_factors"][:2], beam["twist_share"][:2]) if is_kind(sh)]
        s = [m["factor"] for m in shell if m["global"] and is_kind(m["twist_share"])]
        if b:
            rows.append({"mode": kind, "beam": b[0], "shell": s[0] if s else None,
                         "error": (b[0] / s[0] - 1.0) if s else None})
    errors = [abs(r_["error"]) for r_ in rows if r_["error"] is not None]
    return {"rows": rows, "shell_modes": shell, "max_error": max(errors) if errors else None,
            "passed": bool(errors) and len(errors) == len(rows) and max(errors) <= case.target}


def compare_path(case, mesh, frames, beam: Dict, linear: bool = False) -> Dict:
    levels = beam["levels"]
    motions = station_motions(mesh, frames, case.model.stations, linear=linear)
    rows, errors = [], []
    for station, q in COMPARED[case.name]:
        key = min(motions, key=lambda k: abs(float(k) - float(station)))
        bkey = min(beam["stations"], key=lambda k: abs(float(k) - float(station)))
        qi = QUANTITIES.index(q)
        shell_vals, beam_vals = [], []
        for lv, bvals in zip(levels, beam["stations"][bkey]):
            f = int(np.argmin(np.abs(np.asarray(frames.values) - lv)))
            if abs(frames.values[f] - lv) > 1e-6:
                shell_vals.append(None)
            else:
                shell_vals.append(motions[key][f].as_list()[qi])
            beam_vals.append(bvals[qi])
        have = [(s, b) for s, b in zip(shell_vals, beam_vals) if s is not None]
        err = None
        if have:
            # a floor keeps a quantity that stays near zero (u of a linear run) from dividing by noise
            floor = 1e-3 if q == "phi" else 1e-4 * case.model.length
            scale = max(max(abs(s) for s, _ in have), floor)
            err = max(abs(s - b) for s, b in have) / scale
            errors.append(err)
        rows.append({"station": float(station), "quantity": q, "levels": list(levels), "shell": shell_vals,
                     "beam": beam_vals, "error": err})
    complete = all(r_["error"] is not None and None not in r_["shell"] for r_ in rows)
    return {"rows": rows, "max_error": max(errors) if errors else None,
            "passed": complete and max(errors) <= case.target}


def plot_path(case, cmp: Dict, path: Path, solver: str) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    rows = cmp["rows"]
    fig, axes = plt.subplots(1, len(rows), figsize=(4.2 * len(rows), 3.6), squeeze=False)
    labels = {"u": "axial displacement u (m)", "v": "lateral displacement v (m)",
              "w": "vertical displacement w (m)", "phi": "twist (rad)"}
    for ax, r in zip(axes[0], rows):
        lv = np.asarray(r["levels"])
        ax.plot(np.asarray(r["beam"]), lv, "-", color="#1f5fa8", label="1d-fem beam")
        s = [(v, l_) for v, l_ in zip(r["shell"], lv) if v is not None]
        if s:
            ax.plot([v for v, _ in s], [l_ for _, l_ in s], "o", mfc="none", color="#c0392b",
                    label=f"{solver} shell")
        ax.set_xlabel(labels[r["quantity"]])
        ax.set_ylabel("load factor")
        ax.set_title(f"x/L = {r['station']:g}", fontsize=10)
        ax.grid(alpha=0.3)
    axes[0][0].legend(fontsize=8)
    fig.suptitle(case.name.replace("_", " "), fontsize=11)
    fig.tight_layout()
    fig.savefig(path, dpi=120)
    plt.close(fig)


def write_report(summary: Dict, solver: str, path: Path) -> None:
    lines = [f"# Shell validation against {solver}", "",
             f"Shell mesh density factor {summary['fine']:g}. Target: each case within its tolerance "
             "(3 %) of the shell model.", "",
             "| case | compared | 1d-fem | shell | difference | passed |",
             "|---|---|---|---|---|---|"]
    for name, c in summary["cases"].items():
        if "error" in c and isinstance(c["error"], str):
            lines.append(f"| {name} | {c['error']} | | | | no |")
            continue
        for r in c["rows"]:
            if "mode" in r:
                shell = f"{r['shell']:.5g}" if r["shell"] is not None else "not found"
                diff = f"{100 * r['error']:+.2f} %" if r["error"] is not None else ""
                lines.append(f"| {name} | {r['mode']} buckling factor | {r['beam']:.5g} | {shell} | {diff} | "
                             f"{'yes' if c['passed'] else 'no'} |")
            else:
                last_b = r["beam"][-1]
                last_s = r["shell"][-1]
                diff = f"{100 * r['error']:.2f} % of max" if r["error"] is not None else ""
                lines.append(f"| {name} | {r['quantity']} at x/L = {r['station']:g}, final level | {last_b:.5g} | "
                             f"{'' if last_s is None else f'{last_s:.5g}'} | {diff} | "
                             f"{'yes' if c['passed'] else 'no'} |")
    lines += ["", "Path differences are the largest difference over the load levels, relative to the largest "
              "shell value of that quantity (at least 1e-3 rad for a twist and 1e-4 L for a displacement)."]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


# ------------------------------------------------------------------ main
def beam_results(names: List[str], rerun: bool, work: Path, linear: bool = False) -> Dict:
    """1d-fem results per case, cached in ``beam_results.json`` (linear ones under ``<case>/linear``)."""
    from post_processing.validation_visualisers.abaqus.shell_validation.onedfem import run_beam

    cache = json.loads(BEAM_CACHE.read_text()) if BEAM_CACHE.exists() else {}
    out = {}
    for name in names:
        lin = linear and build_case(name).beam["kind"] == "nonlinear"
        key = f"{name}/linear" if lin else name
        if rerun or key not in cache:
            cache[key] = run_beam(build_case(name), str(work / "beam"), linear=lin)
            BEAM_CACHE.write_text(json.dumps(cache, indent=1, sort_keys=True) + "\n")
        out[name] = cache[key]
    return out


def main(argv: Optional[List[str]] = None) -> Dict:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--solver", choices=("abaqus", "calculix"), default="abaqus")
    ap.add_argument("--cases", nargs="*", default=list(CASES))
    ap.add_argument("--fine", type=float, default=1.0, help="shell mesh density factor (1 = validation mesh)")
    ap.add_argument("--element", default=None, help="shell element (default S4R for Abaqus, S8R for CalculiX)")
    ap.add_argument("--cpus", type=int, default=4)
    ap.add_argument("--generate-only", action="store_true")
    ap.add_argument("--compare-only", action="store_true")
    ap.add_argument("--rerun-beam", action="store_true")
    ap.add_argument("--linear", action="store_true",
                    help="run the nonlinear cases without NLGEOM and compare with the linear beams at load 1 "
                         "(a check of the decks; CalculiX's expanded shells do not converge with NLGEOM)")
    args = ap.parse_args(argv)
    element = args.element or ("S4R" if args.solver == "abaqus" else "S8R")
    base = RESULTS / args.solver
    summary = {"solver": args.solver, "fine": args.fine, "element": element, "linear": args.linear, "cases": {}}
    tag = f"{args.solver}_linear" if args.linear else args.solver
    meshes = {}
    for name in args.cases:
        case = build_case(name, args.fine)
        case.model.element_type = element
        if args.linear:
            case.model.step = dataclasses.replace(case.model.step, nlgeom=False)
        d = base / name
        d.mkdir(parents=True, exist_ok=True)
        inp = d / f"{name}.inp"
        if not args.compare_only:
            meshes[name] = write_inp(case.model, str(inp), dialect=args.solver)
            print(f"wrote {inp}")
        else:
            meshes[name] = case.model.mesh()
    if args.generate_only:
        return summary
    beams = beam_results(args.cases, args.rerun_beam, RESULTS, args.linear)
    for name in args.cases:
        case = build_case(name, args.fine)
        case.model.element_type = element
        d = base / name
        inp = d / f"{name}.inp"
        try:
            if args.solver == "abaqus":
                out = inp.with_suffix(".csv") if args.compare_only else run_abaqus(inp, args.cpus)
                frames = read_abaqus_csv(str(out))
            else:
                out = inp.with_suffix(".dat") if args.compare_only else run_calculix(inp)
                frames = read_calculix_dat(str(out), case.model.step.kind)
        except (OSError, subprocess.CalledProcessError, ValueError) as exc:
            summary["cases"][name] = {"error": f"no shell results ({exc.__class__.__name__})", "passed": False}
            print(f"{name}: no shell results: {exc}")
            continue
        if case.model.step.kind == "buckle":
            cmp = compare_buckling(case, meshes[name], frames, beams[name])
        else:
            cmp = compare_path(case, meshes[name], frames, beams[name], linear=args.linear)
            if not args.linear:
                plot_path(case, cmp, HERE / f"{name}_{args.solver}.png", args.solver)
        summary["cases"][name] = cmp
        err = cmp["max_error"]
        print(f"{name}: {'passed' if cmp['passed'] else 'FAILED'}"
              + (f", largest difference {100 * err:.2f} %" if err is not None else ""))
    (HERE / f"summary_{tag}.json").write_text(json.dumps(summary, indent=1) + "\n")
    write_report(summary, tag, HERE / f"report_{tag}.md")
    return summary


if __name__ == "__main__":
    main()

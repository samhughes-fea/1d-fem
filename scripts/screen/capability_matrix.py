"""Probe 1 of the 1d-fem screen: run every registered element through every analysis.

Each cell is one small cantilever job: an element type from ``ElementFactory``, an
analysis type, a section model and a geometry. The job runs in its own process with a
time limit and the cell records what happened:

- ``ok``: the job finished and its key result is finite;
- ``refused``: the job stopped with a clear ``NotImplementedError`` or ``ValueError``;
- ``crash``: any other exception, a non-finite result or the time limit;
- ``differs``: it ran, but its key result is more than ``TOL`` away from the reference
  element (``LinearTimoshenkoBeamElement3D``, or the linear counterpart for
  ``static_nonlinear``) on the same section and geometry.

The cantilever is slender (L / h = 20) and lightly loaded, so all beam theories should
agree to well within ``TOL``. Bars and trusses carry axial load only, so only their axial
response is compared.

Usage::

    python scripts/screen/capability_matrix.py --out screen_out --workers 6
    python scripts/screen/capability_matrix.py --only LinearTimoshenkoBeamElement3D

Writes ``capability_matrix.csv`` (one row per cell) and ``capability_matrix.json``.
"""

from __future__ import annotations

import argparse
import contextlib
import csv
import io
import itertools
import json
import logging
import os
import shutil
import subprocess
import sys
import tempfile
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
ANALYSES = ("static", "static_nonlinear", "eigen", "buckling", "harmonic", "transient", "layup_sensitivity")
SECTIONS = ("closed_form", "thin_walled", "solid_2d")
GEOMETRIES = ("straight", "framed", "curved")
AXIAL_ONLY = ("LinearBarElement3D", "LinearTrussElement3D")
TOL = 0.03
TIMEOUT = 180

L, N_EL = 4.0, 8
B, H, T = 0.3, 0.2, 0.005            # box section, metres
E, NU, RHO = 70e9, 0.3, 2700.0
E2 = 50e9                            # laminate transverse modulus, so ply angles matter
G = E / (2 * (1 + NU))
DOFS = ("UX", "UY", "UZ", "RX", "RY", "RZ")


def element_types():
    sys.path.insert(0, str(ROOT / "src"))
    from onedfem.pre_processing.element_library.element_factory import ElementFactory

    return sorted(ElementFactory.ELEMENT_CLASS_MAP)


# ----------------------------------------------------------------------------- job writing
def _box_closed_form():
    bo, ho, bi, hi = B + T, H + T, B - T, H - T
    A = bo * ho - bi * hi
    Iy = (bo * ho**3 - bi * hi**3) / 12.0
    Iz = (ho * bo**3 - hi * bi**3) / 12.0
    J = 4 * (B * H) ** 2 * T / (2 * (B + H))
    return A, Iy, Iz, J


def _axis(geometry):
    if geometry == "framed":
        d = np.array([1.0, 1.0, 1.0]) / np.sqrt(3.0)
    else:
        d = np.array([1.0, 0.0, 0.0])
    return np.outer(np.linspace(0.0, L, N_EL + 1), d), d


def _settings(analysis):
    solver = "\n[Solver]\ntype = direct\ntolerance = 1e-12\nmax_iterations = 500\n"
    if analysis in ("static", "layup_sensitivity"):
        return "[Simulation]\n[Type]\nstatic\n" + solver
    if analysis == "static_nonlinear":
        return ("[Simulation]\n[Type]\nstatic_nonlinear\n\n[Newton]\ntolerance = 1e-9\nrelative_tolerance = 1e-10\n"
                "max_iterations = 40\n" + solver + "\n[Nonlinear]\nload_steps = 2\n")
    if analysis == "eigen":
        return "[Simulation]\n[Type]\neigen\n\n[Eigen]\nenabled = true\nnum_modes = 4\n"
    if analysis == "buckling":
        return ("[Simulation]\n[Type]\nbuckling\n\n[Buckling]\nenabled = true\nnum_modes = 2\n"
                "buckling_prestress = linear_static\nbuckling_load_factor = 1.0\n")
    if analysis == "harmonic":
        return ("[Simulation]\n[Type]\nharmonic\n\n[Harmonic]\nenabled = true\nfrequency_min_hz = 1.0\n"
                "frequency_max_hz = 2.0\nnum_frequency_points = 2\nmodal_damping_ratio = 0.02\n")
    if analysis == "transient":
        return ("[Simulation]\n[Type]\ntransient\n\n[Transient]\nenabled = true\ntime_step = 0.002\nend_time = 0.01\n"
                "load_scale = 1.0\n\n[PostProcessing]\nrun_secondary_tertiary_dynamic = false\n")
    raise ValueError(analysis)


def _tip_load(etype, analysis, d):
    if etype in AXIAL_ONLY:
        return np.concatenate([1e3 * d, np.zeros(3)])
    if analysis == "buckling":
        return np.concatenate([-1e3 * d, np.zeros(3)])
    # flap force normal to the axis and a torque about it
    n = np.array([0.0, 0.0, 1.0]) if abs(d[2]) < 0.9 else np.array([1.0, 0.0, 0.0])
    n = n - d * (n @ d)
    n /= np.linalg.norm(n)
    return np.concatenate([1e2 * d + 1e3 * n, 2e2 * d])


def _thin_walled_files(job):
    pts = {1: (-B / 2, -H / 2), 2: (B / 2, -H / 2), 3: (B / 2, H / 2), 4: (-B / 2, H / 2)}
    (job / "laminate.txt").write_text(
        "[Ply_Material]\n[material_id] [E1] [E2] [G12] [nu12] [G13] [G23] [rho]\n"
        f"1 {E} {E2} {G} {NU} {G} {G} {RHO}\n"
        "[Laminate]\n[laminate_id] [ply] [material_id] [angle_deg] [thickness]\n"
        f"1 1 1 30 {T / 2}\n1 2 1 -30 {T / 2}\n")
    (job / "section_geometry.txt").write_text(
        "[Section_Point]\n[section_id] [point_id] [y] [z]\n"
        + "".join(f"1 {p} {y} {z}\n" for p, (y, z) in pts.items())
        + "[Section_Wall]\n[section_id] [wall_id] [point_start] [point_end] [laminate_id]\n"
        + "1 1 1 2 1\n1 2 2 3 1\n1 3 3 4 1\n1 4 4 1 1\n"
        + "[Section_Mesh]\n[section_id] [n_sub]\n1 2\n")


def _solid_files(job):
    sys.path.insert(0, str(ROOT / "src"))
    from onedfem.pre_processing.section_integration.solid.from_geometry import write_section_mesh_2d
    from onedfem.pre_processing.section_integration.solid.mesh import SectionMesh
    from onedfem.pre_processing.section_integration.thin_walled import ThinWalledSection, Wall, WallLaminate
    from onedfem.pre_processing.section_integration.thin_walled.laminate import Ply, PlyMaterial

    mat = PlyMaterial(E, E2, G, NU, G, G, RHO)
    lam = WallLaminate([Ply(mat, 30.0, T / 2), Ply(mat, -30.0, T / 2)])
    pts = {1: (-B / 2, -H / 2), 2: (B / 2, -H / 2), 3: (B / 2, H / 2), 4: (-B / 2, H / 2)}
    sec = ThinWalledSection(pts, [Wall(1, 1, 2, 1), Wall(2, 2, 3, 1), Wall(3, 3, 4, 1), Wall(4, 4, 1, 1)])
    mesh = SectionMesh.from_thin_walled(sec, {1: lam}, elements_along=4)
    write_section_mesh_2d(str(job / "section_mesh_2d.txt"), {1: mesh})


def write_job(job: Path, etype: str, analysis: str, section: str, geometry: str) -> None:
    job.mkdir(parents=True)
    xyz, d = _axis(geometry)
    A, Iy, Iz, J = _box_closed_form()
    warping = int("Benscoter" in etype or "GBT" in etype)
    (job / "grid.txt").write_text("[Grid]\n[node_id] [x] [y] [z]\n"
                                  + "".join(f"{i} {x!r} {y!r} {z!r}\n" for i, (x, y, z) in enumerate(xyz.tolist())))
    (job / "element.txt").write_text(
        "[Element]\n[element_id] [node1] [node2] [element_type] [axial_order] [bending_y_order] [bending_z_order] "
        "[shear_y_order] [shear_z_order] [torsion_order] [load_order] [warping]\n"
        + "".join(f"{e} {e} {e + 1} {etype} 3 3 3 2 2 3 3 {warping}\n" for e in range(N_EL)))
    (job / "material.txt").write_text("[Material]\n[element_id] [E] [G] [nu] [rho]\n"
                                      + "".join(f"{e} {E} {G} {NU} {RHO}\n" for e in range(N_EL)))
    (job / "section.txt").write_text("[Section]\n[element_id] [A] [I_x] [I_y] [I_z] [J_t] [kappa] [alpha]\n"
                                     + "".join(f"{e} {A} {Iy + Iz} {Iy} {Iz} {J} 0.5 0\n" for e in range(N_EL)))
    k = (0.0, 0.02, 0.0) if geometry == "curved" else (0.0, 0.0, 0.0)
    (job / "precurvature.txt").write_text("[Precurvature]\n[element_id] [k_x0] [k_y0] [k_z0]\n"
                                          + "".join(f"{e} {k[0]} {k[1]} {k[2]}\n" for e in range(N_EL)))
    (job / "prescribed_displacement.txt").write_text(
        "[Prescribed Displacement]\n[id] [node_id] [dof] [value] [type] [comment]\n"
        + "".join(f"{i} 0 {dof} 0.0 displacement\n" for i, dof in enumerate(DOFS + ("Q*",))))
    P = _tip_load(etype, analysis, d)
    tip = xyz[-1]
    (job / "point_load.txt").write_text("[Point load]\n[x] [y] [z] [F_x] [F_y] [F_z] [M_x] [M_y] [M_z]\n"
                                        + " ".join(repr(float(v)) for v in (*tip, *P)) + "\n")
    (job / "simulation_settings.txt").write_text(_settings(analysis))
    if section in ("thin_walled", "solid_2d"):
        if section == "thin_walled":
            _thin_walled_files(job)
        else:
            _solid_files(job)
        modes = " [modes] [shear_modes]" if "GBT" in etype else ""
        row = " 1 0" if "GBT" in etype else ""
        (job / "element_section.txt").write_text(f"[Element_Section]\n[element_id] [section_id]{modes}\n"
                                                 + "".join(f"{e} 1{row}\n" for e in range(N_EL)))
    if analysis == "layup_sensitivity":
        params = "1 1 angle\n" + ("" if section == "solid_2d" else "1 1 thickness\n")
        (job / "layup_sensitivity.txt").write_text(
            "[Layup_Parameter]\n[laminate_id] [ply] [kind]\n" + params +
            f"[Response]\n[type] [node_id] [dof] [mode]\ndisplacement {N_EL} UZ -\ncompliance - - -\n")


# ----------------------------------------------------------------------------- one cell
def _read_csv(path):
    with open(path, newline="", encoding="utf-8") as f:
        return list(csv.reader(f))


def _find(res: Path, names):
    for name in names:
        hits = sorted(res.rglob(name))
        if hits:
            return hits[0]
    return None


def key_result(res: Path, analysis: str, etype: str):
    """The cell's key number: tip displacement norm, first frequency, load factor, ..."""
    modal = res / "primary_results" / "modal_results"
    if analysis == "eigen":
        return float(open(modal / "job_frequencies.txt").read().split()[0])
    if analysis == "buckling":
        return float(open(modal / "job_buckling_load_factors.txt").read().split()[0])
    if analysis == "harmonic":
        vals = np.loadtxt(res / "primary_results" / "harmonic_results" / "job_displacement_abs.txt")
        return float(np.max(np.abs(vals)))
    if analysis == "transient":
        f = _find(res / "primary_results", ["*displacement*.txt", "*displacement*.csv", "*U_*.csv", "*U*.txt"])
        vals = np.genfromtxt(f, delimiter="," if f.suffix == ".csv" else None)
        vals = vals[np.isfinite(vals)]
        return float(np.max(np.abs(vals)))
    if analysis == "layup_sensitivity":
        f = _find(res, ["layup_sensitivity.csv"])
        with open(f, newline="", encoding="utf-8") as fh:
            rows = list(csv.DictReader(fh))
        return float(max(abs(float(r["derivative"])) for r in rows))
    f = res / "primary_results" / "global" / "U_global.csv"
    U = np.array([float(r[1]) for r in _read_csv(f)[1:]])
    tip = U.reshape(N_EL + 1, -1)[-1, :6]
    return float(abs(tip[0])) if etype in AXIAL_ONLY else float(np.linalg.norm(tip[:3]))


def _is_float(x):
    try:
        float(x)
        return True
    except ValueError:
        return False


def run_one(etype, analysis, section, geometry) -> dict:
    sys.path.insert(0, str(ROOT / "src"))
    logging.disable(logging.WARNING)
    from onedfem.workflow_orchestrator.run_job import process_job

    work = Path(tempfile.mkdtemp(prefix="screen_"))
    job, res = work / "job", work / "res"
    try:
        write_job(job, etype, analysis, section, geometry)
        for sub in ("element_stiffness_matrices", "element_force_vectors", "primary_results", "secondary_results",
                    "logs", "maps", "diagnostics"):
            (res / sub).mkdir(parents=True)
        try:
            with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
                process_job(str(job), str(res), {}, {}, force_serial=True, max_processes_per_job=1)
        except (NotImplementedError, ValueError) as exc:
            return {"status": "refused", "message": f"{type(exc).__name__}: {exc}"[:300]}
        except Exception as exc:  # noqa: BLE001 - every other failure is a finding
            return {"status": "crash", "message": f"{type(exc).__name__}: {exc}"[:300]}
        tb = res / "logs" / "traceback.log"
        if tb.exists() and tb.read_text().strip():
            # process_job logs the failure and returns normally
            last = [ln for ln in tb.read_text().splitlines() if ln.strip()][-1].strip()
            kind = last.split(":")[0]
            status = "refused" if kind in ("NotImplementedError", "ValueError") else "crash"
            return {"status": status, "message": ("(logged, not raised) " + last)[:300]}
        try:
            value = key_result(res, analysis, etype)
        except Exception as exc:  # noqa: BLE001
            return {"status": "crash", "message": f"no readable result: {type(exc).__name__}: {exc}"[:300]}
        if not np.isfinite(value):
            return {"status": "crash", "message": "non-finite result", "value": None}
        return {"status": "ok", "value": value}
    finally:
        shutil.rmtree(work, ignore_errors=True)


# ----------------------------------------------------------------------------- driver
def _cell_subprocess(cell):
    t0 = time.time()
    try:
        p = subprocess.run([sys.executable, __file__, "--cell", *cell], capture_output=True, text=True,
                           timeout=TIMEOUT, cwd=ROOT)
        line = [ln for ln in p.stdout.splitlines() if ln.startswith("{")]
        out = json.loads(line[-1]) if line else {"status": "crash",
                                                  "message": (p.stderr.strip().splitlines() or ["no output"])[-1][:300]}
    except subprocess.TimeoutExpired:
        out = {"status": "crash", "message": f"time limit {TIMEOUT} s"}
    out.update(element=cell[0], analysis=cell[1], section=cell[2], geometry=cell[3], seconds=round(time.time() - t0, 1))
    return out


def reference_type(etype, analysis):
    return "LinearTimoshenkoBeamElement3D"


def compare(rows):
    by = {(r["element"], r["analysis"], r["section"], r["geometry"]): r for r in rows}
    for r in rows:
        r.setdefault("reference_value", None)
        r.setdefault("relative_difference", None)
        if r["status"] != "ok" or r["element"] in AXIAL_ONLY:
            continue
        ref_analysis = "static" if r["analysis"] == "static_nonlinear" else r["analysis"]
        ref = by.get((reference_type(r["element"], r["analysis"]), ref_analysis, r["section"], r["geometry"]))
        if ref is None or ref["status"] != "ok" or not ref.get("value"):
            continue
        d = abs(r["value"] - ref["value"]) / abs(ref["value"])
        r["reference_value"], r["relative_difference"] = ref["value"], d
        if d > TOL:
            r["status"] = "differs"
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cell", nargs=4)
    ap.add_argument("--out", default=str(ROOT / "screen_out"))
    ap.add_argument("--workers", type=int, default=os.cpu_count() or 2)
    ap.add_argument("--only", nargs="*")
    ap.add_argument("--analyses", nargs="*", default=list(ANALYSES))
    a = ap.parse_args()
    if a.cell:
        print(json.dumps(run_one(*a.cell)))
        return
    types = a.only or element_types()
    cells = list(itertools.product(types, a.analyses, SECTIONS, GEOMETRIES))
    if "LinearTimoshenkoBeamElement3D" not in types:
        cells += list(itertools.product(["LinearTimoshenkoBeamElement3D"], a.analyses, SECTIONS, GEOMETRIES))
    with ThreadPoolExecutor(a.workers) as ex:
        rows = []
        for i, r in enumerate(ex.map(_cell_subprocess, cells), 1):
            rows.append(r)
            if i % 25 == 0:
                print(f"{i}/{len(cells)}", flush=True)
    rows = compare(rows)
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    cols = ["element", "analysis", "section", "geometry", "status", "value", "reference_value", "relative_difference",
            "seconds", "message"]
    with open(out / "capability_matrix.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    (out / "capability_matrix.json").write_text(json.dumps(rows, indent=1))
    counts = {}
    for r in rows:
        counts[r["status"]] = counts.get(r["status"], 0) + 1
    print(json.dumps(counts))


if __name__ == "__main__":
    main()

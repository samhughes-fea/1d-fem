# post_processing/validation_visualisers/abaqus/shell_validation/onedfem.py
"""The 1d-fem side of the shell validation: write and run the beam jobs of a case.

The jobs go through the normal job runner (``process_job``) with a thin-walled section,
so the comparison checks the whole path a user takes: section files, section analysis,
element, supports, loads and solver.

Public API
----------
section_analysis, tension_centre, run_beam.
"""

from __future__ import annotations

import csv
import glob
import logging
import os
import shutil
from typing import Dict, Sequence, Tuple

import numpy as np

from post_processing.validation_visualisers.abaqus.shell_validation.shell_model import ShellSection

STANDARD_DOFS = ("UX", "UY", "UZ", "RX", "RY", "RZ")


def section_analysis(sec: ShellSection, n_sub: int = 2):
    """1d-fem's ``ThinWalledSectionAnalysis`` of a shell-model section."""
    from pre_processing.section_integration.thin_walled import (
        ThinWalledSection,
        ThinWalledSectionAnalysis,
        Wall,
        WallLaminate,
    )
    from pre_processing.section_integration.thin_walled.laminate import Ply, PlyMaterial

    mats = {}
    lams = {}
    for lid, lam in sec.laminates.items():
        plies = []
        for p in lam.plies:
            m = p.material
            if m.name not in mats:
                mats[m.name] = PlyMaterial(m.E1, m.E2, m.G12, m.nu12, m.G13, m.G23, m.rho)
            plies.append(Ply(mats[m.name], p.angle_deg, p.thickness))
        lams[lid] = WallLaminate(plies)
    walls = [Wall(wid, a, b, lid) for wid, a, b, lid in sec.walls]
    return ThinWalledSectionAnalysis(ThinWalledSection(dict(sec.points), walls), lams, n_sub=n_sub)


def tension_centre(sec: ShellSection) -> Tuple[float, float]:
    tc = section_analysis(sec).tension_centre
    return float(tc[0]), float(tc[1])


def _write(path: str, text: str) -> None:
    with open(path, "w", encoding="utf-8") as f:
        f.write(text)


def write_beam_job(job_dir: str, sec: ShellSection, beam: Dict, settings: str, scale: float = 1.0) -> None:
    """Write a 1d-fem job for the beam description of a case (loads times ``scale``)."""
    if os.path.isdir(job_dir):
        shutil.rmtree(job_dir)
    os.makedirs(job_dir)
    n, L, et, warping = beam["n_el"], beam["length"], beam["element"], beam["warping"]
    _write(os.path.join(job_dir, "grid.txt"),
           "[Grid]\n[node_id] [x] [y] [z]\n" + "".join(f"{i} {L * i / n:.12g} 0 0\n" for i in range(n + 1)))
    _write(os.path.join(job_dir, "element.txt"),
           "[Element]\n[element_id] [node1] [node2] [element_type] [axial_order] [bending_y_order] "
           "[bending_z_order] [shear_y_order] [shear_z_order] [torsion_order] [load_order] [warping]\n"
           + "".join(f"{e} {e} {e + 1} {et} 2 2 2 2 2 2 2 {warping}\n" for e in range(n)))
    # the classical inputs are required by the parsers but unused with a thin-walled section
    _write(os.path.join(job_dir, "material.txt"),
           "[Material]\n[element_id] [E] [G] [nu] [rho]\n" + "".join(f"{e} 7e10 2.7e10 0.3 0\n" for e in range(n)))
    _write(os.path.join(job_dir, "section.txt"),
           "[Section]\n[element_id] [A] [I_x] [I_y] [I_z] [J_t] [kappa] [alpha]\n"
           + "".join(f"{e} 1e-3 0.0 1e-8 1e-8 1e-9 0.8333 0.0\n" for e in range(n)))
    rows = []
    for x, vals in beam["loads"]:
        rows.append(f"{x:.12g} 0 0 " + " ".join(f"{scale * v:.12g}" for v in vals) + "\n")
    _write(os.path.join(job_dir, "point_load.txt"),
           "[Point load]\n[x] [y] [z] [F_x] [F_y] [F_z] [M_x] [M_y] [M_z]\n" + "".join(rows))
    if beam["supports"] == "fork":
        fixed = [(0, d) for d in ("UX", "UY", "UZ", "RX")] + [(n, d) for d in ("UY", "UZ", "RX")]
    elif beam["supports"] == "clamp":
        fixed = [(0, d) for d in STANDARD_DOFS] + ([(0, "Q1")] if warping else [])
    else:
        raise ValueError(f"unknown supports {beam['supports']!r}")
    _write(os.path.join(job_dir, "prescribed_displacement.txt"),
           "[Prescribed Displacement]\n[id] [node_id] [dof] [value] [type] [comment]\n"
           + "".join(f"{i} {nd} {d} 0.0 displacement\n" for i, (nd, d) in enumerate(fixed)))
    _write(os.path.join(job_dir, "simulation_settings.txt"), settings)

    mats: Dict[str, int] = {}
    lam_rows = []
    mat_rows = []
    for lid, lam in sec.laminates.items():
        for k, p in enumerate(lam.plies):
            m = p.material
            if m.name not in mats:
                mats[m.name] = len(mats) + 1
                mat_rows.append(f"{mats[m.name]} {m.E1:.12g} {m.E2:.12g} {m.G12:.12g} {m.nu12:.12g} "
                                f"{m.G13:.12g} {m.G23:.12g} {m.rho:.12g}\n")
            lam_rows.append(f"{lid} {k + 1} {mats[m.name]} {p.angle_deg:.12g} {p.thickness:.12g}\n")
    _write(os.path.join(job_dir, "laminate.txt"),
           "[Ply_Material]\n[material_id] [E1] [E2] [G12] [nu12] [G13] [G23] [rho]\n" + "".join(mat_rows)
           + "[Laminate]\n[laminate_id] [ply] [material_id] [angle_deg] [thickness]\n" + "".join(lam_rows))
    _write(os.path.join(job_dir, "section_geometry.txt"),
           "[Section_Point]\n[section_id] [point_id] [y] [z]\n"
           + "".join(f"1 {p} {y:.12g} {z:.12g}\n" for p, (y, z) in sec.points.items())
           + "[Section_Wall]\n[section_id] [wall_id] [point_start] [point_end] [laminate_id]\n"
           + "".join(f"1 {wid} {a} {b} {lid}\n" for wid, a, b, lid in sec.walls)
           + "[Section_Mesh]\n[section_id] [n_sub]\n1 2\n")
    _write(os.path.join(job_dir, "element_section.txt"),
           "[Element_Section]\n[element_id] [section_id]\n" + "".join(f"{e} 1\n" for e in range(n)))


def _results_dir(base: str) -> str:
    if os.path.isdir(base):
        shutil.rmtree(base)
    for sub in ("element_stiffness_matrices", "element_force_vectors", "primary_results", "secondary_results",
                "logs", "maps", "diagnostics"):
        os.makedirs(os.path.join(base, sub))
    return base


def _process(job_dir: str, res_dir: str) -> None:
    from workflow_orchestrator.run_job import process_job

    logging.getLogger().setLevel(logging.WARNING)
    process_job(job_dir, _results_dir(res_dir), {}, {}, force_serial=True, max_processes_per_job=1)


def _read_U(res_dir: str, n_el: int) -> np.ndarray:
    with open(os.path.join(res_dir, "primary_results", "global", "U_global.csv"), newline="") as fh:
        U = np.array([float(r[1]) for r in list(csv.reader(fh))[1:]])
    return U.reshape(n_el + 1, -1)


BUCKLING_SETTINGS = ("[Simulation]\n[Type]\nbuckling\n\n[Solver]\ntype = direct\ntolerance = 1e-12\n\n"
                     "[Buckling]\nenabled = true\nnum_modes = {n}\nbuckling_prestress = linear_static\n"
                     "buckling_load_factor = 1.0\n")
LINEAR_SETTINGS = ("[Simulation]\n[Type]\nstatic\n\n[Solver]\ntype = direct\ntolerance = 1e-12\n"
                   "max_iterations = 500\n")
NONLINEAR_SETTINGS = ("[Simulation]\n[Type]\nstatic_nonlinear\n\n[Newton]\ntolerance = {tol:.6g}\n"
                      "relative_tolerance = 1e-8\nrelative_reference = first_residual\nmax_iterations = 40\n"
                      "tolerance_delta_u = 1e-10\n\n[Solver]\ntype = direct\ntolerance = 1e-12\n\n"
                      "[Nonlinear]\nnum_increments = {inc}\nline_search = false\n")


def run_beam(case, work_dir: str, n_modes: int = 4, linear: bool = False) -> Dict:
    """Run the 1d-fem jobs of ``case`` under ``work_dir``; returns the beam results.

    ``linear`` runs a nonlinear case once at load level 1 with the linear element of the
    same theory (the check of a shell deck run without NLGEOM).

    Buckling: ``{"load_factors": [...], "twist_share": [...]}``, the factors sorted and,
    per mode, the share of the twist in the mode (to tell flexural from torsional modes).
    Nonlinear: ``{"levels": [...], "stations": {x/L: [[u, v, w, phi], ...]}}`` with one row
    per load level at the beam nodes nearest the case's stations.
    """
    beam = case.beam
    sec = case.model.section
    n = beam["n_el"]
    if beam["kind"] == "buckling":
        job = os.path.join(work_dir, "jobs", case.name)
        res = os.path.join(work_dir, "results", case.name)
        write_beam_job(job, sec, beam, BUCKLING_SETTINGS.format(n=n_modes))
        _process(job, res)
        f = glob.glob(os.path.join(res, "primary_results", "modal_results", "*_buckling_load_factors.txt"))[0]
        lam = np.atleast_1d(np.loadtxt(f))
        modes = np.loadtxt(f.replace("load_factors", "mode_shapes")).reshape(-1, lam.size)
        stride = modes.shape[0] // (n + 1)
        order = np.argsort(lam)
        r = _radius(sec)
        share = []
        for k in order:
            m = modes[:, k].reshape(n + 1, stride)
            t = np.abs(m[:, 3]).max() * r
            share.append(float(t / max(t, np.abs(m[:, 1:3]).max(), 1e-300)))
        return {"load_factors": lam[order].tolist(), "twist_share": share}

    if beam["kind"] != "nonlinear":
        raise ValueError(f"unknown beam kind {beam['kind']!r}")
    stations = case.model.stations
    nodes = sorted({int(round(s * n)) for s in stations})
    levels = [1.0] if linear else list(beam["levels"])
    if linear:
        beam = dict(beam, element=beam["element"].replace("Corotational", "Linear"))
    out = {"levels": levels, "stations": {f"{k / n:.6g}": [] for k in nodes}}
    for lv in levels:
        name = f"{case.name}_{'linear' if linear else f'{lv:.6f}'}"
        job = os.path.join(work_dir, "jobs", name)
        res = os.path.join(work_dir, "results", name)
        F_ref = max(abs(v) for _, vals in beam["loads"] for v in vals) * lv
        settings = (LINEAR_SETTINGS if linear else
                    NONLINEAR_SETTINGS.format(tol=1e-8 * max(F_ref, 1.0), inc=beam["increments"]))
        write_beam_job(job, sec, beam, settings, scale=lv)
        _process(job, res)
        if not linear:
            with open(os.path.join(res, "logs", "newton_history.csv"), newline="") as fh:
                last = list(csv.DictReader(fh))[-1]
            if float(last["load_factor"]) < 1.0 - 1e-12 or last["residual_ok"] != "1":
                raise RuntimeError(f"{name}: Newton did not reach the full load (last row {last})")
        U = _read_U(res, n)
        for k in nodes:
            out["stations"][f"{k / n:.6g}"].append([float(U[k, 0]), float(U[k, 1]), float(U[k, 2]), float(U[k, 3])])
    return out


def _radius(sec: ShellSection) -> float:
    """Half the largest section dimension: turns a twist into a comparable displacement."""
    yz = np.array(list(sec.points.values()), dtype=float)
    return 0.5 * float(np.max(yz.max(axis=0) - yz.min(axis=0)))


__all__: Sequence[str] = ["section_analysis", "tension_centre", "write_beam_job", "run_beam"]

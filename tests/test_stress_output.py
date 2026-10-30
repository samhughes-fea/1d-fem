"""Stress output: stations, corotational and modal runs, and ply-interface interlaminar values.

References:

- A [0/90/0/90] strip under V_z (nu = 0) has the layered-beam flexure shear
  ``tau_xz(z) = V int_z^top E (z' - z_n) dz' / (EI b)``, continuous across the plies.
- At a small load a corotational job is the linear job, stresses included.
- A mode shape's wall stresses are those of a static state equal to the mode shape.
"""

from __future__ import annotations

import csv
import logging
import os
import sys

import numpy as np
import pytest

from pre_processing.section_integration.solid.material import SolidMaterial
from pre_processing.section_integration.solid.mesh import SectionMesh
from pre_processing.section_integration.solid.section_analysis import SolidSectionAnalysis
from pre_processing.section_integration.thin_walled import ThinWalledSection, Wall, WallLaminate
from pre_processing.section_integration.thin_walled.laminate import Ply, PlyMaterial

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "examples", "thin_walled"))
import common  # noqa: E402

MAT = SolidMaterial(135e9, 10e9, 10e9, 5e9, 5e9, 3.8e9, 0.0, 0.0, 0.0)
PLY = PlyMaterial(E1=135e9, E2=10e9, G12=5e9, nu12=0.0, G13=5e9, G23=3.8e9)
T, B = 0.5e-3, 0.02
ANGLES = (0, 90, 0, 90)


def _strip(rows_per_ply):
    sec = ThinWalledSection({1: (-B / 2, 0.0), 2: (B / 2, 0.0)}, [Wall(1, 1, 2, 1)])
    lam = {1: WallLaminate([Ply(PLY, a, T) for a in ANGLES])}
    return SectionMesh.from_thin_walled(sec, lam, elements_along=4, rows_per_ply=rows_per_ply, materials={PLY: MAT})


def _layered_tau(V):
    bounds = -2 * T + T * np.arange(5)
    Ez = np.array([135e9 if a == 0 else 10e9 for a in ANGLES])
    z_n = sum(Ez * T * 0.5 * (bounds[:-1] + bounds[1:])) / sum(Ez * T)
    EI = sum(Ez * ((bounds[1:] - z_n) ** 3 - (bounds[:-1] - z_n) ** 3) / 3)

    def tau(z):
        k = np.searchsorted(bounds, z + 1e-15, side="right") - 1
        top = sum(Ez[j] * ((bounds[j + 1] - z_n) ** 2 - (max(bounds[j], z) - z_n) ** 2) / 2 for j in range(k, 4))
        return V * top / (EI * B)

    return bounds, tau


# --------------------------------------------------------------- interfaces
def test_ply_interface_shear_is_exact_with_two_rows_per_ply():
    V = 10.0
    bounds, tau = _layered_tau(V)
    for rows in (2, 3):
        a = SolidSectionAnalysis(_strip(rows))
        out = a.interface_stresses(forces=[0, 0, 0, 0, V, 0])
        assert len(out) == 3 * 4                        # three interfaces, four elements along
        for _, _, y, z, s in out:
            k = int(np.argmin(np.abs(bounds - z)))
            assert abs(bounds[k] - z) < 1e-12 and 0 < k < 4
            assert s[2] == pytest.approx(tau(bounds[k]), rel=1e-6)
    # the raw edge values of two rows per ply are off by up to 15 %
    a = SolidSectionAnalysis(_strip(2))
    e, f = a.ply_interfaces()[-1]
    raw = a.stress_at(e, 0.0, 1.0, forces=[0, 0, 0, 0, V, 0]).interlaminar[2]
    assert abs(raw / tau(bounds[3]) - 1) > 0.1


def test_one_row_per_ply_falls_back_to_the_element_values():
    a = SolidSectionAnalysis(_strip(1))
    e, f = a.ply_interfaces()[0]
    _, _, s = a.interface_stress(e, f, 0.0, forces=[0, 0, 0, 0, 10.0, 0])
    lo = a.stress_at(e, 0.0, 1.0, forces=[0, 0, 0, 0, 10.0, 0]).interlaminar[2]
    hi = a.stress_at(f, 0.0, -1.0, forces=[0, 0, 0, 0, 10.0, 0]).interlaminar[2]
    assert s[2] == pytest.approx(0.5 * (lo + hi), rel=1e-9)


# --------------------------------------------------------------- jobs
def _run(job_dir, res_dir):
    from workflow_orchestrator.run_job import process_job

    res = common.results_dir(res_dir)
    logging.getLogger().setLevel(logging.WARNING)
    process_job(job_dir, res, {}, {}, force_serial=True, max_processes_per_job=1)
    return res


def _rows(res, *path):
    with open(os.path.join(res, "secondary_results", *path), newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def _box(etype, load, n=6):
    b, h, t = 0.3, 0.12, 2e-3
    pts = {1: (-b / 2, -h / 2), 2: (b / 2, -h / 2), 3: (b / 2, h / 2), 4: (-b / 2, h / 2)}
    walls = [(1, 1, 2, 1), (2, 2, 3, 1), (3, 3, 4, 1), (4, 4, 1, 1)]
    mat = {1: common.Material(70e9, 70e9, 26.9e9, 0.3, 26.9e9, 26.9e9, 2700.0)}
    sec = common.SectionSpec(pts, walls, mat, {1: [(1, 0.0, t)]}, n_sub=2)
    return common.CantileverSpec(1.5, n, etype, sec, tip_load=load)


def _write(job, spec, settings=None, stations=None):
    common.write_job(job, spec)
    common.write_supports(job, spec)
    if settings:
        with open(os.path.join(job, "simulation_settings.txt"), "w", encoding="utf-8") as f:
            f.write(settings)
    if stations:
        with open(os.path.join(job, "stress_stations.txt"), "w", encoding="utf-8") as f:
            f.write("[Stress_Station]\n[x] [y] [z]\n" + "".join(" ".join(map(str, s)) + "\n" for s in stations))


def _block(rows, x):
    return np.array([[float(r[k]) for k in ("sigma_xx", "tau_xs", "tau_xn")] for r in rows
                     if abs(float(r["x"]) - x) < 1e-9])


def test_stations_pick_points_along_the_beam(tmp_path):
    spec = _box("LinearTimoshenkoBeamElement3D", [0, 300.0, -800.0, 50.0, 0, 0])
    _write(str(tmp_path / "mid"), spec)
    # stations: the mid-point of element 2 (by X), and a point next to the axis at X = 0.6
    _write(str(tmp_path / "st"), spec, stations=[(0.625,), (0.6, 0.01, -0.02)])
    mid = _rows(_run(str(tmp_path / "mid"), str(tmp_path / "mid_res")), "wall_stresses", "wall_stresses.csv")
    st = _rows(_run(str(tmp_path / "st"), str(tmp_path / "st_res")), "wall_stresses", "wall_stresses.csv")
    assert sorted({round(float(r["x"]), 9) for r in st}) == [0.6, 0.625]
    np.testing.assert_allclose(_block(st, 0.625), _block(mid, 0.625), rtol=1e-12, atol=1e-6)
    xi = {round(float(r["xi"]), 9) for r in st if abs(float(r["x"]) - 0.6) < 1e-9}
    assert xi == {round(2 * (0.6 - 0.5) / 0.25 - 1, 9)}          # element 2 runs from 0.5 to 0.75


def test_corotational_job_writes_the_linear_stresses_at_small_load(tmp_path):
    load = [0, 30.0, -80.0, 5.0, 0, 0]
    nl = ("[Simulation]\n[Type]\nstatic_nonlinear\n\n[Newton]\ntolerance = 1e-6\nrelative_tolerance = 1e-9\n"
          "relative_reference = first_residual\nmax_iterations = 40\ntolerance_delta_u = 1e-12\n\n[Solver]\n"
          "type = direct\ntolerance = 1e-12\n\n[Nonlinear]\nnum_increments = 2\nline_search = false\n")
    _write(str(tmp_path / "lin"), _box("LinearTimoshenkoBeamElement3D", load))
    _write(str(tmp_path / "cr"), _box("CorotationalTimoshenkoBeamElement3D", load), settings=nl)
    lin = _rows(_run(str(tmp_path / "lin"), str(tmp_path / "lin_res")), "wall_stresses", "wall_stresses.csv")
    cr = _rows(_run(str(tmp_path / "cr"), str(tmp_path / "cr_res")), "wall_stresses", "wall_stresses.csv")
    a, b = _block(cr, 0.125), _block(lin, 0.125)
    assert a.shape == b.shape and np.abs(b).max() > 0
    np.testing.assert_allclose(a, b, rtol=0, atol=1e-4 * np.abs(b).max())


def test_eigen_job_writes_stresses_per_mode(tmp_path):
    spec = _box("LinearTimoshenkoBeamElement3D", [0, 0, 0, 0, 0, 0])
    _write(str(tmp_path / "eig"), spec, settings="[Simulation]\n[Type]\neigen\n\n[Eigen]\nenabled = true\nnum_modes = 3\n")
    res = _run(str(tmp_path / "eig"), str(tmp_path / "eig_res"))
    rows = _rows(res, "wall_stresses", "wall_stresses.csv")
    assert sorted({int(r["mode"]) for r in rows}) == [1, 2, 3]
    for k in (1, 2, 3):
        s = np.array([float(r["sigma_xx"]) for r in rows if int(r["mode"]) == k])
        assert np.isfinite(s).all() and np.abs(s).max() > 0


def test_two_d_section_job_writes_ply_interface_stresses(tmp_path):
    """A cross-ply strip cantilever with a 2D section (two rows per ply): the interface shear
    at mid-span is the layered-beam value."""
    from workflow_orchestrator.run_job import process_job

    n, L, P = 8, 0.5, 10.0
    mesh = _strip(2)
    job = tmp_path / "job"
    job.mkdir()
    (job / "grid.txt").write_text("[Grid]\n[node_id] [x] [y] [z]\n"
                                  + "".join(f"{i} {L * i / n:.15g} 0 0\n" for i in range(n + 1)))
    (job / "element.txt").write_text(
        "[Element]\n[element_id] [node1] [node2] [element_type] [axial_order] [bending_y_order] [bending_z_order] "
        "[shear_y_order] [shear_z_order] [torsion_order] [load_order]\n"
        + "".join(f"{e} {e} {e + 1} LinearTimoshenkoBeamElement3D 2 2 2 1 1 2 2\n" for e in range(n)))
    (job / "material.txt").write_text("[Material]\n[element_id] [E] [G] [nu] [rho]\n"
                                      + "".join(f"{e} 7e10 2.6e10 0.3 1550\n" for e in range(n)))
    (job / "section.txt").write_text("[Section]\n[element_id] [A] [I_x] [I_y] [I_z] [J_t] [kappa] [alpha]\n"
                                     + "".join(f"{e} 8e-4 1e-7 1e-7 3e-8 1e-8 0.8333 0.0\n" for e in range(n)))
    rows_n = "".join(f"5 {i + 1} {y:.17g} {z:.17g}\n" for i, (y, z) in enumerate(mesh.nodes))
    rows_e = "".join(f"5 {k + 1} Q8 {' '.join(str(c + 1) for c in conn)} 3 {th:g} {al:g}\n"
                     for k, (conn, th, al) in enumerate(zip(mesh.elements, mesh.theta, mesh.alpha)))
    m = MAT
    (job / "section_mesh_2d.txt").write_text(
        "[Section_Node]\n[section_id] [node_id] [y] [z]\n" + rows_n
        + "[Section_Element]\n[section_id] [element_id] [type] [nodes] [material_id] [theta] [alpha]\n" + rows_e
        + "[Solid_Material]\n[material_id] [E1] [E2] [E3] [G12] [G13] [G23] [nu12] [nu13] [nu23] [rho]\n"
        f"3 {m.E1} {m.E2} {m.E3} {m.G12} {m.G13} {m.G23} {m.nu12} {m.nu13} {m.nu23} 1550\n")
    (job / "element_section.txt").write_text("[Element_Section]\n[element_id] [section_id]\n"
                                             + "".join(f"{e} 5\n" for e in range(n)))
    (job / "prescribed_displacement.txt").write_text(
        "[Prescribed Displacement]\n[id] [node_id] [dof] [value] [type] [comment]\n"
        + "".join(f"{i} 0 {d} 0.0 displacement\n" for i, d in enumerate(("UX", "UY", "UZ", "RX", "RY", "RZ"))))
    (job / "point_load.txt").write_text("[Point load]\n[x] [y] [z] [F_x] [F_y] [F_z] [M_x] [M_y] [M_z]\n"
                                        f"{L} 0 0 0 0 {P} 0 0 0\n")
    (job / "simulation_settings.txt").write_text(
        "[Simulation]\n[Type]\nstatic\n\n[Solver]\ntype = direct\ntolerance = 1e-12\nmax_iterations = 500\n")
    # the middle of element 3: a linear Timoshenko element's shear strain is right at its mid-point
    (job / "stress_stations.txt").write_text("[Stress_Station]\n[x]\n0.21875\n")
    res = common.results_dir(str(tmp_path / "res"))
    logging.getLogger().setLevel(logging.WARNING)
    process_job(str(job), res, {}, {}, force_serial=True, max_processes_per_job=1)
    rows = _rows(res, "section_stresses", "interface_stresses.csv")
    bounds, tau = _layered_tau(P)
    assert len(rows) == 12
    for r in rows:
        assert float(r["x"]) == pytest.approx(0.21875)
        z = float(r["z"])
        k = int(np.argmin(np.abs(bounds - z)))
        # the Timoshenko element's shear force is exact under a tip load
        assert abs(float(r["tau_xn"])) == pytest.approx(tau(bounds[k]), rel=1e-4)

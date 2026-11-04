"""Rotating-frame beams: centrifugal stiffening, spin softening, Coriolis and the propeller moment."""

import importlib
import os

import numpy as np
import pytest
import scipy.sparse as sp

from pre_processing.element_library.rotating_frame import Rotation, section_point_masses
from processing.eigen.rotating_frame import campbell, rotating_frame_vibration


def _cantilever(tmp_path, n, etype="LinearEulerBernoulliBeamElement3D", L=1.0, E=1e9, A=1e-3, Iy=1e-7, Iz=4e-7,
                rho=1000.0):
    from pre_processing.element_library.element_factory import ElementFactory

    cls = getattr(importlib.import_module(ElementFactory.ELEMENT_CLASS_MAP[etype]), etype)
    xyz = np.array([[i * L / n, 0.0, 0.0] for i in range(n + 1)])
    sec = {"A": np.full(n, A), "I_x": np.full(n, Iy + Iz), "I_y": np.full(n, Iy), "I_z": np.full(n, Iz),
           "J_t": np.full(n, Iy + Iz)}
    mat = {"E": np.full(n, E), "G": np.full(n, E / 2.6), "nu": np.full(n, 0.3), "rho": np.full(n, rho)}
    orders = {k: np.full(n, 4, dtype=np.int64)
              for k in ("axial", "bending_y", "bending_z", "shear_y", "shear_z", "torsion", "load")}
    ed = {"ids": np.arange(n), "connectivity": np.array([[i, i + 1] for i in range(n)]),
          "types": np.array([etype] * n), "integration_orders": orders}
    gd = {"coordinates": xyz, "ids": np.arange(n + 1)}
    jd = str(tmp_path / f"job_{n}_{A}_{Iy}_{Iz}")
    for s in ("element_stiffness_matrices", "element_force_vectors"):
        os.makedirs(os.path.join(jd, s), exist_ok=True)
    els = [cls(element_id=i, element_dictionary=ed, grid_dictionary=gd, material_dictionary=mat,
               section_dictionary=sec, point_load_array=np.empty((0, 9)), distributed_load_array=np.empty((0, 9)),
               job_results_dir=jd) for i in range(n)]
    nd = 6 * (n + 1)
    K, M = np.zeros((nd, nd)), np.zeros((nd, nd))
    for e in els:
        d = np.asarray(e.assemble_global_dof_indices())
        K[np.ix_(d, d)] += e.element_stiffness_matrix().K_e
        M[np.ix_(d, d)] += e.element_mass_matrix().M_e
    return els, xyz, ed["connectivity"], sp.csr_matrix(K), sp.csr_matrix(M), np.arange(6)


def _mode(r, component):
    """Angular frequency of the lowest mode dominated by DOF ``component`` (0-5)."""
    for j in range(r.modes.shape[1]):
        if np.abs(r.modes[component::6, j]).max() > 0.5 * np.abs(r.modes[:, j]).max():
            return r.frequencies_hz[j] * 2.0 * np.pi
    raise AssertionError("no such mode")


def test_point_masses_reproduce_the_section_mass():
    mp = {"m": 2.0, "y_m": 0.1, "z_m": -0.05, "I_yy": 0.03, "I_zz": 0.08, "I_yz": 0.01}
    yz, m = section_point_masses(mp)
    assert m.sum() == pytest.approx(2.0)
    assert m @ yz[:, 0] == pytest.approx(2.0 * 0.1)
    assert m @ yz[:, 1] == pytest.approx(2.0 * -0.05)
    assert m @ yz[:, 1] ** 2 == pytest.approx(0.03)
    assert m @ yz[:, 0] ** 2 == pytest.approx(0.08)
    assert m @ (yz[:, 0] * yz[:, 1]) == pytest.approx(0.01)


# Hodges & Rutkowski (1981), uniform cantilever rotating about its root: first flap
# frequency, both nondimensionalised by sqrt(EI / (m L^4))
HODGES_FLAP = {0: 3.5160, 2: 4.1373, 4: 5.5850, 6: 7.3604, 8: 9.2568, 10: 11.2023, 12: 13.1702}


def test_flap_frequencies_match_hodges(tmp_path):
    """Hodges' beam is inextensible: a slender section (I / (A L^2) = 1e-6) keeps the
    centrifugal stretch, Omega_bar^2 I / (A L^2) of the tension, negligible."""
    E, Iy, rho, A = 1e9, 1e-9, 1000.0, 1e-3
    els, xyz, con, K, M, bc = _cantilever(tmp_path, 20, E=E, A=A, Iy=Iy, Iz=100 * Iy)
    s = np.sqrt(E * Iy / (rho * A))
    for Ob, lam in HODGES_FLAP.items():
        r = rotating_frame_vibration(els, xyz, con, K, M, bc, Rotation(Ob * s + 1e-12, (0, 0, 1)), 8,
                                     gyroscopic=False)
        assert _mode(r, 2) / s == pytest.approx(lam, rel=6e-4)


def test_lag_is_flap_less_the_spin(tmp_path):
    """In-plane: omega_lag^2 = omega_flap^2 - Omega^2 for the same bending stiffness."""
    A = 1e-1
    a = _cantilever(tmp_path, 40, A=A, Iy=4e-7, Iz=1e-6)     # flap on EI = 4e-7 E
    b = _cantilever(tmp_path, 40, A=A, Iy=1e-6, Iz=4e-7)     # lag on EI = 4e-7 E
    s = np.sqrt(1e9 * 4e-7 / (1000.0 * A))
    for Ob in (2.0, 6.0, 10.0):
        rot = Rotation(Ob * s, (0, 0, 1))
        wf = _mode(rotating_frame_vibration(*a, rot, 8, gyroscopic=False), 2)
        wl = _mode(rotating_frame_vibration(*b, rot, 8, gyroscopic=False), 1)
        assert wl ** 2 == pytest.approx(wf ** 2 - (Ob * s) ** 2, rel=1e-3)


def test_shaft_spinning_about_its_axis_splits_into_forward_and_backward_whirl(tmp_path):
    """Coriolis terms: in the rotating frame a symmetric shaft has omega_n -+ Omega."""
    sys = _cantilever(tmp_path, 30, A=1e-1, Iy=1e-6, Iz=1e-6)
    w1 = rotating_frame_vibration(*sys, Rotation(1e-9, (1, 0, 0)), 4).frequencies_hz[0]
    for frac in (0.2, 0.5):
        r = rotating_frame_vibration(*sys, Rotation(frac * 2 * np.pi * w1, (1, 0, 0)), 4)
        assert r.frequencies_hz[:2] / w1 == pytest.approx([1 - frac, 1 + frac], rel=5e-4)


def test_propeller_moment(tmp_path):
    """A flat section turning about the beam axis: omega^2 = omega_0^2 + Omega^2 (I_z - I_y) / (I_y + I_z)
    (no centrifugal tension, which adds the tension-torsion stiffening on top)."""
    Iy, Iz = 1e-8, 1e-6
    sys = _cantilever(tmp_path, 30, Iy=Iy, Iz=Iz)
    w0 = _mode(rotating_frame_vibration(*sys, Rotation(1e-9, (0, 0, 1)), 12, gyroscopic=False), 3)
    for Om in (0.5 * w0, w0):
        r = rotating_frame_vibration(*sys, Rotation(Om, (0, 0, 1)), 12, gyroscopic=False, prestress=False)
        assert _mode(r, 3) ** 2 == pytest.approx(w0 ** 2 + Om ** 2 * (Iz - Iy) / (Iz + Iy), rel=1e-6)
        # the centrifugal tension stiffens torsion further
        r = rotating_frame_vibration(*sys, Rotation(Om, (0, 0, 1)), 12, gyroscopic=False)
        assert _mode(r, 3) ** 2 > 1.1 * (w0 ** 2 + Om ** 2 * (Iz - Iy) / (Iz + Iy))


def test_campbell_rows_and_gauss_data(tmp_path):
    sys = _cantilever(tmp_path, 10)
    rows = campbell(*sys, Rotation(1.0, (0, 0, 1)), [0.0, 50.0], 3, gyroscopic=True)
    assert rows.shape == (2, 4) and rows[1, 1] > rows[0, 1]       # flap stiffens
    r = rotating_frame_vibration(*sys, Rotation(50.0, (0, 0, 1)), 3)
    o = r.element_objects[0]
    assert len(o.gauss_data) == 4 and o.gauss_data[0].Phi.shape == (4, 3, 12)
    # centrifugal load: total radial force is m Omega^2 int x dx
    m = 1000.0 * 1e-3
    assert r.f_centrifugal[0::6].sum() == pytest.approx(m * 50.0 ** 2 * 0.5, rel=1e-9)
    assert np.allclose(o.G, -o.G.T) and np.allclose(o.K_soft, o.K_soft.T)


def test_eigen_job_with_rotation(tmp_path):
    """``[Rotation]`` in an eigen job: rotating-frame frequencies, a Campbell table and the Gauss data."""
    import logging
    import shutil

    from workflow_orchestrator.run_job import process_job

    root = os.path.join(os.path.dirname(__file__), "..")
    job = tmp_path / "job"
    shutil.copytree(os.path.join(root, "jobs", "job_smoke_eigen"), job)
    with open(job / "simulation_settings.txt", "a", encoding="utf-8") as f:
        f.write("\n[Rotation]\nomega = 40\naxis = 0, 0, 1\norigin = 0, 0, 0\nomega_sweep = 0, 20, 40\n")
    res = tmp_path / "results"
    for sub in ("element_stiffness_matrices", "element_force_vectors", "primary_results", "secondary_results",
                "logs", "maps", "diagnostics"):
        os.makedirs(res / sub)
    level = logging.getLogger().level
    try:
        process_job(str(job), str(res), {}, {}, force_serial=True, max_processes_per_job=1)
    finally:
        logging.getLogger().setLevel(level)
    rot = res / "primary_results" / "rotation"
    rows = np.loadtxt(rot / "campbell.csv", delimiter=",", skiprows=1)
    assert rows.shape == (3, 5)
    assert rows[2, 1] > rows[1, 1] > rows[0, 1]           # the first (flap) mode stiffens with speed
    assert (rot / "rotating_frame_gauss.csv").is_file()
    U = np.loadtxt(rot / "U_centrifugal.csv", delimiter=",")
    assert U[-6] > 0.0                                       # the tip moves outward


def test_static_job_with_rotation_matches_the_exact_centrifugal_stretch(tmp_path):
    """A bar spinning about z through its root, axial only: u(x) = m Omega^2 x (3 L^2 - x^2) / (6 EA')
    with EA' softened by the spin; at the tip u = m Omega^2 L^3 / (3 EA) to first order."""
    import logging
    import shutil

    from workflow_orchestrator.run_job import process_job

    root = os.path.join(os.path.dirname(__file__), "..")
    job = tmp_path / "job"
    shutil.copytree(os.path.join(root, "jobs", "job_smoke_eigen"), job)
    (job / "point_load.txt").write_text("[Point load]\n[x] [y] [z] [F_x] [F_y] [F_z] [M_x] [M_y] [M_z]\n")
    (job / "simulation_settings.txt").write_text(
        "[Simulation]\n[Type]\nstatic\n\n[Rotation]\nomega = 100\naxis = 0, 0, 1\n")
    res = tmp_path / "results"
    for sub in ("element_stiffness_matrices", "element_force_vectors", "primary_results", "secondary_results",
                "logs", "maps", "diagnostics"):
        os.makedirs(res / sub)
    level = logging.getLogger().level
    try:
        process_job(str(job), str(res), {}, {}, force_serial=True, max_processes_per_job=1)
    finally:
        logging.getLogger().setLevel(level)
    import csv

    with open(res / "primary_results" / "global" / "U_global.csv", newline="", encoding="utf-8") as f:
        U = np.array([float(r[1]) for r in list(csv.reader(f))[1:]])
    E, A, rho, L, Om = 2.1e11, 0.00131, 7850.0, 2.0, 100.0
    k = np.sqrt(rho * Om ** 2 / E)                 # spin-softened bar: u'' + k^2 u = -k^2 x
    tip = (np.sin(k * L) / (k * np.cos(k * L)) - L)
    assert U[4 * 6] == pytest.approx(tip, rel=2e-3)
    assert abs(U[4 * 6 + 1]) < 1e-12 and abs(U[4 * 6 + 2]) < 1e-12

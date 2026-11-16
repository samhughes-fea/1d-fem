"""Rotating-frame beams: centrifugal stiffening, spin softening, Coriolis and the propeller moment."""

import importlib
import os

import numpy as np
import pytest
import scipy.sparse as sp

from onedfem.pre_processing.element_library.rotating_frame import Rotation, section_point_masses
from onedfem.processing.eigen.rotating_frame import campbell, rotating_frame_vibration


def _cantilever(tmp_path, n, etype="LinearEulerBernoulliBeamElement3D", L=1.0, E=1e9, A=1e-3, Iy=1e-7, Iz=4e-7,
                rho=1000.0):
    from onedfem.pre_processing.element_library.element_factory import ElementFactory

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


def _fake_beam(n_dof=12, restrained_warping=False, cross_section_modes=False, extra_dofs=0,
               warping_inertia=0.0, n_warping=0, n_modes=0, n_shear=0, gbt_mass=None,
               curved=False, mass=None):
    """Minimal straight beam for ``element_rotating_frame`` (duck typed)."""
    mp = mass or {"m": 2.0, "y_m": 0.0, "z_m": 0.0, "I_yy": 0.08, "I_zz": 0.02, "I_yz": 0.0}
    I_w = float(warping_inertia)
    Ms = gbt_mass
    n_w, n_m, n_s = n_warping, n_modes, n_shear
    n = n_dof
    warp, gbt, extra, is_curved = restrained_warping, cross_section_modes, extra_dofs, curved

    class _SA:
        mass = mp
        warping_inertia = I_w

        def gbt_matrices(self, n_modes, n_shear=0):
            if Ms is None:
                raise ValueError("no GBT matrices")
            return np.zeros_like(Ms), np.asarray(Ms, dtype=np.float64)

    class _El:
        element_id = 0
        material_array = np.array([0.0, 0.0, 0.0, 1000.0])
        A, I_y, I_z = 1e-3, 1e-7, 4e-7
        bending_interpolation = "linear"
        has_thin_walled_section = True
        section_analysis = _SA()
        _frame_arc = object() if is_curved else None
        _frame_R0 = np.eye(3)
        _frame_R1 = np.eye(3)
        extra_dofs = extra
        restrained_warping = warp
        cross_section_modes = gbt
        n_warping = n_w
        n_modes = n_m
        n_shear = n_s

        def assemble_global_dof_indices(self):
            return np.arange(n)

    return _El()


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


def test_straight_beam_phi_stays_twelve_columns():
    from onedfem.pre_processing.element_library.rotating_frame import element_rotating_frame

    el = _fake_beam(12)
    o = element_rotating_frame(el, [0, 0, 0], [1, 0, 0], Rotation(10.0, (0, 0, 1)))
    assert o.gauss_data[0].Phi.shape == (4, 3, 12)


def test_reddy_slopes_get_no_rotating_frame_terms():
    from onedfem.pre_processing.element_library.rotating_frame import element_rotating_frame

    el = _fake_beam(16, extra_dofs=2, warping_inertia=1e-6)
    o = element_rotating_frame(el, [0, 0, 0], [1, 0, 0], Rotation(10.0, (0, 0, 1)))
    assert o.gauss_data[0].Phi.shape == (4, 3, 16)
    assert np.allclose(o.gauss_data[0].Phi[:, :, 12:], 0.0)
    assert np.allclose(o.K_soft[12:, :], 0.0) and np.allclose(o.G[12:, :], 0.0)


def test_warping_dofs_take_spin_softening_of_the_warping_inertia():
    from onedfem.pre_processing.element_library.rotating_frame import element_rotating_frame

    I_w, L, Om = 4.0e-6, 2.0, 25.0
    el12 = _fake_beam(12)
    el14 = _fake_beam(14, restrained_warping=True, warping_inertia=I_w)
    rot = Rotation(Om, (0, 0, 1))
    a = element_rotating_frame(el12, [0, 0, 0], [L, 0, 0], rot)
    b = element_rotating_frame(el14, [0, 0, 0], [L, 0, 0], rot)
    assert b.gauss_data[0].Phi.shape == (4, 3, 14)
    assert np.allclose(b.K_soft[:12, :12], a.K_soft, atol=1e-12 * max(np.abs(a.K_soft).max(), 1.0))
    psi = np.zeros(14)
    psi[12:] = 1.0
    # constant ψ: u_p = ω_p e_x, so u^T K_soft u = I_ω L (e_x · W² · e_x) = -Ω² I_ω L
    assert psi @ b.K_soft @ psi == pytest.approx(-Om ** 2 * I_w * L, rel=1e-9)
    assert np.allclose(b.G, -b.G.T)
    assert abs(psi @ b.G @ psi) < 1e-18
    # extra columns are axial: orthogonal to a rigid translation, so no 12-block leak
    assert np.linalg.norm(b.K_soft[:12, 12:]) < 1e-12 * max(np.abs(b.K_soft).max(), 1.0)


def test_gbt_extras_take_radial_and_warping_inertia():
    from onedfem.pre_processing.element_library.rotating_frame import element_rotating_frame

    I_w, I_a, I_ap, L, Om = 3.0e-6, 0.5, 1.0e-6, 1.0, 10.0
    Ms = np.diag([2.0, 2.0, 2.0, 0.1, 0.08, 0.02, I_w, I_a, I_ap])
    el = _fake_beam(18, restrained_warping=True, cross_section_modes=True, n_warping=1, n_modes=1,
                    gbt_mass=Ms, warping_inertia=I_w)
    o_z = element_rotating_frame(el, [0, 0, 0], [L, 0, 0], Rotation(Om, (0, 0, 1)))
    assert o_z.gauss_data[0].Phi.shape == (4, 3, 18)
    b = np.zeros(18)
    b[12] = b[12 + 3] = 1.0
    assert b @ o_z.K_soft @ b == pytest.approx(-Om ** 2 * I_w * L, rel=1e-9)
    o_x = element_rotating_frame(el, [0, 0, 0], [L, 0, 0], Rotation(Om, (1, 0, 0)))
    a = np.zeros(18)
    a[13] = a[13 + 3] = 1.0   # a_k at both nodes, a' = 0
    assert a @ o_x.K_soft @ a == pytest.approx(-Om ** 2 * I_a * L, rel=1e-9)


def test_curved_elements_still_raise():
    from onedfem.pre_processing.element_library.rotating_frame import element_rotating_frame

    el = _fake_beam(14, restrained_warping=True, warping_inertia=1e-6, curved=True)
    with pytest.raises(NotImplementedError, match="curved"):
        element_rotating_frame(el, [0, 0, 0], [1, 0, 0], Rotation(1.0, (0, 0, 1)))


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

    from onedfem.workflow_orchestrator.run_job import process_job

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

    from onedfem.workflow_orchestrator.run_job import process_job

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

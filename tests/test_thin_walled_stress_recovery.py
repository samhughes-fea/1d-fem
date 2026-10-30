"""Wall, ply and interlaminar stress recovery for thin-walled sections."""

from __future__ import annotations

import importlib
import importlib.util
import os
import shutil
import tempfile

import numpy as np
import pytest

from pre_processing.element_library.dof_layout import mesh_dof_per_node
from pre_processing.section_integration.section_analysis import SectionAnalysis
from pre_processing.section_integration.thin_walled import (
    ThinWalledSection,
    ThinWalledSectionAnalysis,
    Wall,
    WallLaminate,
)
from pre_processing.section_integration.thin_walled.element_recovery import element_wall_stress
from pre_processing.section_integration.thin_walled.laminate import Ply, PlyMaterial
from pre_processing.section_integration.thin_walled.stress_recovery import (
    SaintVenantRecovery,
    mid_surface_strains,
    through_thickness,
    wall_resultants,
)
from test_gbt_cross_section_modes import BS, E, G, GBT, NU, _ed, box, i_beam

CFRP = PlyMaterial(E1=140e9, E2=10e9, G12=5e9, nu12=0.3, G13=5e9, G23=3.5e9, rho=1600.0)


def laminate(angles, t_ply=0.25e-3):
    return WallLaminate([Ply(CFRP, a, t_ply) for a in angles])


def composite_box(b=0.4, h=0.2, angles=(0, 45, -45, 90), n_sub=1, hoop="restrained"):
    pts = {1: (-b / 2, -h / 2), 2: (b / 2, -h / 2), 3: (b / 2, h / 2), 4: (-b / 2, h / 2)}
    walls = [Wall(1, 1, 2, 1), Wall(2, 2, 3, 1), Wall(3, 3, 4, 1), Wall(4, 4, 1, 1)]
    return ThinWalledSectionAnalysis(ThinWalledSection(pts, walls), {1: laminate(angles)}, hoop=hoop, n_sub=n_sub)


# --------------------------------------------------------------------------- #
# Saint-Venant recovery against hand values
# --------------------------------------------------------------------------- #
B_, H_, T_ = 0.4, 0.2, 2e-3
I_Y = B_ * T_ * H_**2 / 2 + T_ * H_**3 / 6  # thin-wall box about y


def test_bending_stress_is_m_z_over_i():
    an = box(B_, H_, T_)
    # free hoop: the wall bends in uniaxial stress (restrained hoop adds the plate factor 1/(1 - nu^2))
    an = ThinWalledSectionAnalysis(an.section, an.laminates, hoop="free")
    rec = SaintVenantRecovery(an)
    M = 1e3
    p = rec.at(3, B_ / 4, [0, M, 0, 0, 0, 0])  # top flange (z = h/2)
    assert p.resultants[0] == pytest.approx(M * (H_ / 2) / I_Y * T_, rel=1e-4)
    for q in p.plies:
        assert q.stress_wall[0] == pytest.approx(M * (H_ / 2 - q.n) / I_Y, rel=2e-4)  # wall 3 runs -y, so n = -z
    assert abs(p.resultants[1]) < 1e-9 * abs(p.resultants[0])  # N_ss = 0


def test_torsion_shear_flow_is_bredt():
    rec = SaintVenantRecovery(box(B_, H_, T_))
    T = 100.0
    for wall in (1, 2, 3, 4):
        assert rec.at(wall, 0.05, [0, 0, 0, 0, 0, T]).resultants[2] == pytest.approx(T / (2 * B_ * H_), rel=2e-4)


def test_shear_flow_at_the_web_middle():
    rec = SaintVenantRecovery(box(B_, H_, T_))
    V = 1e3
    q = rec.at(2, H_ / 2, [0, 0, 0, 0, V, 0]).resultants[2]
    assert q == pytest.approx(V * T_ / I_Y * (B_ * H_ / 4 + H_**2 / 8), rel=1e-4)


# --------------------------------------------------------------------------- #
# Wall and ply level
# --------------------------------------------------------------------------- #
def test_released_resultants_are_zero():
    lam = laminate((0, 30, -60, 90, 0))  # unsymmetric
    k4 = np.array([1e-4, 2e-4, 0.3, -0.1])
    e = mid_surface_strains(lam, k4, "restrained")
    R = wall_resultants(lam, e)
    assert abs(R[1]) < 1e-10 * np.abs(R).max() and e[4] == 0.0
    e = mid_surface_strains(lam, k4, "free")
    R = wall_resultants(lam, e)
    assert abs(R[1]) < 1e-10 * np.abs(R).max() and abs(R[4]) < 1e-10 * np.abs(R[3:]).max()
    e5 = mid_surface_strains(lam, [1e-4, 2e-4, 0.3, 0.05, -0.1])
    assert e5[4] == 0.05 and abs(wall_resultants(lam, e5)[1]) < 1e-10 * np.abs(R).max()


def test_ply_stresses_integrate_to_the_wall_resultants():
    lam = laminate((0, 30, -60, 90, 0))
    e = mid_surface_strains(lam, [1e-4, 2e-4, 0.3, -0.1])
    R = wall_resultants(lam, e)
    pts = through_thickness(lam, e, points_per_ply=2)  # stresses are linear in each ply
    N = np.zeros(3)
    M = np.zeros(3)
    for a, b in zip(pts[::2], pts[1::2]):
        dn = b.n - a.n
        N += 0.5 * dn * (a.stress_wall + b.stress_wall)
        M += dn / 6 * ((2 * a.n + b.n) * a.stress_wall + (a.n + 2 * b.n) * b.stress_wall)
    np.testing.assert_allclose(N, R[:3], rtol=1e-10, atol=1e-10 * np.abs(R[:3]).max())
    np.testing.assert_allclose(M, R[3:], rtol=1e-10, atol=1e-10 * np.abs(R[3:]).max())


def test_ply_axes_of_a_90_degree_ply():
    lam = laminate((90,))
    e = np.array([1e-4, -3e-5, 2e-4, 0, 0, 0])
    q = through_thickness(lam, e)[0]
    np.testing.assert_allclose(q.strain_ply, [e[1], e[0], -e[2]], atol=1e-18)
    np.testing.assert_allclose(q.stress_ply, [q.stress_wall[1], q.stress_wall[0], -q.stress_wall[2]], rtol=1e-12)


# --------------------------------------------------------------------------- #
# Interlaminar stresses
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("forces", [[0, 0, 0, 0, 1e3, 0], [0, 0, 0, 1e3, 0, 0], [0, 0, 0, 0, 0, 100.0]])
@pytest.mark.parametrize("wall, s", [(1, 0.1), (2, 0.07), (3, 0.3)])
def test_interlaminar_shear_is_free_at_the_top_and_integrates_to_q(forces, wall, s):
    rec = SaintVenantRecovery(composite_box(angles=(0, 45, -45, 90, 0, 30)))
    p = rec.at(wall, s, forces, points_per_ply=5)
    t_xn = np.array([q.tau_xn for q in p.plies])
    n = np.array([q.n for q in p.plies])
    scale = max(np.abs(t_xn).max(), 1.0)
    assert abs(t_xn[0]) == 0.0
    assert abs(t_xn[-1]) < 1e-6 * scale  # wall membrane equilibrium dN_xx/dx + dN_xs/ds = 0
    # int tau_xn dn = Q_x (integrating by parts with free faces); Simpson is exact per ply
    integral = sum((b - a) / 6 * (fa + 4 * fm + fb) for a, b, fa, fm, fb in
                   zip(n[0::5], n[4::5], t_xn[0::5], t_xn[2::5], t_xn[4::5]))
    assert integral == pytest.approx(p.transverse_shear[0], rel=1e-6, abs=1e-9 * scale * 1e-3)


def test_homogeneous_wall_shear_is_parabolic():
    rec = SaintVenantRecovery(box(B_, H_, T_))
    p = rec.at(3, B_ / 4, [0, 0, 0, 0, 1e3, 0])
    mid = [q for q in p.plies if abs(q.n) < 1e-12][0]
    assert mid.tau_xn == pytest.approx(1.5 * p.transverse_shear[0] / T_, rel=1e-8)


# --------------------------------------------------------------------------- #
# Element recovery
# --------------------------------------------------------------------------- #
def build_and_solve(etype, sa, L, n_el, tip_load, modes=0, shear=0, root_free=()):
    """Cantilever by direct assembly; returns (elements, nodal U (n_nodes, stride))."""
    from pre_processing.element_library.element_factory import ElementFactory

    Cls = getattr(importlib.import_module(ElementFactory.ELEMENT_CLASS_MAP[etype]), etype)
    ed = _ed(etype, n_el, modes, sa, shear)
    if etype == "LinearEulerBernoulliBeamElement3D":
        ed["warping"] = np.zeros(n_el, dtype=np.int8)
    ed["section_analysis_per_element"] = tuple([sa] * n_el)
    x = np.linspace(0.0, L, n_el + 1)
    coords = np.column_stack([x, 0 * x, 0 * x])
    jr = os.path.join(tempfile.mkdtemp(), "job")
    for sub in ("element_stiffness_matrices", "element_force_vectors", "logs"):
        os.makedirs(os.path.join(jr, sub), exist_ok=True)
    try:
        mat = {k: np.full(n_el, v) for k, v in dict(E=E, G=G, nu=NU, rho=0.0).items()}
        sec = {k: np.full(n_el, v) for k, v in dict(A=1e-3, I_x=0.0, I_y=1e-8, I_z=1e-8, J_t=1e-9, kappa=5 / 6,
                                                    alpha=0.0, y_sc=0.0, z_sc=0.0, Gamma=1e-9).items()}
        if etype == "LinearEulerBernoulliBeamElement3D":
            sec = {k: sec[k] for k in ("A", "I_x", "I_y", "I_z", "J_t")}
        dpn = mesh_dof_per_node(ed)
        n = dpn * (n_el + 1)
        K = np.zeros((n, n))
        els = []
        for e in range(n_el):
            el = Cls(element_id=e, element_dictionary=ed, grid_dictionary={"coordinates": coords},
                     material_dictionary=mat, section_dictionary=sec,
                     point_load_array=np.empty((0, 9)), distributed_load_array=np.empty((0, 9)), job_results_dir=jr)
            el.logger_operator = None
            g = el.assemble_global_dof_indices()
            K[np.ix_(g, g)] += el.element_stiffness_matrix().K_e
            els.append((el, g))
        F = np.zeros(n)
        for local, val in tip_load.items():
            F[dpn * n_el + local] = val
        fixed = set(range(dpn)) - set(root_free)
        free = np.array(sorted(set(range(n)) - fixed))
        U = np.zeros(n)
        U[free] = np.linalg.solve(K[np.ix_(free, free)], F[free])
        return [(el, U[g]) for el, g in els]
    finally:
        shutil.rmtree(os.path.dirname(jr), ignore_errors=True)


@pytest.mark.parametrize("etype", ["LinearEulerBernoulliBeamElement3D", "LinearTimoshenkoBeamElement3D"])
def test_beam_element_recovery_matches_exact_section_forces(etype):
    an = composite_box()
    sa = SectionAnalysis.from_thin_walled(an)
    L, n_el, P = 2.0, 8, 1e3
    els = build_and_solve(etype, sa, L, n_el, {2: P})
    el, U = els[3]
    x = 0.5 * (el.x_start + el.x_end)
    exact = SaintVenantRecovery(an).at(3, 0.1, [0, P * (x - L), 0, 0, P, 0])  # M_y' = V_z, M_y(L) = 0
    got = element_wall_stress(el, U, 0.0, an, 3, 0.1)
    np.testing.assert_allclose(got.resultants[[0, 2]], exact.resultants[[0, 2]], rtol=2e-3)
    got_t = [q.tau_xn for q in got.plies]
    np.testing.assert_allclose(got_t, [q.tau_xn for q in exact.plies], rtol=1e-6, atol=1e-6 * np.abs(got_t).max())


def test_benscoter_with_free_warping_is_saint_venant_and_restraint_adds_axial_stress():
    an = i_beam()
    sa = SectionAnalysis.from_thin_walled(an)
    L, T = 3.0, 50.0
    free = build_and_solve(BS, sa, L, 12, {3: T}, root_free=(6,))
    el, U = free[2]
    got = element_wall_stress(el, U, 0.0, an, 1, 0.05)
    exact = SaintVenantRecovery(an).at(1, 0.05, [0, 0, 0, 0, 0, T])
    # open section: the torque is carried by wall twisting moments M_xs
    np.testing.assert_allclose(got.resultants[[2, 5]], exact.resultants[[2, 5]], rtol=1e-6, atol=1e-6 * abs(exact.resultants[5]))
    assert abs(got.resultants[0]) < 1e-6 * T / 0.2
    fixed = build_and_solve(BS, sa, L, 12, {3: T})
    el, U = fixed[0]
    flange_tip = element_wall_stress(el, U, 0.0, an, 1, 0.0)
    assert abs(flange_tip.resultants[0]) > 1e3 * abs(got.resultants[0]) + 1.0  # bimoment stresses


def test_gbt_with_free_warping_matches_saint_venant():
    an = composite_box(n_sub=2)
    sa = SectionAnalysis.from_thin_walled(an, gbt=True)
    n_sv = sa.gbt_sv_modes
    L, n_el, P = 2.0, 32, 1e3
    els = build_and_solve(GBT, sa, L, n_el, {2: P}, root_free=tuple(range(6, 6 + n_sv)))
    el, U = els[16]
    x = 0.5 * (el.x_start + el.x_end)
    exact = SaintVenantRecovery(an).at(3, 0.1, [0, P * (x - L), 0, 0, P, 0])
    got = element_wall_stress(el, U, 0.0, an, 3, 0.1)
    np.testing.assert_allclose(got.resultants[[0, 2]], exact.resultants[[0, 2]], rtol=5e-3)


def test_gbt_recovery_carries_the_hoop_curvature_of_distortion():
    an = i_beam()
    sa = SectionAnalysis.from_thin_walled(an, gbt_modes=2)
    els = build_and_solve(GBT, sa, 3.0, 12, {3: 50.0}, modes=2)
    el, U = els[6]
    p = element_wall_stress(el, U, 0.0, an, 5, 0.1)  # web
    assert abs(p.strains[4]) > 0  # kappa_ss from the distortional modes
    assert abs(p.resultants[1]) < 1e-9 * np.abs(p.resultants).max()  # N_ss = 0


def test_process_job_writes_wall_stresses(tmp_path):
    """End to end: a Timoshenko cantilever on a composite box writes wall_stresses.csv."""
    import csv

    from workflow_orchestrator.run_job import process_job, setup_job_results_directory

    L, n_el, P = 2.0, 4, 1e3
    b, h = 0.4, 0.2
    job = tmp_path / "job_ws"
    job.mkdir()
    TS = "LinearTimoshenkoBeamElement3D"
    (job / "grid.txt").write_text(
        "[Grid]\n[node_id] [x] [y] [z]\n" + "\n".join(f"{i} {L * i / n_el} 0 0" for i in range(n_el + 1)) + "\n",
        encoding="utf-8")
    (job / "element.txt").write_text(
        "[Element]\n[element_id] [node1] [node2] [element_type] [axial_order] [bending_y_order] "
        "[bending_z_order] [shear_y_order] [shear_z_order] [torsion_order] [load_order] [warping]\n"
        + "\n".join(f"{e} {e} {e + 1} {TS} 2 2 2 2 2 2 2 0" for e in range(n_el)) + "\n", encoding="utf-8")
    (job / "material.txt").write_text(
        "[Material]\n[element_id] [E] [G] [nu] [rho]\n" + "\n".join(f"{e} 7e10 2.7e10 0.3 0" for e in range(n_el)) + "\n",
        encoding="utf-8")
    (job / "section.txt").write_text(
        "[Section]\n[element_id] [A] [I_x] [I_y] [I_z] [J_t] [kappa] [alpha]\n"
        + "\n".join(f"{e} 1e-3 0.0 1e-8 1e-8 1e-9 0.8333 0.0" for e in range(n_el)) + "\n", encoding="utf-8")
    (job / "prescribed_displacement.txt").write_text(
        "[Prescribed Displacement]\n[id] [node_id] [dof] [value] [type] [comment]\n"
        + "\n".join(f"{i} 0 {d} 0.0 displacement" for i, d in enumerate(["UX", "UY", "UZ", "RX", "RY", "RZ"])) + "\n",
        encoding="utf-8")
    (job / "point_load.txt").write_text(
        "[Point load]\n[x] [y] [z] [F_x] [F_y] [F_z] [M_x] [M_y] [M_z]\n" f"{L} 0 0 0 0 {P} 0 0 0\n", encoding="utf-8")
    (job / "simulation_settings.txt").write_text(
        "[Simulation]\n[Type]\nstatic\n\n[Solver]\ntype = direct\ntolerance = 1e-12\nmax_iterations = 500\n",
        encoding="utf-8")
    (job / "laminate.txt").write_text(
        "[Ply_Material]\n[material_id] [E1] [E2] [G12] [nu12] [G13] [G23] [rho]\n"
        "1 140e9 10e9 5e9 0.3 5e9 3.5e9 1600\n"
        "[Laminate]\n[laminate_id] [ply] [material_id] [angle_deg] [thickness]\n"
        "1 1 1 0 0.25e-3\n1 2 1 45 0.25e-3\n1 3 1 -45 0.25e-3\n1 4 1 90 0.25e-3\n", encoding="utf-8")
    (job / "section_geometry.txt").write_text(
        "[Section_Point]\n[section_id] [point_id] [y] [z]\n"
        f"1 1 {-b / 2} {-h / 2}\n1 2 {b / 2} {-h / 2}\n1 3 {b / 2} {h / 2}\n1 4 {-b / 2} {h / 2}\n"
        "[Section_Wall]\n[section_id] [wall_id] [point_start] [point_end] [laminate_id]\n"
        "1 1 1 2 1\n1 2 2 3 1\n1 3 3 4 1\n1 4 4 1 1\n", encoding="utf-8")
    (job / "element_section.txt").write_text(
        "[Element_Section]\n[element_id] [section_id]\n" + "\n".join(f"{e} 1" for e in range(n_el)) + "\n",
        encoding="utf-8")
    res = setup_job_results_directory("test_wall_stress_job")
    process_job(str(job), res, {}, {}, force_serial=True, max_processes_per_job=1)
    path = os.path.join(res, "secondary_results", "wall_stresses", "wall_stresses.csv")
    with open(path, newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    assert len(rows) == n_el * 4 * 5 * 4 * 2  # elements x walls x points x plies x faces
    # top flange, element 1 mid-point, centre of the wall: N_xx against the section forces
    an = composite_box(hoop="free")   # the job default
    x = 0.75
    exact = SaintVenantRecovery(an).at(3, b / 2, [0, P * (x - L), 0, 0, P, 0])
    r = [r for r in rows if r["element_id"] == "1" and r["wall_id"] == "3" and abs(float(r["s"]) - b / 2) < 1e-12][0]
    assert float(r["N_xx"]) == pytest.approx(exact.resultants[0], rel=1e-6)
    # N_xs is zero here by symmetry: compare on the scale of N_xx
    assert float(r["N_xs"]) == pytest.approx(exact.resultants[2], rel=1e-5, abs=1e-9 * abs(exact.resultants[0]))


@pytest.mark.parametrize("example", ["box_spar", "three_cell_blade"])
def test_thin_walled_examples_run(example, tmp_path):
    """The worked examples run end to end on a coarse mesh and pass their own checks."""
    path = os.path.join(os.path.dirname(__file__), "..", "examples", "thin_walled", example, "run.py")
    spec = importlib.util.spec_from_file_location(f"example_{example}", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    s = mod.main(n_el=4, plot=False, work_dir=str(tmp_path))
    for q in s["torque_shear_flow"].values():
        assert q["1d-fem"] == pytest.approx(q["bredt"], rel=5e-3)
    assert s["tip"]["Timoshenko"]["twist"] == pytest.approx(s["exact_saint_venant_tip"]["twist"], rel=1e-6)
    assert s["tip"]["Timoshenko"]["w"] == pytest.approx(s["exact_saint_venant_tip"]["w"], rel=2e-2)
    a, b = s["csv_vs_direct_N_xx"]
    assert a == pytest.approx(b, rel=1e-9)

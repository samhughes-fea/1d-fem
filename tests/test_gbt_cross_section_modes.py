"""Cross-section distortional modes (conventional GBT) and the GBT beam element."""

from __future__ import annotations

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
from pre_processing.section_integration.thin_walled.cross_section_modes import CrossSectionModes

E, NU = 70e9, 0.3
G = E / (2 * (1 + NU))
GBT = "LinearGBTBeamElement3D"
BS = "LinearBenscoterBeamElement3D"


def iso(t, rho=0.0):
    return WallLaminate.isotropic(E, NU, t, rho=rho)


def box(b=0.3, h=0.3, t=2e-3, rho=0.0, n_sub=1):
    pts = {1: (-b / 2, -h / 2), 2: (b / 2, -h / 2), 3: (b / 2, h / 2), 4: (-b / 2, h / 2)}
    walls = [Wall(1, 1, 2, 1), Wall(2, 2, 3, 1), Wall(3, 3, 4, 1), Wall(4, 4, 1, 1)]
    return ThinWalledSectionAnalysis(ThinWalledSection(pts, walls), {1: iso(t, rho)}, n_sub=n_sub)


def i_beam(b=0.2, h=0.3, t=2e-3, rho=0.0, n_sub=1):
    pts = {1: (-b / 2, h / 2), 2: (0, h / 2), 3: (b / 2, h / 2),
           4: (-b / 2, -h / 2), 5: (0, -h / 2), 6: (b / 2, -h / 2)}
    walls = [Wall(1, 1, 2, 1), Wall(2, 2, 3, 1), Wall(3, 4, 5, 1), Wall(4, 5, 6, 1), Wall(5, 5, 2, 1)]
    return ThinWalledSectionAnalysis(ThinWalledSection(pts, walls), {1: iso(t, rho)}, n_sub=n_sub)


# --------------------------------------------------------------------------- #
# Modes
# --------------------------------------------------------------------------- #


def test_square_box_distortion_mode_matches_hand_frame_analysis():
    a, t = 0.3, 2e-3
    cm = CrossSectionModes(box(a, a, t))
    assert cm.n_rigid == 3  # closed: no Vlasov torsion mode
    md = cm.mode(0)
    # corners move along the diagonals, joints do not rotate (symmetry)
    np.testing.assert_allclose(np.abs(md["dy"]), 1.0, rtol=1e-6)
    np.testing.assert_allclose(np.abs(md["dz"]), 1.0, rtol=1e-6)
    np.testing.assert_allclose(md["beta"], 0.0, atol=1e-6)
    np.testing.assert_allclose(np.abs(md["u"]), a / 2, rtol=1e-6)  # du/ds = -v with |v| = 1
    fm = cm.frame_matrices(1)
    D_plate = E * t**3 / (12 * (1 - NU**2))
    # four walls, each a fixed-fixed beam with end offset 2: 4 * 12 D (2)^2 / a^3
    assert fm["B"][0, 0] == pytest.approx(192 * D_plate / a**3, rel=1e-6)
    # membrane warping: 4 walls of linear u between -a/2 and a/2, plus the tiny plate term
    membrane = E * t * a**3 / 3
    w2 = 4 * D_plate * a * 17 / 35  # w = -1 + 2 (3 xi^2 - 2 xi^3) on each wall: int w^2 = 17 / 35
    assert fm["C"][0, 0] == pytest.approx(membrane + w2, rel=1e-6)


def test_open_section_has_four_rigid_modes_and_orthogonal_modes():
    cm = CrossSectionModes(i_beam(n_sub=2))
    assert cm.n_rigid == 4
    m = cm.n_modes
    assert m > 4
    fm = cm.frame_matrices(m)
    Cn = fm["C"] / np.sqrt(np.outer(np.diag(fm["C"]), np.diag(fm["C"])))
    np.testing.assert_allclose(Cn, np.eye(m), atol=1e-10)
    Bn = fm["B"] / np.sqrt(np.outer(np.diag(fm["B"]), np.diag(fm["B"])))
    np.testing.assert_allclose(Bn, np.eye(m), atol=1e-10)
    np.testing.assert_allclose(np.diag(fm["B"]) / np.diag(fm["C"]), cm.eigenvalues, rtol=1e-8)
    assert np.all(np.diff(cm.eigenvalues) >= 0)


def test_gbt_section_stiffness_extends_benscoter():
    an = box(0.3, 0.15)
    cm = CrossSectionModes(an)
    K = cm.stiffness(2)
    assert K.shape == (14, 14)
    np.testing.assert_allclose(K[:8, :8], an.benscoter_stiffness)
    np.testing.assert_allclose(K, K.T)
    assert np.linalg.eigvalsh(K).min() > -1e-9 * np.abs(K).max()
    sa = SectionAnalysis.from_thin_walled(an, gbt_modes=2)
    K1, M1 = sa.gbt_matrices(1)
    assert K1.shape == (11, 11) and M1.shape == (9, 9)
    np.testing.assert_allclose(K1[8:, 8:], sa.gbt_stiffness[np.ix_([8, 10, 12], [8, 10, 12])])


def test_gbt_mass_rigid_block_matches_section_mass():
    an = box(0.3, 0.15, rho=1600.0)
    sa = SectionAnalysis.from_thin_walled(an, gbt_modes=1)
    np.testing.assert_allclose(sa.gbt_mass[:6, :6], sa.mass_matrix(), rtol=1e-8, atol=1e-10 * sa.mass["m"])
    assert sa.gbt_mass[6, 6] == pytest.approx(sa.warping_inertia, rel=1e-10)


def test_modes_need_restrained_hoop():
    an = box()
    with pytest.raises(ValueError, match="restrained"):
        CrossSectionModes(ThinWalledSectionAnalysis(an.section, an.laminates, hoop="free"))


# --------------------------------------------------------------------------- #
# Element
# --------------------------------------------------------------------------- #


def _ed(etype, n_el, modes):
    o = {"axial": 2, "bending_y": 2, "bending_z": 2, "shear_y": 2, "shear_z": 2, "torsion": 2, "load": 2}
    ed = {
        "ids": np.arange(n_el),
        "connectivity": np.column_stack([np.arange(n_el), np.arange(1, n_el + 1)]),
        "types": np.array([etype] * n_el),
        "warping": np.ones(n_el, dtype=np.int8),
        "integration_orders": {k: np.full(n_el, v) for k, v in o.items()},
    }
    if modes:
        ed["gbt_modes"] = np.full(n_el, modes)
        ed["extra_dofs"] = np.full(n_el, 1 + 2 * modes)
    return ed


def solve_cantilever(etype, sa, L, n_el, modes, tip_load, fixed_local=(), free_only=None):
    """Direct assembly; root clamped with all extra DOFs fixed. Returns nodal U (n_nodes, stride)."""
    from pre_processing.element_library.element_factory import ElementFactory
    import importlib

    Cls = getattr(importlib.import_module(ElementFactory.ELEMENT_CLASS_MAP[etype]), etype)
    ed = _ed(etype, n_el, modes)
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
        dpn = mesh_dof_per_node(ed)
        n = dpn * (n_el + 1)
        K = np.zeros((n, n))
        for e in range(n_el):
            el = Cls(element_id=e, element_dictionary=ed, grid_dictionary={"coordinates": coords},
                     material_dictionary=mat, section_dictionary=sec,
                     point_load_array=np.empty((0, 9)), distributed_load_array=np.empty((0, 9)), job_results_dir=jr)
            el.logger_operator = None
            g = el.assemble_global_dof_indices()
            K[np.ix_(g, g)] += el.element_stiffness_matrix().K_e
        F = np.zeros(n)
        for local, val in tip_load.items():
            F[dpn * n_el + local] = val
        fixed = set(range(dpn))
        for i in range(n_el + 1):
            for local in fixed_local:
                fixed.add(dpn * i + local)
            if free_only is not None:
                fixed |= {dpn * i + j for j in range(dpn) if j not in free_only}
        free = np.array(sorted(set(range(n)) - fixed))
        U = np.zeros(n)
        U[free] = np.linalg.solve(K[np.ix_(free, free)], F[free])
        return U.reshape(n_el + 1, dpn)
    finally:
        shutil.rmtree(os.path.dirname(jr), ignore_errors=True)


def exact_mode_cantilever(C, D, B, K20, P, L, x):
    """a(x) of C a'''' - (D - 2 K20) a'' + B a = 0, a(0) = a'(0) = 0, C a'' + K20 a = 0 and
    -C a''' + (D - K20) a' = P at x = L."""
    r = np.roots([C, 0, -(D - 2 * K20), 0, B]).astype(complex)

    def rows(xx, k):
        return np.array([ri**k * np.exp(ri * xx) for ri in r])

    A = np.array([rows(0, 0), rows(0, 1), C * rows(L, 2) + K20 * rows(L, 0), -C * rows(L, 3) + (D - K20) * rows(L, 1)])
    c = np.linalg.solve(A, np.array([0, 0, 0, P], dtype=complex))
    return np.real(np.array([rows(xi, 0) @ c for xi in x]))


def test_single_mode_cantilever_matches_exact_gbt_solution():
    an = box(0.3, 0.15)
    sa = SectionAnalysis.from_thin_walled(an, gbt_modes=1)
    K, _ = sa.gbt_matrices(1)
    C, D, B, K20 = K[8, 8], K[9, 9], K[10, 10], K[8, 10]
    assert abs(K[8, 9]) < 1e-9 * C and abs(K[9, 10]) < 1e-9 * np.sqrt(D * B)
    k = (B / C) ** 0.25
    L, P = 4.0 / k, 1.0e3
    x = np.linspace(0, L, 9)
    exact = exact_mode_cantilever(C, D, B, K20, P, L, x)
    # only the mode DOFs (a, a' = local 7, 8) are free
    errs = []
    for n_el in (8, 32):
        U = solve_cantilever(GBT, sa, L, n_el, 1, {7: P}, free_only=(7, 8))
        a = U[:, 7]
        errs.append(abs(a[-1] - exact[-1]) / abs(exact[-1]))
        if n_el == 32:
            np.testing.assert_allclose(a[::4], exact, atol=1e-6 * np.abs(exact).max())
    assert errs[1] < 1e-6 and errs[1] < errs[0]


def test_gbt_with_modes_fixed_is_benscoter_and_freeing_them_softens():
    sa = SectionAnalysis.from_thin_walled(i_beam(), gbt_modes=2)
    L, T = 3.0, 50.0
    U_b = solve_cantilever(BS, sa, L, 12, 0, {3: T})
    U_fixed = solve_cantilever(GBT, sa, L, 12, 2, {3: T}, fixed_local=(7, 8, 9, 10))
    np.testing.assert_allclose(U_fixed[:, :7], U_b, rtol=1e-10, atol=1e-14)
    U_free = solve_cantilever(GBT, sa, L, 12, 2, {3: T}, fixed_local=())
    assert U_free[-1, 3] >= U_b[-1, 3] * (1 - 1e-12)


def test_gbt_element_needs_modes():
    sa = SectionAnalysis.from_thin_walled(box(0.3, 0.15), gbt_modes=1)
    with pytest.raises(ValueError, match="mode"):
        solve_cantilever(GBT, sa, 1.0, 1, 0, {})


def test_gbt_mass_matrix_is_consistent():
    sa = SectionAnalysis.from_thin_walled(box(0.3, 0.15, rho=1600.0), gbt_modes=1)
    from pre_processing.element_library.linear.beam.first_order_shear_deformation_theory.benscoter.linear_gbt_3D import (
        LinearGBTBeamElement3D,
    )

    ed = _ed(GBT, 1, 1)
    ed["section_analysis_per_element"] = (sa,)
    jr = os.path.join(tempfile.mkdtemp(), "job")
    for sub in ("element_stiffness_matrices", "element_force_vectors", "logs"):
        os.makedirs(os.path.join(jr, sub), exist_ok=True)
    try:
        el = LinearGBTBeamElement3D(
            element_id=0, element_dictionary=ed, grid_dictionary={"coordinates": np.array([[0.0, 0, 0], [2.0, 0, 0]])},
            material_dictionary={k: np.array([v]) for k, v in dict(E=E, G=G, nu=NU, rho=0.0).items()},
            section_dictionary={k: np.array([v]) for k, v in dict(A=1e-3, I_x=0.0, I_y=1e-8, I_z=1e-8, J_t=1e-9,
                                                                   kappa=5 / 6, alpha=0.0, y_sc=0.0, z_sc=0.0, Gamma=1e-9).items()},
            point_load_array=np.empty((0, 9)), distributed_load_array=np.empty((0, 9)), job_results_dir=jr)
        el.logger_operator = None
        M = el.element_mass_matrix().M_e
        assert M.shape == (18, 18)
        assert np.allclose(M, M.T) and np.linalg.eigvalsh(M).min() > -1e-12 * np.abs(M).max()
        # a rigid translation in y carries the section mass times the length
        q = np.zeros(18)
        q[[1, 7]] = 1.0
        assert q @ M @ q == pytest.approx(sa.mass["m"] * 2.0, rel=1e-10)
        # a uniform unit mode amplitude carries int rho_A |d|^2 ds times the length
        q = np.zeros(18)
        q[[13, 16]] = 1.0  # a_1 at node 1 and node 2 (after psi)
        assert q @ M @ q == pytest.approx(sa.gbt_mass[7, 7] * 2.0, rel=1e-10)
    finally:
        shutil.rmtree(os.path.dirname(jr), ignore_errors=True)


def test_process_job_with_gbt_modes(tmp_path):
    """End to end: the [modes] column of element_section.txt sets the DOF layout and the GBT element."""
    import csv

    from workflow_orchestrator.run_job import process_job, setup_job_results_directory

    L, n_el, T = 3.0, 6, 50.0
    b, h, t = 0.2, 0.3, 2e-3
    job = tmp_path / "job_gbt"
    job.mkdir()
    nodes = "\n".join(f"{i}  {L * i / n_el:.6f}  0.0  0.0" for i in range(n_el + 1))
    (job / "grid.txt").write_text(f"[Grid]\n[node_id] [x] [y] [z]\n{nodes}\n", encoding="utf-8")
    rows = "\n".join(f"{e} {e} {e + 1} {GBT} 2 2 2 2 2 2 2 1" for e in range(n_el))
    (job / "element.txt").write_text(
        "[Element]\n[element_id] [node1] [node2] [element_type] [axial_order] [bending_y_order] "
        f"[bending_z_order] [shear_y_order] [shear_z_order] [torsion_order] [load_order] [warping]\n{rows}\n",
        encoding="utf-8",
    )
    (job / "material.txt").write_text(
        "[Material]\n[element_id] [E] [G] [nu] [rho]\n" + "\n".join(f"{e} 7e10 2.7e10 0.3 0" for e in range(n_el)) + "\n",
        encoding="utf-8",
    )
    (job / "section.txt").write_text(
        "[Section]\n[element_id] [A] [I_x] [I_y] [I_z] [J_t] [kappa] [alpha]\n"
        + "\n".join(f"{e} 1e-3 0.0 1e-8 1e-8 1e-9 0.8333 0.0" for e in range(n_el)) + "\n",
        encoding="utf-8",
    )
    dofs = ["UX", "UY", "UZ", "RX", "RY", "RZ", "W", "Q2", "Q3", "Q4", "Q5"]
    (job / "prescribed_displacement.txt").write_text(
        "[Prescribed Displacement]\n[id] [node_id] [dof] [value] [type] [comment]\n"
        + "\n".join(f"{i} 0 {d} 0.0 displacement" for i, d in enumerate(dofs)) + "\n",
        encoding="utf-8",
    )
    (job / "point_load.txt").write_text(
        "[Point load]\n[x] [y] [z] [F_x] [F_y] [F_z] [M_x] [M_y] [M_z]\n" f"{L} 0.0 0.0 0.0 0.0 0.0 {T} 0.0 0.0\n",
        encoding="utf-8",
    )
    (job / "simulation_settings.txt").write_text(
        "[Simulation]\n[Type]\nstatic\n\n[Solver]\ntype = direct\ntolerance = 1e-12\nmax_iterations = 500\n",
        encoding="utf-8",
    )
    (job / "laminate.txt").write_text(
        "[Ply_Material]\n[material_id] [E1] [E2] [G12] [nu12] [G13] [G23] [rho]\n"
        f"1 {E} {E} {G} {NU} {G} {G} 0\n"
        "[Laminate]\n[laminate_id] [ply] [material_id] [angle_deg] [thickness]\n"
        f"1 1 1 0.0 {t}\n",
        encoding="utf-8",
    )
    (job / "section_geometry.txt").write_text(
        "[Section_Point]\n[section_id] [point_id] [y] [z]\n"
        f"1 1 {-b / 2} {h / 2}\n1 2 0 {h / 2}\n1 3 {b / 2} {h / 2}\n"
        f"1 4 {-b / 2} {-h / 2}\n1 5 0 {-h / 2}\n1 6 {b / 2} {-h / 2}\n"
        "[Section_Wall]\n[section_id] [wall_id] [point_start] [point_end] [laminate_id]\n"
        "1 1 1 2 1\n1 2 2 3 1\n1 3 4 5 1\n1 4 5 6 1\n1 5 5 2 1\n",
        encoding="utf-8",
    )
    (job / "element_section.txt").write_text(
        "[Element_Section]\n[element_id] [section_id] [modes]\n" + "\n".join(f"{e} 1 2" for e in range(n_el)) + "\n",
        encoding="utf-8",
    )
    res = setup_job_results_directory("test_gbt_job")
    process_job(str(job), res, {}, {}, force_serial=True, max_processes_per_job=1)
    with open(os.path.join(res, "primary_results", "global", "U_global.csv"), newline="", encoding="utf-8") as f:
        rows = list(csv.reader(f))
    U = np.array([float(r[1]) for r in rows[1:]])
    assert U.size == 11 * (n_el + 1)
    sa = SectionAnalysis.from_thin_walled(i_beam(b, h, t), gbt_modes=2)
    expected = solve_cantilever(GBT, sa, L, n_el, 2, {3: T})
    np.testing.assert_allclose(U.reshape(n_el + 1, 11), expected, rtol=1e-9, atol=1e-12)
    assert np.abs(expected[:, 7:]).max() > 0

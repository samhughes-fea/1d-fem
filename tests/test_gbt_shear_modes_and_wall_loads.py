"""Shear-lag warping modes, wall surface loads and their job inputs."""

from __future__ import annotations

import os

import numpy as np
import pytest

from pre_processing.section_integration.section_analysis import SectionAnalysis
from pre_processing.section_integration.thin_walled.cross_section_modes import CrossSectionModes
from test_gbt_cross_section_modes import E, G, GBT, NU, box, i_beam


def test_warping_modes_span_saint_venant_and_are_orthogonal():
    an = box(0.6, 0.15, n_sub=4)
    cm = CrossSectionModes(an)
    n_sv, p = cm.n_sv_warping, cm.n_shear_modes
    assert n_sv == 3 and p > 3
    Om = np.column_stack([cm.warping_mode(j) for j in range(n_sv + p)])
    A = cm._axial_warping_stiffness
    Em = an._matrices["E"]
    # orthogonal to the rigid axial fields 1, y, z and to each other in A
    assert np.abs(an._rigid_modes.T @ A @ Om).max() < 1e-9 * np.abs(A).max() * np.abs(Om).max()
    An = Om.T @ A @ Om
    np.testing.assert_allclose(An / np.sqrt(np.outer(np.diag(An), np.diag(An))), np.eye(n_sv + p), atol=1e-8)
    np.testing.assert_allclose(np.einsum("ik,ij,jk->k", Om, Em, Om) / np.diag(An), cm.warping_eigenvalues, rtol=1e-10)
    # the first n_sv modes span the Saint-Venant warping (less its rigid part)
    W = np.hstack([an._solution["W0"], an._solution["W1"]])
    Q = np.hstack([an._rigid_modes, Om[:, :n_sv]])
    resid = W - Q @ np.linalg.lstsq(Q, W, rcond=None)[0]
    assert np.abs(resid).max() < 1e-8 * np.abs(W).max()


def test_wall_load_resultant_hand_values():
    b, h, t = 0.2, 0.3, 2e-3
    an = i_beam(b, h, t)
    p = 1.5e3
    # wall 2: top flange right half, from (0, h/2) to (b/2, h/2): tangent +y, normal +z
    f = an.wall_load_resultant({2: (0.0, 0.0, p)})
    np.testing.assert_allclose(f[:6], [0, 0, p * b / 2, p * b**2 / 8, 0, 0], atol=1e-9)
    # axial traction on the web (wall 5 from (0, -h/2) to (0, h/2)): N and M_y = int z p_x = 0
    f = an.wall_load_resultant({5: (p, 0.0, 0.0)})
    np.testing.assert_allclose(f[:6], [p * h, 0, 0, 0, 0, 0], atol=1e-9)
    # tangential load on the web: F_z, no torque about the web line
    f = an.wall_load_resultant({5: (0.0, p, 0.0)})
    np.testing.assert_allclose(f[:6], [0, 0, p * h, 0, 0, 0], atol=1e-9)
    # the GBT load starts with the same six beam entries
    cm = CrossSectionModes(an)
    g = cm.wall_load({2: (0.0, 0.0, p), 5: (p, p, 0.0)}, 2, 0)
    np.testing.assert_allclose(g[:6], an.wall_load_resultant({2: (0.0, 0.0, p), 5: (p, p, 0.0)})[:6], atol=1e-12)
    assert g.size == 6 + cm.n_sv_warping + 4
    assert np.abs(g[6 + cm.n_sv_warping:]).max() > 0  # an eccentric flange load drives distortion
    with pytest.raises(ValueError, match="unknown walls"):
        an.wall_load_resultant({9: (1, 0, 0)})


def test_section_mesh_and_wall_load_parsers(tmp_path):
    from pre_processing.parsing.section_geometry_parser import parse_section_mesh
    from pre_processing.parsing.wall_load_parser import parse_wall_load_file

    g = tmp_path / "section_geometry.txt"
    g.write_text(
        "[Section_Point]\n[section_id] [point_id] [y] [z]\n1 1 0 0\n1 2 1 0\n1 3 1 1\n"
        "[Section_Wall]\n[section_id] [wall_id] [point_start] [point_end] [laminate_id]\n1 1 1 2 1\n1 2 2 3 1\n"
        "[Section_Mesh]\n[section_id] [n_sub]\n1 3\n",
        encoding="utf-8",
    )
    assert parse_section_mesh(str(g)) == {1: 3}
    w = tmp_path / "wall_load.txt"
    w.write_text("[Wall_Load]\n[element_id] [wall_id] [p_x] [p_s] [p_n]\n0 1 0 0 2\n0 1 1 0 3\n2 2 0 1 0\n",
                 encoding="utf-8")
    assert parse_wall_load_file(str(w)) == {0: {1: (1.0, 0.0, 5.0)}, 2: {2: (0.0, 1.0, 0.0)}}


def test_process_job_with_shear_modes_and_wall_loads(tmp_path):
    """End to end: [shear_modes], [Section_Mesh] and wall_load.txt through process_job."""
    import csv

    from pre_processing.section_integration.thin_walled.job_input import section_analyses_from_files
    from workflow_orchestrator.run_job import process_job, setup_job_results_directory

    L, n_el, p = 2.0, 8, 2.0e3
    b, h, t = 0.6, 0.15, 2e-3
    job = tmp_path / "job_gbt_loads"
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
    n_sv = CrossSectionModes(box(b, h, t, n_sub=2)).n_sv_warping
    dpn = 6 + n_sv + 1 + 2
    dofs = ["UX", "UY", "UZ", "RX", "RY", "RZ"] + [f"Q{k}" for k in range(1, dpn - 5)]
    (job / "prescribed_displacement.txt").write_text(
        "[Prescribed Displacement]\n[id] [node_id] [dof] [value] [type] [comment]\n"
        + "\n".join(f"{i} 0 {d} 0.0 displacement" for i, d in enumerate(dofs)) + "\n",
        encoding="utf-8",
    )
    (job / "point_load.txt").write_text("[Point load]\n[x] [y] [z] [F_x] [F_y] [F_z] [M_x] [M_y] [M_z]\n", encoding="utf-8")
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
        f"1 1 {-b / 2} {-h / 2}\n1 2 {b / 2} {-h / 2}\n1 3 {b / 2} {h / 2}\n1 4 {-b / 2} {h / 2}\n"
        "[Section_Wall]\n[section_id] [wall_id] [point_start] [point_end] [laminate_id]\n"
        "1 1 1 2 1\n1 2 2 3 1\n1 3 3 4 1\n1 4 4 1 1\n"
        "[Section_Mesh]\n[section_id] [n_sub]\n1 2\n",
        encoding="utf-8",
    )
    (job / "element_section.txt").write_text(
        "[Element_Section]\n[element_id] [section_id] [modes] [shear_modes]\n"
        + "\n".join(f"{e} 1 1 1" for e in range(n_el)) + "\n",
        encoding="utf-8",
    )
    # downward pressure on the top flange (wall 3, normal -z) over the outer half of the span
    (job / "wall_load.txt").write_text(
        "[Wall_Load]\n[element_id] [wall_id] [p_x] [p_s] [p_n]\n"
        + "\n".join(f"{e} 3 0.0 0.0 {p}" for e in range(n_el // 2, n_el)) + "\n",
        encoding="utf-8",
    )
    res = setup_job_results_directory("test_gbt_wall_load_job")
    process_job(str(job), res, {}, {}, force_serial=True, max_processes_per_job=1)
    with open(os.path.join(res, "primary_results", "global", "U_global.csv"), newline="", encoding="utf-8") as f:
        rows = list(csv.reader(f))
    U = np.array([float(r[1]) for r in rows[1:]]).reshape(n_el + 1, dpn)

    # reference: direct assembly with the same generalised loads
    an = section_analyses_from_files(str(job / "laminate.txt"), str(job / "section_geometry.txt"))[1]
    assert an.n_sub == 2
    sa = SectionAnalysis.from_thin_walled(an, gbt_modes=1, gbt_shear_modes=1)
    f = CrossSectionModes(an).wall_load({3: (0.0, 0.0, p)}, 1, 1)
    # total vertical load: p b over half the span, downward (wall 3 runs -y, normal -z)
    assert f[2] == pytest.approx(-p * b, rel=1e-12)
    from pre_processing.element_library.element_factory import ElementFactory  # noqa: F401
    import importlib
    import shutil
    import tempfile

    from pre_processing.element_library.dof_layout import mesh_dof_per_node

    Cls = getattr(importlib.import_module(ElementFactory.ELEMENT_CLASS_MAP[GBT]), GBT)
    import test_gbt_cross_section_modes as T

    ed = T._ed(GBT, n_el, 1, sa, 1)
    ed["section_analysis_per_element"] = tuple([sa] * n_el)
    ed["wall_load_per_element"] = tuple(None if e < n_el // 2 else f for e in range(n_el))
    x = np.linspace(0.0, L, n_el + 1)
    jr = os.path.join(tempfile.mkdtemp(), "job")
    for sub in ("element_stiffness_matrices", "element_force_vectors", "logs"):
        os.makedirs(os.path.join(jr, sub), exist_ok=True)
    try:
        dpn = mesh_dof_per_node(ed)
        Kg = np.zeros((dpn * (n_el + 1),) * 2)
        Fg = np.zeros(dpn * (n_el + 1))
        for e in range(n_el):
            el = Cls(element_id=e, element_dictionary=ed, grid_dictionary={"coordinates": np.column_stack([x, 0 * x, 0 * x])},
                     material_dictionary={k: np.full(n_el, v) for k, v in dict(E=E, G=G, nu=NU, rho=0.0).items()},
                     section_dictionary={k: np.full(n_el, v) for k, v in dict(A=1e-3, I_x=0.0, I_y=1e-8, I_z=1e-8, J_t=1e-9,
                                                                               kappa=5 / 6, alpha=0.0).items()},
                     point_load_array=np.empty((0, 9)), distributed_load_array=np.empty((0, 9)), job_results_dir=jr)
            el.logger_operator = None
            gi = el.assemble_global_dof_indices()
            Kg[np.ix_(gi, gi)] += el.element_stiffness_matrix().K_e
            Fg[gi] += el.element_force_vector().F_e
        free = np.arange(dpn, dpn * (n_el + 1))
        Uref = np.zeros_like(Fg)
        Uref[free] = np.linalg.solve(Kg[np.ix_(free, free)], Fg[free])
    finally:
        shutil.rmtree(os.path.dirname(jr), ignore_errors=True)
    np.testing.assert_allclose(U, Uref.reshape(n_el + 1, dpn), rtol=1e-9, atol=1e-13)
    assert U[-1, 2] < 0


def test_wall_load_on_timoshenko_and_benscoter_elements_is_lumped_by_halves():
    import shutil
    import tempfile

    import test_gbt_cross_section_modes as T
    from pre_processing.element_library.linear.beam.first_order_shear_deformation_theory.benscoter.linear_benscoter_3D import (
        LinearBenscoterBeamElement3D,
    )

    an = i_beam()
    sa = SectionAnalysis.from_thin_walled(an)
    f = an.wall_load_resultant({2: (0.0, 0.0, 1.5e3), 5: (2.0e3, 0.0, 0.0)})
    ed = T._ed(T.BS, 1, 0)
    ed["section_analysis_per_element"] = (sa,)
    ed["wall_load_per_element"] = (f,)
    jr = os.path.join(tempfile.mkdtemp(), "job")
    for sub in ("element_stiffness_matrices", "element_force_vectors", "logs"):
        os.makedirs(os.path.join(jr, sub), exist_ok=True)
    try:
        el = LinearBenscoterBeamElement3D(
            element_id=0, element_dictionary=ed, grid_dictionary={"coordinates": np.array([[0.0, 0, 0], [2.0, 0, 0]])},
            material_dictionary={k: np.array([v]) for k, v in dict(E=E, G=G, nu=NU, rho=0.0).items()},
            section_dictionary={k: np.array([v]) for k, v in dict(A=1e-3, I_x=0.0, I_y=1e-8, I_z=1e-8, J_t=1e-9,
                                                                   kappa=5 / 6, alpha=0.0).items()},
            point_load_array=np.empty((0, 9)), distributed_load_array=np.empty((0, 9)), job_results_dir=jr)
        el.logger_operator = None
        Fe = el.element_force_vector().F_e
    finally:
        shutil.rmtree(os.path.dirname(jr), ignore_errors=True)
    np.testing.assert_allclose(Fe[:6], f[:6], rtol=1e-12)    # half of 2 m of load
    np.testing.assert_allclose(Fe[6:12], f[:6], rtol=1e-12)
    np.testing.assert_allclose(Fe[[12, 13]], [f[6], f[6]], rtol=1e-12, atol=1e-30)

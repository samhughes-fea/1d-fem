"""Benscoter restrained-warping beam element with thin-walled sections.

The restrained-torsion cantilever is checked against the exact solution of the Benscoter
equations with the section's own stiffnesses: for a tip torque T, root clamped with
``psi = 0`` and a free tip,

    k^2 = S GJ / ((GJ + S) E Gamma),  psi = T/GJ (1 - cosh kx + tanh kL sinh kx),
    theta_x' = (T + S psi) / (GJ + S).
"""

from __future__ import annotations

import os
import shutil
import tempfile

import numpy as np
import pytest

from pre_processing.section_integration.section_analysis import SectionAnalysis
from pre_processing.section_integration.thin_walled import (
    Ply,
    PlyMaterial,
    ThinWalledSection,
    ThinWalledSectionAnalysis,
    Wall,
    WallLaminate,
)

E, NU = 70e9, 0.3
G = E / (2 * (1 + NU))
CFRP = PlyMaterial(E1=140e9, E2=10e9, G12=5e9, nu12=0.3, G13=5e9, G23=3.5e9, rho=1600.0)
BS = "LinearBenscoterBeamElement3D"
TS = "LinearTimoshenkoBeamElement3D"


def iso(t, rho=0.0):
    return WallLaminate.isotropic(E, NU, t, rho=rho)


def i_beam(b=0.2, h=0.3, tf=2e-3, tw=2e-3, rho=0.0):
    pts = {1: (-b / 2, h / 2), 2: (0, h / 2), 3: (b / 2, h / 2),
           4: (-b / 2, -h / 2), 5: (0, -h / 2), 6: (b / 2, -h / 2)}
    walls = [Wall(1, 1, 2, 1), Wall(2, 2, 3, 1), Wall(3, 4, 5, 1), Wall(4, 5, 6, 1), Wall(5, 5, 2, 2)]
    return ThinWalledSectionAnalysis(ThinWalledSection(pts, walls), {1: iso(tf, rho), 2: iso(tw, rho)}, hoop="free")


def box(b=0.3, h=0.15, t=2e-3):
    pts = {1: (-b / 2, -h / 2), 2: (b / 2, -h / 2), 3: (b / 2, h / 2), 4: (-b / 2, h / 2)}
    walls = [Wall(1, 1, 2, 1), Wall(2, 2, 3, 1), Wall(3, 3, 4, 1), Wall(4, 4, 1, 1)]
    return ThinWalledSectionAnalysis(ThinWalledSection(pts, walls), {1: iso(t)}, hoop="free")


def cas_box():
    top = WallLaminate([Ply(CFRP, 20.0, 0.5e-3)] * 4)
    bottom = WallLaminate([Ply(CFRP, -20.0, 0.5e-3)] * 4)
    webs = WallLaminate([Ply(CFRP, a, 0.5e-3) for a in (45, -45, -45, 45)])
    b, h = 0.3, 0.15
    pts = {1: (-b / 2, -h / 2), 2: (b / 2, -h / 2), 3: (b / 2, h / 2), 4: (-b / 2, h / 2)}
    walls = [Wall(1, 1, 2, 2), Wall(2, 2, 3, 3), Wall(3, 3, 4, 1), Wall(4, 4, 1, 3)]
    return ThinWalledSectionAnalysis(ThinWalledSection(pts, walls), {1: top, 2: bottom, 3: webs})


def _element_class(etype):
    if etype == BS:
        from pre_processing.element_library.linear.beam.first_order_shear_deformation_theory.benscoter.linear_benscoter_3D import (
            LinearBenscoterBeamElement3D as Cls,
        )
    else:
        from pre_processing.element_library.linear.beam.first_order_shear_deformation_theory.timoshenko.linear_timoshenko_3D import (
            LinearTimoshenkoBeamElement3D as Cls,
        )
    return Cls


def _build(etype, section, L, n_el, warping):
    jr = os.path.join(tempfile.mkdtemp(), "job")
    for sub in ("element_stiffness_matrices", "element_force_vectors", "logs"):
        os.makedirs(os.path.join(jr, sub), exist_ok=True)
    x = np.linspace(0.0, L, n_el + 1)
    coords = np.column_stack([x, np.zeros_like(x), np.zeros_like(x)])
    o = {"axial": 2, "bending_y": 2, "bending_z": 2, "shear_y": 2, "shear_z": 2, "torsion": 2, "load": 2}
    ed = {
        "ids": np.arange(n_el),
        "connectivity": np.column_stack([np.arange(n_el), np.arange(1, n_el + 1)]),
        "types": np.array([etype] * n_el),
        "warping": np.full(n_el, int(warping), dtype=np.int8),
        "integration_orders": {k: np.full(n_el, v) for k, v in o.items()},
    }
    if section is not None:
        ed["section_analysis_per_element"] = tuple([section] * n_el)
    mat = {k: np.full(n_el, v) for k, v in dict(E=E, G=G, nu=NU, rho=0.0).items()}
    sec = {k: np.full(n_el, v) for k, v in dict(A=1e-3, I_x=0.0, I_y=1e-8, I_z=1e-8, J_t=1e-9, kappa=5 / 6,
                                                alpha=0.0, y_sc=0.0, z_sc=0.0, Gamma=1e-9).items()}
    Cls = _element_class(etype)
    els = []
    for e in range(n_el):
        el = Cls(
            element_id=e, element_dictionary=ed, grid_dictionary={"coordinates": coords},
            material_dictionary=mat, section_dictionary=sec,
            point_load_array=np.empty((0, 9)), distributed_load_array=np.empty((0, 9)),
            job_results_dir=jr,
        )
        el.logger_operator = None
        els.append(el)
    return els, jr


def fe_cantilever(etype, section, L, P, n_el=16, restrain_root_warping=True):
    """Tip displacements ``[u, v, w, theta_x, theta_y, theta_z]`` and nodal ``psi`` (Benscoter)."""
    warping = etype == BS
    els, jr = _build(etype, section, L, n_el, warping)
    try:
        dpn = 7 if warping else 6
        n_dof = dpn * (n_el + 1)
        K = np.zeros((n_dof, n_dof))
        for e, el in enumerate(els):
            Ke = el.element_stiffness_matrix().K_e
            # element layout: 12 standard (node-major) then psi1, psi2; global: 7 per node
            gmap = [dpn * e + i for i in range(6)] + [dpn * (e + 1) + i for i in range(6)]
            if warping:
                gmap += [dpn * e + 6, dpn * (e + 1) + 6]
            K[np.ix_(gmap, gmap)] += Ke
        F = np.zeros(n_dof)
        F[dpn * n_el:dpn * n_el + 6] = P
        fixed = list(range(6)) + ([6] if warping and restrain_root_warping else [])
        free = np.setdiff1d(np.arange(n_dof), fixed)
        U = np.zeros(n_dof)
        U[free] = np.linalg.solve(K[np.ix_(free, free)], F[free])
        tip = U[dpn * n_el:dpn * n_el + 6]
        psi = U[6::dpn] if warping else None
        return tip, psi, els
    finally:
        shutil.rmtree(os.path.dirname(jr), ignore_errors=True)


def exact_restrained_torsion(GJ, S, EG, T, L, x):
    k = np.sqrt(S * GJ / ((GJ + S) * EG))
    psi = T / GJ * (1 - np.cosh(k * x) + np.tanh(k * L) * np.sinh(k * x))
    # theta(x) = int_0^x (T + S psi) / (GJ + S)
    ipsi = T / GJ * (x - np.sinh(k * x) / k + np.tanh(k * L) * (np.cosh(k * x) - 1) / k)
    theta = (T * x + S * ipsi) / (GJ + S)
    return theta, psi


# --------------------------------------------------------------------------- #
# Section stiffness
# --------------------------------------------------------------------------- #


def test_benscoter_section_stiffness_of_i_beam():
    b, h, tf = 0.2, 0.3, 2e-3
    sa = i_beam(b, h, tf)
    K8 = sa.benscoter_stiffness
    np.testing.assert_allclose(K8[:6, :6], sa.stiffness, rtol=1e-12, atol=1e-12 * np.abs(sa.stiffness).max())
    assert K8[6, 6] == pytest.approx(G * 2 * b * tf * (h / 2) ** 2, rel=1e-10)
    assert K8[7, 7] == pytest.approx(E * tf * b**3 * h**2 / 24, rel=1e-10)
    off = K8.copy()
    np.fill_diagonal(off, 0.0)
    assert np.abs(off[6:, :]).max() < 1e-9 * np.abs(K8).max()
    assert np.allclose(K8, K8.T)
    assert np.all(np.linalg.eigvalsh(K8) > 0)


def test_benscoter_stiffness_element_signs_and_warping_inertia():
    tw = cas_box()
    sa = SectionAnalysis.from_thin_walled(tw)
    S8 = np.diag([1.0, -1, 1, 1, 1, 1, 1, 1])
    np.testing.assert_allclose(sa.benscoter_stiffness, S8 @ tw.benscoter_stiffness @ S8)
    np.testing.assert_allclose(sa.benscoter_stiffness[:6, :6], sa.stiffness, rtol=1e-12,
                               atol=1e-12 * np.abs(sa.stiffness).max())
    # warping inertia of the isotropic I-beam: rho t int omega^2 ds = rho * Gamma_primary
    b, h, tf, rho = 0.2, 0.3, 2e-3, 2700.0
    assert i_beam(b, h, tf, rho=rho).warping_inertia == pytest.approx(rho * tf * b**3 * h**2 / 24, rel=1e-10)


# --------------------------------------------------------------------------- #
# Element
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize("make", [i_beam, box])
def test_restrained_torsion_cantilever_matches_benscoter_ode(make):
    tw = make()
    sa = SectionAnalysis.from_thin_walled(tw)
    K8 = sa.benscoter_stiffness
    GJ, S, EG = K8[5, 5], K8[6, 6], K8[7, 7]
    T = 100.0
    k = np.sqrt(S * GJ / ((GJ + S) * EG))
    L = 6.0 / k if make is i_beam else 2.0  # the box layer is very short, use a fixed length
    theta_ex, psi_ex = exact_restrained_torsion(GJ, S, EG, T, L, np.linspace(0, L, 129))
    errs = []
    for n_el in (32, 128):
        tip, psi, _ = fe_cantilever(BS, sa, L, [0, 0, 0, T, 0, 0], n_el=n_el)
        errs.append(abs(tip[3] - theta_ex[-1]) / abs(theta_ex[-1]))
        if n_el == 128:
            np.testing.assert_allclose(psi, psi_ex, atol=2e-3 * np.abs(psi_ex).max())
    assert errs[1] < 1e-4
    assert errs[1] < errs[0] / 8 or errs[1] < 1e-9
    # restraint stiffens the open section noticeably, not the closed one
    free = T * L / GJ
    if make is i_beam:
        assert theta_ex[-1] < 0.9 * free
    else:
        assert theta_ex[-1] > 0.98 * free


@pytest.mark.parametrize("load", [[1e4, 0, 0, 0, 0, 0], [0, 0, 0, 2e2, 0, 0], [0, 0, 0, 0, 3e2, 0], [0, 0, 0, 0, 0, 3e2]])
def test_free_warping_reproduces_timoshenko_for_constant_section_forces(load):
    sa = SectionAnalysis.from_thin_walled(cas_box())
    L = 2.0
    tip_b, psi, _ = fe_cantilever(BS, sa, L, load, n_el=4, restrain_root_warping=False)
    tip_t, _, _ = fe_cantilever(TS, sa, L, load, n_el=4)
    np.testing.assert_allclose(tip_b, tip_t, rtol=1e-8, atol=1e-10 * np.abs(tip_t).max())


def test_free_warping_follows_twist_rate():
    sa = SectionAnalysis.from_thin_walled(i_beam())
    T, L = 50.0, 3.0
    tip, psi, _ = fe_cantilever(BS, sa, L, [0, 0, 0, T, 0, 0], n_el=6, restrain_root_warping=False)
    GJ = sa.benscoter_stiffness[5, 5]
    assert tip[3] == pytest.approx(T * L / GJ, rel=1e-10)
    np.testing.assert_allclose(psi, T / GJ, rtol=1e-10)


def test_benscoter_mass_matrix_has_warping_inertia():
    tw = i_beam(rho=2700.0)
    sa = SectionAnalysis.from_thin_walled(tw)
    els, jr = _build(BS, sa, 1.0, 2, True)
    try:
        M = els[0].element_mass_matrix().M_e
        Le = 0.5
        assert M.shape == (14, 14)
        assert M[12, 12] == pytest.approx(sa.warping_inertia * Le / 3, rel=1e-12)
        assert M[12, 13] == pytest.approx(sa.warping_inertia * Le / 6, rel=1e-12)
        assert np.allclose(M, M.T)
        # translational mass is the section mass per length
        assert M[0, 0] + M[0, 6] == pytest.approx(sa.mass["m"] * Le / 2, rel=1e-10)
    finally:
        shutil.rmtree(os.path.dirname(jr), ignore_errors=True)


def test_benscoter_needs_thin_walled_section_and_warping():
    sa = SectionAnalysis.from_thin_walled(i_beam())
    with pytest.raises(ValueError, match="warping"):
        _build(BS, sa, 1.0, 1, False)
    with pytest.raises(ValueError, match="thin-walled"):
        _build(BS, None, 1.0, 1, True)


def test_benscoter_is_registered():
    from pre_processing.element_library.element_factory import ElementFactory
    from pre_processing.element_library.shape_function_registry import SHAPE_FUNCTION_REGISTRY

    assert BS in ElementFactory.ELEMENT_CLASS_MAP
    assert BS in SHAPE_FUNCTION_REGISTRY


def test_process_job_with_benscoter_elements(tmp_path):
    """End to end: element_section.txt plus [warping] = 1 drive a restrained-torsion cantilever."""
    import csv

    from workflow_orchestrator.run_job import process_job, setup_job_results_directory

    L, n_el, T = 20.0, 8, 100.0
    b, h, t = 0.2, 0.3, 2e-3
    job = tmp_path / "job_benscoter"
    job.mkdir()
    nodes = "\n".join(f"{i}  {L * i / n_el:.6f}  0.0  0.0" for i in range(n_el + 1))
    (job / "grid.txt").write_text(f"[Grid]\n[node_id] [x] [y] [z]\n{nodes}\n", encoding="utf-8")
    rows = "\n".join(f"{e} {e} {e + 1} {BS} 2 2 2 2 2 2 2 1" for e in range(n_el))
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
    dofs = ["UX", "UY", "UZ", "RX", "RY", "RZ", "W"]
    (job / "prescribed_displacement.txt").write_text(
        "[Prescribed Displacement]\n[id] [node_id] [dof] [value] [type] [comment]\n"
        + "\n".join(f"{i} 0 {d} 0.0 displacement" for i, d in enumerate(dofs)) + "\n",
        encoding="utf-8",
    )
    (job / "point_load.txt").write_text(
        "[Point load]\n[x] [y] [z] [F_x] [F_y] [F_z] [M_x] [M_y] [M_z]\n"
        f"{L} 0.0 0.0 0.0 0.0 0.0 {T} 0.0 0.0\n",
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
        "[Element_Section]\n[element_id] [section_id]\n" + "\n".join(f"{e} 1" for e in range(n_el)) + "\n",
        encoding="utf-8",
    )
    res = setup_job_results_directory("test_benscoter_job")
    process_job(str(job), res, {}, {}, force_serial=True, max_processes_per_job=1)
    with open(os.path.join(res, "primary_results", "global", "U_global.csv"), newline="", encoding="utf-8") as f:
        rows = list(csv.reader(f))
    U = np.array([float(r[1]) for r in rows[1:]]) if len(rows) > 2 else np.array([float(x) for x in rows[1][1:]])
    assert U.size == 7 * (n_el + 1)
    from pre_processing.section_integration.thin_walled.job_input import section_analyses_from_files

    tw = section_analyses_from_files(str(job / "laminate.txt"), str(job / "section_geometry.txt"))[1]
    expected_tip, expected_psi, _ = fe_cantilever(BS, SectionAnalysis.from_thin_walled(tw), L, [0, 0, 0, T, 0, 0], n_el=n_el)
    np.testing.assert_allclose(U[7 * n_el:7 * n_el + 6], expected_tip, rtol=1e-9, atol=1e-12)
    np.testing.assert_allclose(U[6::7], expected_psi, rtol=1e-9, atol=1e-12)

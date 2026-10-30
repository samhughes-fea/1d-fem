"""2D section solver (Giavotto / BECAS): exact cases, and thin-walled and 2D results on the same sections."""

from __future__ import annotations

import os

import numpy as np
import pytest

from pre_processing.section_integration.section_analysis import ELEMENT_SIGN, SectionAnalysis
from pre_processing.section_integration.solid.material import SolidMaterial
from pre_processing.section_integration.solid.mesh import SectionMesh
from pre_processing.section_integration.solid.section_analysis import SolidSectionAnalysis
from pre_processing.section_integration.thin_walled import ThinWalledSection, Wall, WallLaminate
from pre_processing.section_integration.thin_walled.job_input import attach_section_analyses
from pre_processing.section_integration.thin_walled.laminate import Ply, PlyMaterial
from pre_processing.section_integration.thin_walled.section_analysis import (
    ThinWalledSectionAnalysis,
    stiffness_about,
)

CFRP = PlyMaterial(E1=135e9, E2=10e9, G12=5e9, nu12=0.3, G13=5e9, G23=3.8e9, rho=1550.0)
ORTHO = SolidMaterial(135e9, 10e9, 10e9, 5e9, 5e9, 3.8e9, 0.3, 0.3, 0.32, 1550.0)


def _scaled(Ka, Kb):
    """``(Ka - Kb) / sqrt(Kb_ii Kb_jj)``: differences relative to the diagonal."""
    d = np.sqrt(np.outer(np.diag(Kb), np.diag(Kb)))
    return (Ka - Kb) / d


# ------------------------------------------------------------------ exact cases
def test_isotropic_rectangle():
    E, nu = 70e9, 0.3
    G = E / (2 * (1 + nu))
    b, h = 0.05, 0.1
    a = SolidSectionAnalysis(SectionMesh.rectangle(b, h, 8, 16, SolidMaterial.isotropic(E, nu, rho=2700.0)))
    K = a.stiffness
    A = b * h
    assert K[0, 0] == pytest.approx(E * A, rel=1e-12)
    assert K[1, 1] == pytest.approx(E * b * h**3 / 12, rel=1e-12)
    assert K[2, 2] == pytest.approx(E * h * b**3 / 12, rel=1e-12)
    # Saint-Venant torsion of a 2:1 rectangle: J = 0.228682 h b^3 (series solution)
    assert K[5, 5] == pytest.approx(0.228682 * G * h * b**3, rel=1e-4)
    np.testing.assert_allclose(K - np.diag(np.diag(K)), 0.0, atol=1e-10 * K[0, 0])
    np.testing.assert_allclose(a.shear_centre, 0.0, atol=1e-12)
    np.testing.assert_allclose(a.tension_centre, 0.0, atol=1e-12)
    assert a.mass_properties["m"] == pytest.approx(2700.0 * A)
    assert a.mass_properties["I_p"] == pytest.approx(2700.0 * A * (b * b + h * h) / 12)
    # shear-rigid stiffness: the bending and torsion block of K (uncoupled section)
    np.testing.assert_allclose(a.classical_stiffness, K[np.ix_([0, 1, 2, 5], [0, 1, 2, 5])],
                               rtol=1e-12, atol=1e-12 * K[5, 5])


def test_shear_coefficient_is_five_sixths_without_poisson():
    """With nu = 0 the flexure shear stress is the parabola of elementary theory, k = 5/6."""
    a = SolidSectionAnalysis(SectionMesh.rectangle(1.0, 1.0, 8, 8, SolidMaterial.isotropic(1.0, 0.0)))
    G = 0.5
    assert a.stiffness[3, 3] / G == pytest.approx(5 / 6, rel=1e-4)
    assert a.stiffness[4, 4] / G == pytest.approx(5 / 6, rel=1e-4)


def test_homogeneous_anisotropic_rectangle_is_exact():
    """A homogeneous section under N or M has a stress of elementary theory (uniform or
    linear sigma_xx alone); its strains then give the flexibility exactly on any mesh.
    A shear strain that varies linearly across the section is a twist of half its slope."""
    b, h, c = 0.03, 0.05, (0.01, -0.02)
    mesh = SectionMesh.rectangle(b, h, 3, 4, ORTHO, centre=c, theta=30.0, alpha=20.0)
    a = SolidSectionAnalysis(mesh)
    S = np.linalg.inv(mesh.element_stiffness(0))         # [xx, yy, zz, yz, xz, xy]
    A, I_y, I_z = b * h, b * h**3 / 12, h * b**3 / 12
    F = np.linalg.inv(stiffness_about(a.stiffness, *c))  # about the centre
    sxx, sxz, sxy = S[0, 0], S[4, 0], S[5, 0]
    np.testing.assert_allclose(F[:, 0], np.array([sxx, 0, 0, sxy, sxz, 0]) / A, rtol=1e-9, atol=1e-9 * sxx / A)
    # M_y: sigma = M z / I_y, so gamma_xy = sxy M z / I_y and phi = -sxy M / (2 I_y)
    np.testing.assert_allclose(F[:, 1], np.array([0, sxx, 0, 0, 0, -0.5 * sxy]) / I_y, atol=1e-9 * sxx / I_y)
    # M_z: sigma = -M y / I_z, gamma_xz = -sxz M y / I_z and phi = -sxz M / (2 I_z)
    np.testing.assert_allclose(F[:, 2], np.array([0, 0, sxx, 0, 0, -0.5 * sxz]) / I_z, atol=1e-9 * sxx / I_z)
    np.testing.assert_allclose(a.tension_centre, c, atol=1e-12)
    np.testing.assert_allclose(a.shear_centre, c, atol=1e-12)


def test_offset_isotropic_rectangle_moves_rigidly():
    mat = SolidMaterial.isotropic(1.0, 0.25)
    K0 = SolidSectionAnalysis(SectionMesh.rectangle(0.2, 0.1, 4, 2, mat)).stiffness
    a = SolidSectionAnalysis(SectionMesh.rectangle(0.2, 0.1, 4, 2, mat, centre=(0.3, -0.1)))
    np.testing.assert_allclose(stiffness_about(a.stiffness, 0.3, -0.1), K0, atol=1e-10 * K0[0, 0])
    np.testing.assert_allclose(a.shear_centre, (0.3, -0.1), atol=1e-12)
    np.testing.assert_allclose(a.tension_centre, (0.3, -0.1), atol=1e-12)


def test_point_strains_and_warping():
    a = SolidSectionAnalysis(SectionMesh.rectangle(0.04, 0.04, 4, 4, ORTHO, theta=45.0))
    F = np.array([10.0, 1.0, -2.0, 3.0, 4.0, 0.5])
    pts = a.point_strains(F)
    assert len(pts) == 16 * 9
    # the stresses carry the section forces
    N = sum(p.weight * (p.C @ e)[0] for p, (_, _, _, e) in zip(a._points, pts))
    V_z = sum(p.weight * (p.C @ e)[4] for p, (_, _, _, e) in zip(a._points, pts))
    assert N == pytest.approx(F[0], rel=1e-9)
    assert V_z == pytest.approx(F[4], rel=1e-9)
    assert a.nodal_warping(F).shape == (a.mesh.nodes.shape[0], 3)
    assert a.warping_function().shape == (a.mesh.nodes.shape[0],)


# ------------------------------------------------------------------ thin-walled and 2D
def _i_section():
    pts = {1: (-0.075, 0.15), 2: (0.0, 0.15), 3: (0.075, 0.15), 4: (-0.075, -0.15), 5: (0.0, -0.15), 6: (0.075, -0.15)}
    walls = [Wall(1, 1, 2, 1), Wall(2, 2, 3, 1), Wall(3, 4, 5, 1), Wall(4, 5, 6, 1), Wall(5, 5, 2, 2)]
    lams = {1: WallLaminate([Ply(CFRP, a, 0.5e-3) for a in (0, 45, -45, 0)]),
            2: WallLaminate([Ply(CFRP, a, 0.5e-3) for a in (45, -45, -45, 45)])}
    return ThinWalledSection(pts, walls), lams


def _cas_box():
    b = 0.1
    pts = {1: (-b / 2, -b / 2), 2: (b / 2, -b / 2), 3: (b / 2, b / 2), 4: (-b / 2, b / 2)}
    walls = [Wall(1, 1, 2, 1), Wall(2, 2, 3, 1), Wall(3, 3, 4, 1), Wall(4, 4, 1, 1)]
    return ThinWalledSection(pts, walls), {1: WallLaminate([Ply(CFRP, 30.0, 0.25e-3)] * 4)}


@pytest.fixture(scope="module")
def i_section():
    sec, lams = _i_section()
    # the 2D section bends its walls freely across their width: the thin-walled hoop="free"
    tw = ThinWalledSectionAnalysis(sec, lams, hoop="free")
    return tw, SolidSectionAnalysis(SectionMesh.from_thin_walled(sec, lams, elements_along=8))


@pytest.fixture(scope="module")
def cas_box():
    sec, lams = _cas_box()
    tw = ThinWalledSectionAnalysis(sec, lams, hoop="free")
    return tw, SolidSectionAnalysis(SectionMesh.from_thin_walled(sec, lams, elements_along=8))


def test_i_section_matches_thin_walled(i_section):
    tw, so = i_section
    Kt, Ks = tw.stiffness, so.stiffness
    D = _scaled(Ks, Kt)
    # extension, bending and torsion (the laminated walls' D66, with D26 relaxed) within 0.2 %
    for i in (0, 1, 2, 5):
        assert abs(D[i, i]) < 2e-3
    # shear along the flanges: the 2D web also shears through its thickness, which
    # Kirchhoff walls cannot (stiffer by a few per cent)
    assert 0.0 < D[3, 3] < 0.08
    assert abs(D[4, 4]) < 0.015
    off = D - np.diag(np.diag(D))
    # V_z-T is the open section's shear centre offset times sqrt(GA / GJ), a large factor:
    # compare the centres instead
    off[4, 5] = off[5, 4] = 0.0
    assert np.max(np.abs(off)) < 1e-3
    assert so.warping_constant == pytest.approx(tw.warping_constant, rel=0.02)
    assert so.mass_properties["m"] == pytest.approx(tw.mass_properties["m"], rel=1e-12)
    assert so.mass_properties["I_p"] == pytest.approx(tw.mass_properties["I_p"], rel=1e-12)
    # the 2D web's +-45 plies shear through the thickness: the centre moves by 0.01 % of the depth
    np.testing.assert_allclose(so.shear_centre, tw.shear_centre, atol=2e-4 * 0.3)
    np.testing.assert_allclose(so.tension_centre, tw.tension_centre, atol=1e-4 * 0.3)
    np.testing.assert_allclose(so.geometric_terms, tw.geometric_terms, atol=5e-3)


def test_cas_box_couplings_match_thin_walled(cas_box):
    """Circumferentially asymmetric stiffness box: extension-twist and bending-shear couplings."""
    tw, so = cas_box
    Kt, Ks = tw.stiffness, so.stiffness
    assert abs(Kt[0, 5]) > 0.1 * np.sqrt(Kt[0, 0] * Kt[5, 5])          # strongly coupled
    np.testing.assert_allclose(_scaled(Ks, Kt), 0.0, atol=8e-3)
    # the overlapping junctions of the 2D strips are tied one way round: a 0.1 % asymmetry
    np.testing.assert_allclose(so.shear_centre, 0.0, atol=2e-3 * 0.1)


# ------------------------------------------------------------------ SectionAnalysis and jobs
def test_section_analysis_from_solid(cas_box):
    _, so = cas_box
    sa = SectionAnalysis.from_solid(so, 7)
    assert sa.source == "solid" and sa.section_id == 7
    np.testing.assert_allclose(sa.stiffness, ELEMENT_SIGN @ so.stiffness @ ELEMENT_SIGN)
    assert sa.is_coupled
    assert sa.benscoter_stiffness.shape == (8, 8) and sa.gbt_stiffness is None
    assert sa.geometric_terms.shape == (7, 6)
    assert sa.mass_matrix()[0, 0] == pytest.approx(so.mass_properties["m"])
    # pretwist needs the axial quadrature: sum k r^2 is int E_x r^2
    twisted = sa.with_pretwist(0.5)
    y, z, k = sa.axial_points.T
    assert twisted.stiffness[0, 5] - sa.stiffness[0, 5] == pytest.approx(0.5 * k @ (y * y + z * z))


MESH_2D = (
    "[Section_Node]\n[section_id] [node_id] [y] [z]\n"
    + "".join(f"2 {i + 1} {y} {z}\n" for i, (y, z) in enumerate(
        [(-0.01, -0.02), (0.01, -0.02), (0.01, 0.02), (-0.01, 0.02), (0, -0.02), (0.01, 0), (0, 0.02), (-0.01, 0)]))
    + "[Section_Element]\n[section_id] [element_id] [type] [nodes] [material_id] [theta] [alpha]\n"
    "2 1 Q8 1 2 3 4 5 6 7 8 1 0 0\n"
)
LAMINATE = (
    "[Ply_Material]\n[material_id] [E1] [E2] [G12] [nu12] [G13] [G23] [rho]\n"
    "1 135e9 10e9 5e9 0.3 5e9 3.8e9 1550\n"
    "[Laminate]\n[laminate_id] [ply] [material_id] [angle_deg] [thickness]\n1 1 1 0.0 0.002\n"
)
GEOMETRY = (
    "[Section_Point]\n[section_id] [point_id] [y] [z]\n"
    "1 1 -0.2 -0.1\n1 2 0.2 -0.1\n1 3 0.2 0.1\n1 4 -0.2 0.1\n"
    "[Section_Wall]\n[section_id] [wall_id] [point_start] [point_end] [laminate_id]\n"
    "1 1 1 2 1\n1 2 2 3 1\n1 3 3 4 1\n1 4 4 1 1\n"
)


def _write(d, name, text):
    with open(os.path.join(d, name), "w", encoding="utf-8") as fh:
        fh.write(text)


def test_job_with_2d_sections(tmp_path):
    ed = {"ids": np.array([0, 1, 2])}
    _write(tmp_path, "section_mesh_2d.txt", MESH_2D)
    with pytest.raises(FileNotFoundError, match="element_section.txt"):
        attach_section_analyses(str(tmp_path), ed)
    _write(tmp_path, "element_section.txt", "[Element_Section]\n[element_id] [section_id]\n0 2\n")
    # material 1 is in neither [Solid_Material] nor laminate.txt
    with pytest.raises(ValueError, match="no material"):
        attach_section_analyses(str(tmp_path), ed)
    _write(tmp_path, "laminate.txt", LAMINATE)
    out = attach_section_analyses(str(tmp_path), ed)
    sa = out["section_analysis_per_element"][0]
    assert sa.source == "solid" and sa.section_id == 2
    assert sa.stiffness[0, 0] == pytest.approx(135e9 * 0.02 * 0.04, rel=1e-9)
    assert out["thin_walled_analysis_per_element"] == (None, None, None)
    assert isinstance(out["solid_analysis_per_element"][0], SolidSectionAnalysis)
    # thin-walled sections alongside
    _write(tmp_path, "section_geometry.txt", GEOMETRY)
    _write(tmp_path, "element_section.txt", "[Element_Section]\n[element_id] [section_id]\n0 2\n1 1\n")
    out = attach_section_analyses(str(tmp_path), ed)
    assert [s.source for s in out["section_analysis_per_element"][:2]] == ["solid", "thin_walled"]
    assert out["thin_walled_analysis_per_element"][0] is None
    # GBT needs thin walls
    ed_gbt = {**ed, "types": np.array(["LinearGBTBeamElement3D", "LinearTimoshenkoBeamElement3D",
                                        "LinearTimoshenkoBeamElement3D"])}
    with pytest.raises(ValueError, match="GBT needs a thin-walled section"):
        attach_section_analyses(str(tmp_path), ed_gbt)
    # a section id in both files is refused
    _write(tmp_path, "section_mesh_2d.txt", MESH_2D.replace("\n2 ", "\n1 "))
    with pytest.raises(ValueError, match="both"):
        attach_section_analyses(str(tmp_path), ed)


def test_cantilever_job_with_a_coupled_2d_section(tmp_path):
    """A Timoshenko cantilever on an off-axis solid section: tip deflection from the 2D
    section's coupled flexibility, ``w(L) = int (gamma_xz - theta_y) dx``."""
    import csv
    import shutil

    from workflow_orchestrator.run_job import process_job, setup_job_results_directory

    n, L, P = 40, 1.0, 100.0          # linear elements: tip error 1 / (4 n^2)
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
    mesh = SectionMesh.rectangle(0.02, 0.04, 3, 6, ORTHO, theta=30.0)
    rows_n = "".join(f"5 {i + 1} {y:.17g} {z:.17g}\n" for i, (y, z) in enumerate(mesh.nodes))
    rows_e = "".join(f"5 {k + 1} Q8 {' '.join(str(c + 1) for c in conn)} 3 30 0\n" for k, conn in enumerate(mesh.elements))
    m = ORTHO
    (job / "section_mesh_2d.txt").write_text(
        "[Section_Node]\n[section_id] [node_id] [y] [z]\n" + rows_n
        + "[Section_Element]\n[section_id] [element_id] [type] [nodes] [material_id] [theta] [alpha]\n" + rows_e
        + "[Solid_Material]\n[material_id] [E1] [E2] [E3] [G12] [G13] [G23] [nu12] [nu13] [nu23] [rho]\n"
        f"3 {m.E1} {m.E2} {m.E3} {m.G12} {m.G13} {m.G23} {m.nu12} {m.nu13} {m.nu23} {m.rho}\n")
    (job / "element_section.txt").write_text("[Element_Section]\n[element_id] [section_id]\n"
                                             + "".join(f"{e} 5\n" for e in range(n)))
    (job / "prescribed_displacement.txt").write_text(
        "[Prescribed Displacement]\n[id] [node_id] [dof] [value] [type] [comment]\n"
        + "".join(f"{i} 0 {d} 0.0 displacement\n" for i, d in enumerate(("UX", "UY", "UZ", "RX", "RY", "RZ"))))
    (job / "point_load.txt").write_text("[Point load]\n[x] [y] [z] [F_x] [F_y] [F_z] [M_x] [M_y] [M_z]\n"
                                        f"{L} 0 0 0 0 {P} 0 0 0\n")
    (job / "simulation_settings.txt").write_text(
        "[Simulation]\n[Type]\nstatic\n\n[Solver]\ntype = direct\ntolerance = 1e-12\nmax_iterations = 500\n")
    res = setup_job_results_directory("test_cantilever_job_with_a_coupled_2d_section")
    try:
        process_job(str(job), res, {}, {}, force_serial=True, max_processes_per_job=1)
        with open(os.path.join(res, "primary_results", "global", "U_global.csv"), newline="") as fh:
            U = np.array([float(r[1]) for r in list(csv.reader(fh))[1:]]).reshape(n + 1, -1)
    finally:
        shutil.rmtree(res, ignore_errors=True)
    flex = SolidSectionAnalysis(mesh).flexibility
    assert abs(flex[1, 5]) > 0.1 * np.sqrt(flex[1, 1] * flex[5, 5])      # bending-twist coupled
    # right-handed section forces at x: M_y = -P (L - x), V_z = P
    x = np.linspace(0.0, L, 2001)
    psi = np.outer(-P * (L - x), flex[:, 1]) + P * flex[:, 4]
    theta_y = np.concatenate([[0.0], np.cumsum(0.5 * (psi[1:, 1] + psi[:-1, 1]) * np.diff(x))])
    w_prime = psi[:, 4] - theta_y
    w = np.trapezoid(w_prime, x) if hasattr(np, "trapezoid") else np.trapz(w_prime, x)
    phi = flex[5, 1] * (-P * L * L / 2) + flex[5, 4] * P * L
    assert U[n, 2] == pytest.approx(w, rel=5e-4)
    assert U[n, 3] == pytest.approx(phi, rel=5e-4)

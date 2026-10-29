"""SectionAnalysis: element sign convention, isotropic equivalence and job wiring."""

from __future__ import annotations

import os
import shutil
import tempfile

import numpy as np
import pytest

from pre_processing.element_library.shape_function_registry import get_shape_function_operator
from pre_processing.parsing.element_section_parser import parse_element_section_file
from pre_processing.section_integration.section_analysis import (
    ELEMENT_SIGN,
    SectionAnalysis,
    to_element_voigt,
)
from pre_processing.section_integration.thin_walled import (
    Ply,
    PlyMaterial,
    ThinWalledSection,
    ThinWalledSectionAnalysis,
    Wall,
    WallLaminate,
)
from pre_processing.section_integration.thin_walled.job_input import attach_section_analyses

E, NU = 70e9, 0.3
G = E / (2 * (1 + NU))
CFRP = PlyMaterial(E1=140e9, E2=10e9, G12=5e9, nu12=0.3, G13=5e9, G23=3.5e9, rho=1600.0)


def _tw(points, walls, lams, **kw):
    return ThinWalledSectionAnalysis(ThinWalledSection(points, walls), lams, **kw)


# --------------------------------------------------------------------------- #
# Element sign convention
# --------------------------------------------------------------------------- #


def _element_B(element_type, L=2.0, xi=0.3):
    from pre_processing.element_library.linear.beam.first_order_shear_deformation_theory.timoshenko.utilities.B_matrix import (
        StrainDisplacementOperator as TimoshenkoB,
    )
    from pre_processing.element_library.linear.beam.zero_order_shear_deformation_theory.euler_bernoulli.utilities.B_matrix import (
        StrainDisplacementOperator as EulerBernoulliB,
    )

    N, dN, d2N = get_shape_function_operator(element_type, L).natural_coordinate_form(np.array([xi]))
    if element_type == "LinearTimoshenkoBeamElement3D":
        return TimoshenkoB(L).physical_coordinate_form(dN, d2N, N)[0]
    return EulerBernoulliB(L).physical_coordinate_form(dN, d2N)[0]


@pytest.mark.parametrize("etype", ["LinearEulerBernoulliBeamElement3D", "LinearTimoshenkoBeamElement3D"])
def test_element_strains_are_sign_mapped_right_handed_strains(etype):
    """Nodal fields with known right-handed strains give S @ psi through the element B."""
    L = 2.0
    # mid-point: linear Timoshenko interpolation represents these fields exactly there
    B = _element_B(etype, L, xi=0.0)
    # u_z = x^2 / 2, theta_y = -u_z' = -x (right hand): kappa_y_rh = -1
    U = np.zeros(12)
    U[8], U[10] = L * L / 2, -L
    psi_rh = np.array([0, -1.0, 0, 0, 0, 0])
    np.testing.assert_allclose(B @ U, ELEMENT_SIGN @ psi_rh, atol=1e-12)
    # u_y = x^2 / 2, theta_z = u_y' = x: kappa_z_rh = +1
    U = np.zeros(12)
    U[7], U[11] = L * L / 2, L
    np.testing.assert_allclose(B @ U, ELEMENT_SIGN @ np.array([0, 0, 1.0, 0, 0, 0]), atol=1e-12)
    # rigid right-handed rotations give no strain
    for U in (np.r_[0, 0, 0, 0, 1, 0, 0, 0, -L, 0, 1, 0], np.r_[0, 0, 0, 0, 0, 1, 0, L, 0, 0, 0, 1]):
        np.testing.assert_allclose(B @ U, 0.0, atol=1e-12)


def test_timoshenko_shear_strain_sign():
    L = 2.0
    B = _element_B("LinearTimoshenkoBeamElement3D", L)
    # u_z = x with theta_y = 0: gamma_xz_rh = w' + theta_y = 1
    U = np.zeros(12)
    U[8] = L
    np.testing.assert_allclose(B @ U, [0, 0, 0, 0, 1.0, 0], atol=1e-12)


def test_offset_strip_axial_bending_coupling_sign():
    """A box centred at z = z0: eps_xx = eps - z kappa_y_element, so K_e[0, 1] = -EA z0."""
    z0, b, h, t = 0.05, 0.3, 0.1, 2e-3
    pts = {1: (-b / 2, z0 - h / 2), 2: (b / 2, z0 - h / 2), 3: (b / 2, z0 + h / 2), 4: (-b / 2, z0 + h / 2)}
    sa = _tw(pts, [Wall(i, i, i % 4 + 1, 1) for i in range(1, 5)], {1: WallLaminate.isotropic(E, NU, t)}, hoop="free")
    se = SectionAnalysis.from_thin_walled(sa, section_id=7)
    EA = E * 2 * (b + h) * t
    assert se.stiffness[0, 1] == pytest.approx(-EA * z0, rel=1e-10)
    assert se.classical_stiffness[0, 1] == pytest.approx(-EA * z0, rel=1e-10)
    assert sa.stiffness[0, 1] == pytest.approx(EA * z0, rel=1e-10)
    assert se.source == "thin_walled" and se.section_id == 7
    assert se.is_coupled


def test_sign_map_preserves_energy():
    lam = WallLaminate([Ply(CFRP, 30.0, 1e-3), Ply(CFRP, -10.0, 2e-3)])
    pts = {1: (-0.2, -0.1), 2: (0.2, -0.1), 3: (0.2, 0.1), 4: (-0.2, 0.1)}
    walls = [Wall(i, i, i % 4 + 1, 1) for i in range(1, 5)]
    sa = _tw(pts, walls, {1: lam})
    Ke = to_element_voigt(sa.stiffness)
    rng = np.random.default_rng(1)
    for _ in range(3):
        psi = rng.standard_normal(6)
        assert (ELEMENT_SIGN @ psi) @ Ke @ (ELEMENT_SIGN @ psi) == pytest.approx(psi @ sa.stiffness @ psi, rel=1e-12)


# --------------------------------------------------------------------------- #
# Isotropic equivalence
# --------------------------------------------------------------------------- #


def test_isotropic_section_analysis_matches_timoshenko_D():
    from pre_processing.element_library.linear.beam.first_order_shear_deformation_theory.timoshenko.utilities.D_matrix import (
        MaterialStiffnessOperator,
    )

    props = dict(A=0.00131, I_y=3.234e-7, I_z=2.08769e-6, J_t=2.60673e-8)
    D = MaterialStiffnessOperator(
        youngs_modulus=2.1e11, shear_modulus=8.1e10, cross_section_area=props["A"],
        moment_inertia_y=props["I_y"], moment_inertia_z=props["I_z"], torsion_constant=props["J_t"],
        shear_correction_factor=0.8,
    ).assembly_form()
    sa = SectionAnalysis.from_isotropic(E=2.1e11, G=8.1e10, kappa=0.8, **props)
    np.testing.assert_array_equal(sa.stiffness, D)
    assert not sa.is_coupled


def test_thin_walled_isotropic_rectangle_tube_matches_line_formulas():
    b, h, t = 0.4, 0.2, 2e-3
    pts = {1: (-b / 2, -h / 2), 2: (b / 2, -h / 2), 3: (b / 2, h / 2), 4: (-b / 2, h / 2)}
    walls = [Wall(i, i, i % 4 + 1, 1) for i in range(1, 5)]
    se = SectionAnalysis.from_thin_walled(_tw(pts, walls, {1: WallLaminate.isotropic(E, NU, t)}, hoop="free"))
    I_y = 2 * b * t * (h / 2) ** 2 + 2 * t * h**3 / 12 + 2 * b * t**3 / 12
    assert se.stiffness[1, 1] == pytest.approx(E * I_y, rel=1e-10)
    assert se.stiffness[0, 0] == pytest.approx(E * 2 * (b + h) * t, rel=1e-12)
    assert not se.is_coupled


def test_section_analysis_validates_shapes():
    with pytest.raises(ValueError):
        SectionAnalysis(np.eye(5), np.eye(4))
    K = np.eye(6)
    K[0, 1] = 1.0
    with pytest.raises(ValueError):
        SectionAnalysis(K, np.eye(4))


# --------------------------------------------------------------------------- #
# Inputs and wiring
# --------------------------------------------------------------------------- #

LAMINATE = (
    "[Ply_Material]\n[material_id] [E1] [E2] [G12] [nu12] [G13] [G23]\n"
    f"1 {E} {E} {G} {NU} {G} {G}\n"
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


def test_parse_element_section_file(tmp_path):
    _write(tmp_path, "element_section.txt", "# x\n[Element_Section]\n[element_id] [section_id]\n0 1\n3 2\n")
    assert parse_element_section_file(str(tmp_path / "element_section.txt")) == {0: 1, 3: 2}
    _write(tmp_path, "bad.txt", "[Element_Section]\n[element_id] [section_id]\n0 1\n0 2\n")
    with pytest.raises(ValueError, match="Duplicate"):
        parse_element_section_file(str(tmp_path / "bad.txt"))


def test_attach_section_analyses(tmp_path):
    ed = {"ids": np.array([0, 1, 2])}
    assert attach_section_analyses(str(tmp_path), ed) is ed
    _write(tmp_path, "laminate.txt", LAMINATE)
    with pytest.raises(FileNotFoundError, match="section_geometry.txt"):
        attach_section_analyses(str(tmp_path), ed)
    _write(tmp_path, "section_geometry.txt", GEOMETRY)
    _write(tmp_path, "element_section.txt", "[Element_Section]\n[element_id] [section_id]\n0 1\n2 1\n")
    out = attach_section_analyses(str(tmp_path), ed)
    per = out["section_analysis_per_element"]
    assert per[1] is None
    assert per[0] is per[2]
    assert per[0].stiffness[0, 0] == pytest.approx(E * 1.2 * 0.002, rel=1e-6)
    _write(tmp_path, "element_section.txt", "[Element_Section]\n[element_id] [section_id]\n5 1\n")
    with pytest.raises(ValueError, match="unknown elements"):
        attach_section_analyses(str(tmp_path), ed)
    _write(tmp_path, "element_section.txt", "[Element_Section]\n[element_id] [section_id]\n0 9\n")
    with pytest.raises(ValueError, match="unknown sections"):
        attach_section_analyses(str(tmp_path), ed)


def _element_kwargs(job_results_dir, etype):
    for sub in ("element_stiffness_matrices", "element_force_vectors", "logs"):
        os.makedirs(os.path.join(job_results_dir, sub), exist_ok=True)
    ed = {
        "ids": np.array([0]),
        "connectivity": np.array([[0, 1]]),
        "types": np.array([etype]),
        "integration_orders": {k: np.array([2]) for k in
                               ("axial", "bending_y", "bending_z", "shear_y", "shear_z", "torsion", "load")},
    }
    return dict(
        element_id=0,
        element_dictionary=ed,
        grid_dictionary={"coordinates": np.array([[0.0, 0, 0], [1.0, 0, 0]])},
        material_dictionary={k: np.array([v]) for k, v in dict(E=2.1e11, G=8.1e10, nu=0.3, rho=7850.0).items()},
        section_dictionary={k: np.array([v]) for k, v in dict(
            A=1e-3, I_x=0.0, I_y=1e-8, I_z=2e-8, J_t=1e-9, kappa=0.8, alpha=0.0).items()},
        point_load_array=np.empty((0, 9)),
        distributed_load_array=np.empty((0, 9)),
        job_results_dir=job_results_dir,
    )


def test_base_class_isotropic_section_analysis_and_unsupported_assignment():
    from pre_processing.element_library.linear.beam.first_order_shear_deformation_theory.timoshenko.linear_timoshenko_3D import (
        LinearTimoshenkoBeamElement3D,
    )

    tmp = tempfile.mkdtemp()
    try:
        kw = _element_kwargs(os.path.join(tmp, "job"), "LinearTimoshenkoBeamElement3D")
        el = LinearTimoshenkoBeamElement3D(**kw)
        assert not el.has_thin_walled_section
        np.testing.assert_allclose(el.section_analysis.stiffness, el.D_matrix(), rtol=1e-15)
        assert el.section_analysis.mass["m"] == pytest.approx(7850.0 * 1e-3)
        se = SectionAnalysis.from_isotropic(E=1.0, G=1.0, A=1.0, I_y=1.0, I_z=1.0, J_t=1.0)
        se = SectionAnalysis(se.stiffness, se.classical_stiffness, source="thin_walled", section_id=3)
        kw["element_dictionary"] = {**kw["element_dictionary"], "section_analysis_per_element": (se,)}
        if LinearTimoshenkoBeamElement3D.supports_section_analysis:
            assert LinearTimoshenkoBeamElement3D(**kw).section_analysis is se
        else:
            with pytest.raises(NotImplementedError, match="section 3"):
                LinearTimoshenkoBeamElement3D(**kw)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

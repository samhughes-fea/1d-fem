"""
Unit tests for bar and truss utility suite (shape_functions, B_matrix, D_matrix, interpolate_loads).
"""

import sys
from pathlib import Path

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from pre_processing.element_library.linear.bar.utilities import (
    ShapeFunctionOperator as BarShapeFunctionOperator,
    StrainDisplacementOperator as BarStrainDisplacementOperator,
    MaterialStiffnessOperator as BarMaterialStiffnessOperator,
    LoadInterpolationOperator as BarLoadInterpolationOperator,
)
from pre_processing.element_library.linear.truss.utilities import (
    ShapeFunctionOperator as TrussShapeFunctionOperator,
    StrainDisplacementOperator as TrussStrainDisplacementOperator,
    MaterialStiffnessOperator as TrussMaterialStiffnessOperator,
    LoadInterpolationOperator as TrussLoadInterpolationOperator,
)


# --- Bar utilities ---


def test_bar_shape_functions_linear_lagrange():
    """Bar shape functions: linear Lagrange at ξ=-1,1; only axial and torsion DOF non-zero."""
    L = 1.0
    op = BarShapeFunctionOperator(element_length=L)
    xi = np.array([-1.0, 0.0, 1.0])
    N, dN_dξ, d2N_dξ2 = op.natural_coordinate_form(xi)
    assert N.shape == (3, 12, 6)
    assert dN_dξ.shape == (3, 12, 6)
    assert np.all(d2N_dξ2 == 0.0)
    # At ξ=-1: N1=1, N2=0 for axial (DOF 0,6) and torsion (3,9)
    np.testing.assert_allclose(N[0, 0, 0], 1.0)
    np.testing.assert_allclose(N[0, 6, 0], 0.0)
    np.testing.assert_allclose(N[0, 3, 3], 1.0)
    np.testing.assert_allclose(N[0, 9, 3], 0.0)
    # At ξ=1: N1=0, N2=1
    np.testing.assert_allclose(N[2, 0, 0], 0.0)
    np.testing.assert_allclose(N[2, 6, 0], 1.0)
    np.testing.assert_allclose(N[2, 3, 3], 0.0)
    np.testing.assert_allclose(N[2, 9, 3], 1.0)
    # dN/dξ constant: -0.5 for node1, +0.5 for node2
    np.testing.assert_allclose(dN_dξ[:, 0, 0], -0.5)
    np.testing.assert_allclose(dN_dξ[:, 6, 0], 0.5)
    assert op.dξ_dx == 2.0 / L


def test_bar_B_matrix_constant_and_axial_torsion():
    """Bar B-matrix is (2, 12); row 0 axial, row 1 torsion; constant in ξ."""
    L = 2.0
    axial = np.array([1.0, 0.0, 0.0])
    op = BarStrainDisplacementOperator(element_length=L, axial=axial)
    dN_dξ = np.zeros((1, 12, 6))
    dN_dξ[0, 0, 0], dN_dξ[0, 6, 0] = -0.5, 0.5
    dN_dξ[0, 3, 3], dN_dξ[0, 9, 3] = -0.5, 0.5
    B = op.physical_coordinate_form(dN_dξ)[0]
    assert B.shape == (2, 12)
    # ε_axial: (1/L)(-ux1 + ux2) -> B[0,0]=-1/L, B[0,6]=1/L
    np.testing.assert_allclose(B[0, 0], -1.0 / L)
    np.testing.assert_allclose(B[0, 6], 1.0 / L)
    np.testing.assert_allclose(B[0, 1:6], 0.0)
    np.testing.assert_allclose(B[0, 7:12], 0.0)
    # φ_torsion: B[1,3]=-1/L, B[1,9]=1/L
    np.testing.assert_allclose(B[1, 3], -1.0 / L)
    np.testing.assert_allclose(B[1, 9], 1.0 / L)
    np.testing.assert_allclose(op.jacobian, L / 2)


def test_bar_D_matrix_diagonal_and_stress_resultants():
    """Bar D-matrix is 2×2 diag(EA, GJ_t); stress resultants and energy."""
    E, G = 2.1e11, 8.1e10
    A, J_t = 0.001, 1e-8
    op = BarMaterialStiffnessOperator(
        youngs_modulus=E,
        shear_modulus=G,
        cross_section_area=A,
        torsion_constant=J_t,
    )
    D = op.assembly_form()
    assert D.shape == (2, 2)
    np.testing.assert_allclose(D[0, 0], E * A)
    np.testing.assert_allclose(D[1, 1], G * J_t)
    np.testing.assert_allclose(D[0, 1], 0.0)
    np.testing.assert_allclose(D[1, 0], 0.0)
    assert op.postprocessing_form() is not op.assembly_form()
    np.testing.assert_allclose(op.postprocessing_form(), D)
    strain = np.array([1e-4, 0.01])
    stress = op.compute_stress_resultants(strain)
    np.testing.assert_allclose(stress[0], E * A * strain[0])
    np.testing.assert_allclose(stress[1], G * J_t * strain[1])
    energy = op.energy_density_components(strain)
    assert "total" in energy and "axial" in energy and "torsion" in energy


def test_bar_load_interpolation():
    """Bar LoadInterpolationOperator: (N, 9) array, interpolate returns (M, 6)."""
    loads = np.array([
        [0.0, 0.0, 0.0, 10.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0, 20.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    ])
    op = BarLoadInterpolationOperator(
        distributed_loads_array=loads,
        boundary_mode="clamp",
        interpolation_order="linear",
        n_gauss_points=1,
    )
    q = op.interpolate(0.5)
    np.testing.assert_allclose(q[0], 15.0)
    q_multi = op.interpolate(np.array([0.0, 0.5, 1.0]))
    assert q_multi.shape == (3, 6)


# --- Truss utilities ---


def test_truss_shape_functions_linear_lagrange():
    """Truss shape functions: linear Lagrange for axial, transverse, torsion DOFs."""
    L = 1.0
    op = TrussShapeFunctionOperator(element_length=L)
    xi = np.array([-1.0, 1.0])
    N, dN_dξ, d2N_dξ2 = op.natural_coordinate_form(xi)
    assert N.shape == (2, 12, 6)
    assert np.all(d2N_dξ2 == 0.0)
    # At ξ=-1: node1=1, node2=0 for u_x(0,6), u_y(1,7), u_z(2,8), θ_x(3,9)
    np.testing.assert_allclose(N[0, [0, 1, 2, 3], [0, 1, 2, 3]], [1, 1, 1, 1])
    np.testing.assert_allclose(N[0, [6, 7, 8, 9], [0, 1, 2, 3]], [0, 0, 0, 0])
    np.testing.assert_allclose(N[1, [0, 1, 2, 3], [0, 1, 2, 3]], [0, 0, 0, 0])
    np.testing.assert_allclose(N[1, [6, 7, 8, 9], [0, 1, 2, 3]], [1, 1, 1, 1])


def test_truss_B_matrix_three_rows():
    """Truss B-matrix is (3, 12): axial, transverse, torsion rows."""
    L = 2.0
    axial = np.array([1.0, 0.0, 0.0])
    transverse = np.array([0.0, 1.0, 0.0])
    op = TrussStrainDisplacementOperator(
        element_length=L, axial=axial, transverse=transverse
    )
    dN_dξ = np.zeros((1, 12, 6))
    B = op.physical_coordinate_form(dN_dξ)[0]
    assert B.shape == (3, 12)
    np.testing.assert_allclose(B[0, 0], -1.0 / L)
    np.testing.assert_allclose(B[0, 6], 1.0 / L)
    np.testing.assert_allclose(B[1, 1], -1.0 / L)
    np.testing.assert_allclose(B[1, 7], 1.0 / L)
    np.testing.assert_allclose(B[2, 3], -1.0 / L)
    np.testing.assert_allclose(B[2, 9], 1.0 / L)


def test_truss_D_matrix_three_by_three():
    """Truss D-matrix is 3×3 diag(EA, κGA, GJ_t)."""
    E, G = 2.1e11, 8.1e10
    A, J_t = 0.001, 1e-8
    kappa = 5.0 / 6.0
    op = TrussMaterialStiffnessOperator(
        youngs_modulus=E,
        shear_modulus=G,
        cross_section_area=A,
        torsion_constant=J_t,
        shear_correction_factor=kappa,
    )
    D = op.assembly_form()
    assert D.shape == (3, 3)
    np.testing.assert_allclose(D[0, 0], E * A)
    np.testing.assert_allclose(D[1, 1], kappa * G * A)
    np.testing.assert_allclose(D[2, 2], G * J_t)
    strain = np.array([1e-4, 0.001, 0.01])
    stress = op.compute_stress_resultants(strain)
    np.testing.assert_allclose(stress, D @ strain)
    energy = op.energy_density_components(strain)
    assert "total" in energy and "transverse" in energy


def test_truss_load_interpolation():
    """Truss LoadInterpolationOperator same API as bar."""
    loads = np.array([
        [0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0],
    ])
    op = TrussLoadInterpolationOperator(
        distributed_loads_array=loads,
        boundary_mode="clamp",
        interpolation_order="linear",
        n_gauss_points=1,
    )
    q = op.interpolate(0.5)
    assert q.shape == (6,)
    np.testing.assert_allclose(q[1], 1.5)


def test_bar_element_stiffness_gauss_data_uses_operators():
    """Bar element_stiffness_matrix gauss_data has B, D, N from utilities."""
    import os
    import tempfile
    from pre_processing.element_library.linear.bar.linear_bar_3D import LinearBarElement3D

    L = 1.0
    node_coords = np.array([[0.0, 0.0, 0.0], [L, 0.0, 0.0]], dtype=np.float64)
    element_dictionary = {
        "ids": np.array([0]),
        "connectivity": np.array([[0, 1]]),
        "types": np.array(["LinearBarElement3D"]),
        "integration_orders": {k: np.array([1]) for k in ["axial", "bending_y", "bending_z", "shear_y", "shear_z", "torsion", "load"]},
    }
    grid_dictionary = {"coordinates": node_coords}
    section_dictionary = {"A": np.array([0.001]), "I_x": np.array([0.0]), "I_y": np.array([0.0]), "I_z": np.array([0.0]), "J_t": np.array([1e-8])}
    material_dictionary = {"E": np.array([2.1e11]), "G": np.array([8.1e10]), "nu": np.array([0.3]), "rho": np.array([7850.0])}
    job_dir = tempfile.mkdtemp()
    os.makedirs(os.path.join(job_dir, "element_stiffness_matrices"), exist_ok=True)
    os.makedirs(os.path.join(job_dir, "element_force_vectors"), exist_ok=True)

    el = LinearBarElement3D(
        element_id=0,
        element_dictionary=element_dictionary,
        grid_dictionary=grid_dictionary,
        section_dictionary=section_dictionary,
        material_dictionary=material_dictionary,
        point_load_array=np.empty((0, 9)),
        distributed_load_array=np.empty((0, 9)),
        job_results_dir=job_dir,
    )
    obj = el.element_stiffness_matrix()
    assert len(obj.gauss_data) == 1
    gp = obj.gauss_data[0]
    assert gp.B_matrix.shape == (2, 12)
    assert gp.D_matrix.shape == (2, 2)
    assert gp.shape_functions is not None and gp.shape_functions.shape == (12, 6)
    assert gp.shape_derivatives is not None and gp.shape_derivatives.shape == (12, 6)
    assert obj.evaluate_shape_functions is not None
    N, dN, _ = obj.evaluate_shape_functions(np.array([0.0]))
    assert N.shape == (1, 12, 6)
    assert obj.shape_function_N_coefficients is not None
    assert obj.shape_function_N_coefficients.shape == (12, 6, 4)
    assert obj.integration_scheme == "Gauss-Legendre"


def test_bar_Ke_equals_Lt_Klocal_L():
    """Bar K_e from ∫ Bᵀ D B dx (quadrature) matches analytical Lᵀ K_local L."""
    import os
    import tempfile
    from pre_processing.element_library.linear.bar.linear_bar_3D import LinearBarElement3D
    from pre_processing.element_library.linear.bar.utilities import build_L_matrix_4x12, direction_cosines

    L = 2.0
    E, G = 2.1e11, 8.1e10
    A, J_t = 0.001, 1e-8
    node_coords = np.array([[0.0, 0.0, 0.0], [L, 0.0, 0.0]], dtype=np.float64)
    element_dictionary = {
        "ids": np.array([0]),
        "connectivity": np.array([[0, 1]]),
        "types": np.array(["LinearBarElement3D"]),
        "integration_orders": {k: np.array([1]) for k in ["axial", "bending_y", "bending_z", "shear_y", "shear_z", "torsion", "load"]},
    }
    grid_dictionary = {"coordinates": node_coords}
    section_dictionary = {"A": np.array([A]), "I_x": np.array([0.0]), "I_y": np.array([0.0]), "I_z": np.array([0.0]), "J_t": np.array([J_t])}
    material_dictionary = {"E": np.array([E]), "G": np.array([G]), "nu": np.array([0.3]), "rho": np.array([7850.0])}
    job_dir = tempfile.mkdtemp()
    os.makedirs(os.path.join(job_dir, "element_stiffness_matrices"), exist_ok=True)
    os.makedirs(os.path.join(job_dir, "element_force_vectors"), exist_ok=True)

    el = LinearBarElement3D(
        element_id=0,
        element_dictionary=element_dictionary,
        grid_dictionary=grid_dictionary,
        section_dictionary=section_dictionary,
        material_dictionary=material_dictionary,
        point_load_array=np.empty((0, 9)),
        distributed_load_array=np.empty((0, 9)),
        job_results_dir=job_dir,
    )
    obj = el.element_stiffness_matrix()
    Ke_quad = obj.K_e

    axial = direction_cosines(node_coords)
    L_mat = build_L_matrix_4x12(axial)
    k_axial = (E * A / L) * np.array([[1, -1], [-1, 1]], dtype=np.float64)
    k_torsion = (G * J_t / L) * np.array([[1, -1], [-1, 1]], dtype=np.float64)
    K_local = np.zeros((4, 4))
    K_local[0:2, 0:2] = k_axial
    K_local[2:4, 2:4] = k_torsion
    Ke_ref = L_mat.T @ K_local @ L_mat
    np.testing.assert_allclose(Ke_quad, Ke_ref, rtol=1e-10, err_msg="Bar K_e from ∫ Bᵀ D B dx should match Lᵀ K_local L")


def test_truss_Ke_equals_Lt_Klocal_L():
    """Truss K_e from ∫ Bᵀ D B dx (quadrature) matches analytical Lᵀ K_local L."""
    import os
    import tempfile
    from pre_processing.element_library.linear.truss.linear_truss_3D import LinearTrussElement3D
    from pre_processing.element_library.linear.truss.utilities import (
        direction_cosines_and_transverse,
        build_L_matrix_6x12,
    )

    L = 2.0
    E, G = 2.1e11, 8.1e10
    A, J_t = 0.001, 1e-8
    kappa = 5.0 / 6.0
    node_coords = np.array([[0.0, 0.0, 0.0], [L, 0.0, 0.0]], dtype=np.float64)
    element_dictionary = {
        "ids": np.array([0]),
        "connectivity": np.array([[0, 1]]),
        "types": np.array(["LinearTrussElement3D"]),
        "integration_orders": {k: np.array([1]) for k in ["axial", "bending_y", "bending_z", "shear_y", "shear_z", "torsion", "load"]},
    }
    grid_dictionary = {"coordinates": node_coords}
    section_dictionary = {"A": np.array([A]), "I_x": np.array([0.0]), "I_y": np.array([0.0]), "I_z": np.array([0.0]), "J_t": np.array([J_t])}
    material_dictionary = {"E": np.array([E]), "G": np.array([G]), "nu": np.array([0.3]), "rho": np.array([7850.0])}
    job_dir = tempfile.mkdtemp()
    os.makedirs(os.path.join(job_dir, "element_stiffness_matrices"), exist_ok=True)
    os.makedirs(os.path.join(job_dir, "element_force_vectors"), exist_ok=True)

    el = LinearTrussElement3D(
        element_id=0,
        element_dictionary=element_dictionary,
        grid_dictionary=grid_dictionary,
        section_dictionary=section_dictionary,
        material_dictionary=material_dictionary,
        point_load_array=np.empty((0, 9)),
        distributed_load_array=np.empty((0, 9)),
        job_results_dir=job_dir,
    )
    obj = el.element_stiffness_matrix()
    Ke_quad = obj.K_e

    axial, transverse = direction_cosines_and_transverse(node_coords)
    L_mat = build_L_matrix_6x12(axial, transverse)
    k_axial = (E * A / L) * np.array([[1, -1], [-1, 1]], dtype=np.float64)
    k_trans = (kappa * G * A / L) * np.array([[1, -1], [-1, 1]], dtype=np.float64)
    k_torsion = (G * J_t / L) * np.array([[1, -1], [-1, 1]], dtype=np.float64)
    K_local = np.zeros((6, 6))
    K_local[0:2, 0:2] = k_axial
    K_local[2:4, 2:4] = k_trans
    K_local[4:6, 4:6] = k_torsion
    Ke_ref = L_mat.T @ K_local @ L_mat
    np.testing.assert_allclose(Ke_quad, Ke_ref, rtol=1e-10, err_msg="Truss K_e from ∫ Bᵀ D B dx should match Lᵀ K_local L")

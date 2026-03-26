"""
Unit tests for B2 shape function coefficient evaluator.

Validates that evaluate_shape_functions_from_coefficients returns (N, dN_dξ, d2N_dξ2)
consistent with the monomial convention and matches EB ShapeFunctionOperator for known coefficients.
"""

import numpy as np
import pytest

from pre_processing.element_library.utilities.shape_function_coefficient_evaluator import (
    evaluate_shape_functions_from_coefficients,
)


def test_evaluator_n1_xi_hand_values():
    """N1(ξ) = 0.5 - 0.75*ξ + 0.25*ξ³: coeffs (1,1) = [0.5, -0.75, 0, 0.25]. Check at ξ=-1,0,1."""
    N_coeffs = np.zeros((12, 6, 4))
    dN_coeffs = np.zeros((12, 6, 4))
    d2N_coeffs = np.zeros((12, 6, 4))
    # N1 at dof 1, comp 1: 0.5 - 0.75*ξ + 0.25*ξ³
    N_coeffs[1, 1, 0] = 0.5
    N_coeffs[1, 1, 1] = -0.75
    N_coeffs[1, 1, 3] = 0.25
    # dN1/dξ = -0.75 + 0.75*ξ²  (coeffs for ξ^0 and ξ^2)
    dN_coeffs[1, 1, 0] = -0.75
    dN_coeffs[1, 1, 2] = 0.75
    # d²N1/dξ² = 1.5*ξ
    d2N_coeffs[1, 1, 1] = 1.5

    xi = np.array([-1.0, 0.0, 1.0])
    N, dN_dxi, d2N_dxi2 = evaluate_shape_functions_from_coefficients(
        N_coeffs, dN_coeffs, d2N_coeffs, xi
    )
    # N1(-1)=0.5+0.75-0.25=1, N1(0)=0.5, N1(1)=0.5-0.75+0.25=0
    np.testing.assert_allclose(N[:, 1, 1], [1.0, 0.5, 0.0])
    # dN1/dξ(-1)=-0.75+0.75=0, dN1/dξ(0)=-0.75, dN1/dξ(1)=-0.75+0.75=0
    np.testing.assert_allclose(dN_dxi[:, 1, 1], [0.0, -0.75, 0.0])
    # d²N1/dξ²(-1)=-1.5, d²N1/dξ²(0)=0, d²N1/dξ²(1)=1.5
    np.testing.assert_allclose(d2N_dxi2[:, 1, 1], [-1.5, 0.0, 1.5])


def test_evaluator_matches_eb_operator():
    """Full EB coefficient arrays (from EB export) should match ShapeFunctionOperator at sample ξ."""
    from pre_processing.element_library.linear.euler_bernoulli.linear_euler_bernoulli_3D import (
        LinearEulerBernoulliBeamElement3D,
    )
    import tempfile
    import os
    import shutil

    L = 1.0
    element_dictionary = {
        "ids": np.array([0]),
        "connectivity": np.array([[0, 1]]),
        "types": np.array(["LinearEulerBernoulliBeamElement3D"]),
        "integration_orders": {
            "axial": np.array([3]),
            "bending_y": np.array([3]),
            "bending_z": np.array([3]),
            "shear_y": np.array([2]),
            "shear_z": np.array([2]),
            "torsion": np.array([3]),
            "load": np.array([2]),
        },
    }
    grid_dictionary = {"coordinates": np.array([[0.0, 0.0, 0.0], [L, 0.0, 0.0]])}
    material_dictionary = {
        "E": np.array([2.1e11]),
        "G": np.array([8.1e10]),
        "nu": np.array([0.3]),
        "rho": np.array([7850.0]),
    }
    section_dictionary = {
        "A": np.array([1e-3]),
        "I_x": np.array([1e-9]),
        "I_y": np.array([1e-8]),
        "I_z": np.array([1e-8]),
        "J_t": np.array([1e-9]),
    }
    temp_dir = tempfile.mkdtemp()
    job_results_dir = os.path.join(temp_dir, "test_results")
    os.makedirs(job_results_dir, exist_ok=True)
    os.makedirs(os.path.join(job_results_dir, "element_stiffness_matrices"), exist_ok=True)
    os.makedirs(os.path.join(job_results_dir, "element_force_vectors"), exist_ok=True)
    try:
        element = LinearEulerBernoulliBeamElement3D(
            element_id=0,
            element_dictionary=element_dictionary,
            grid_dictionary=grid_dictionary,
            material_dictionary=material_dictionary,
            section_dictionary=section_dictionary,
            point_load_array=np.empty((0, 9)),
            distributed_load_array=np.array([[0, 0, 0, 0, 0, 0, 0, 0, 0], [L, 0, 0, 0, 0, 0, 0, 0, 0]]),
            job_results_dir=job_results_dir,
        )
        elem_obj = element.element_stiffness_matrix()
        assert elem_obj.shape_function_N_coefficients is not None
        assert elem_obj.shape_function_dN_dxi_coefficients is not None
        assert elem_obj.shape_function_d2N_dxi2_coefficients is not None

        xi = np.array([-0.8, -0.3, 0.0, 0.5, 1.0])
        N_coeff, dN_coeff, d2N_coeff = evaluate_shape_functions_from_coefficients(
            elem_obj.shape_function_N_coefficients,
            elem_obj.shape_function_dN_dxi_coefficients,
            elem_obj.shape_function_d2N_dxi2_coefficients,
            xi,
        )
        N_op, dN_op, d2N_op = element.shape_function_operator.natural_coordinate_form(xi)
        np.testing.assert_allclose(N_coeff, N_op, rtol=1e-12, atol=1e-12)
        np.testing.assert_allclose(dN_coeff, dN_op, rtol=1e-12, atol=1e-12)
        np.testing.assert_allclose(d2N_coeff, d2N_op, rtol=1e-12, atol=1e-12)
    finally:
        if os.path.isdir(temp_dir):
            shutil.rmtree(temp_dir, ignore_errors=True)

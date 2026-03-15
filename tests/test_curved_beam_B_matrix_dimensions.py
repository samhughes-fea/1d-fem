"""
Curved beam B-matrix: shape (6, 12), DOF mapping, and curvature terms when κ0 != 0.
"""

import os
import sys
from pathlib import Path

import numpy as np
import pytest

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from pre_processing.element_library.linear.curved_beam.linear_curved_timoshenko_3D import (
    LinearCurvedTimoshenkoBeamElement3D,
)


def test_curved_beam_B_matrix_dimensions():
    """Curved element B at a Gauss point has shape (6, 12) and correct DOF mapping."""
    L = 1.0
    element_dictionary = {
        "ids": np.array([0]),
        "connectivity": np.array([[0, 1]]),
        "types": np.array(["LinearCurvedTimoshenkoBeamElement3D"]),
        "integration_orders": {
            "axial": np.array([2]),
            "bending_y": np.array([2]),
            "bending_z": np.array([2]),
            "shear_y": np.array([2]),
            "shear_z": np.array([2]),
            "torsion": np.array([2]),
            "load": np.array([2]),
        },
        "curvature": np.array([0.5]),
    }
    grid_dictionary = {"coordinates": np.array([[0.0, 0.0, 0.0], [L, 0.0, 0.0]])}
    material_dictionary = {
        "E": np.array([2.1e11]),
        "G": np.array([8.1e10]),
        "nu": np.array([0.3]),
        "rho": np.array([7850.0]),
    }
    section_dictionary = {
        "A": np.array([0.001]),
        "I_x": np.array([0.0]),
        "I_y": np.array([1e-8]),
        "I_z": np.array([1e-8]),
        "J_t": np.array([1e-9]),
        "kappa": np.array([5.0 / 6.0]),
    }
    point_load_array = np.empty((0, 9))
    distributed_load_array = np.empty((0, 9))

    import tempfile
    temp_dir = tempfile.mkdtemp()
    job_results_dir = os.path.join(temp_dir, "curved_b")
    os.makedirs(job_results_dir, exist_ok=True)
    os.makedirs(os.path.join(job_results_dir, "element_stiffness_matrices"), exist_ok=True)
    os.makedirs(os.path.join(job_results_dir, "element_force_vectors"), exist_ok=True)

    try:
        element = LinearCurvedTimoshenkoBeamElement3D(
            element_id=0,
            element_dictionary=element_dictionary,
            grid_dictionary=grid_dictionary,
            material_dictionary=material_dictionary,
            section_dictionary=section_dictionary,
            point_load_array=point_load_array,
            distributed_load_array=distributed_load_array,
            job_results_dir=job_results_dir,
        )
        xi = np.array([0.0])
        N, dN_dξ, d2N_dξ2 = element.shape_function_operator.natural_coordinate_form(xi)
        B = element.strain_displacement_operator.physical_coordinate_form(dN_dξ, d2N_dξ2, N)[0]
        assert B.shape == (6, 12), f"B shape {B.shape} != (6, 12)"
        assert element.curvature == 0.5
        assert np.any(B[0, [1, 7]] != 0), "Axial row should have curvature coupling for u_y (DOFs 1,7)"
        assert np.any(B[3, [0, 6]] != 0), "Shear row should have curvature coupling for u_x (DOFs 0,6)"
    finally:
        import shutil
        if os.path.isdir(temp_dir):
            shutil.rmtree(temp_dir, ignore_errors=True)


def test_curved_beam_B_matrix_straight_limit():
    """When κ0=0, curved B matches straight Timoshenko (axial row no u_y coupling, shear row no u_x coupling)."""
    L = 1.0
    element_dictionary = {
        "ids": np.array([0]),
        "connectivity": np.array([[0, 1]]),
        "types": np.array(["LinearCurvedTimoshenkoBeamElement3D"]),
        "integration_orders": {
            "axial": np.array([2]),
            "bending_y": np.array([2]),
            "bending_z": np.array([2]),
            "shear_y": np.array([2]),
            "shear_z": np.array([2]),
            "torsion": np.array([2]),
            "load": np.array([2]),
        },
        "curvature": np.array([0.0]),
    }
    grid_dictionary = {"coordinates": np.array([[0.0, 0.0, 0.0], [L, 0.0, 0.0]])}
    material_dictionary = {
        "E": np.array([2.1e11]),
        "G": np.array([8.1e10]),
        "nu": np.array([0.3]),
        "rho": np.array([7850.0]),
    }
    section_dictionary = {
        "A": np.array([0.001]),
        "I_x": np.array([0.0]),
        "I_y": np.array([1e-8]),
        "I_z": np.array([1e-8]),
        "J_t": np.array([1e-9]),
        "kappa": np.array([5.0 / 6.0]),
    }
    point_load_array = np.empty((0, 9))
    distributed_load_array = np.empty((0, 9))

    import tempfile
    temp_dir = tempfile.mkdtemp()
    job_results_dir = os.path.join(temp_dir, "curved_straight")
    os.makedirs(job_results_dir, exist_ok=True)
    os.makedirs(os.path.join(job_results_dir, "element_stiffness_matrices"), exist_ok=True)
    os.makedirs(os.path.join(job_results_dir, "element_force_vectors"), exist_ok=True)

    try:
        element = LinearCurvedTimoshenkoBeamElement3D(
            element_id=0,
            element_dictionary=element_dictionary,
            grid_dictionary=grid_dictionary,
            material_dictionary=material_dictionary,
            section_dictionary=section_dictionary,
            point_load_array=point_load_array,
            distributed_load_array=distributed_load_array,
            job_results_dir=job_results_dir,
        )
        from pre_processing.element_library.linear.timoshenko.utilities.B_matrix import StrainDisplacementOperator
        straight_op = StrainDisplacementOperator(element_length=L)
        xi = np.array([0.0])
        N, dN_dξ, d2N_dξ2 = element.shape_function_operator.natural_coordinate_form(xi)
        B_curved = element.strain_displacement_operator.physical_coordinate_form(dN_dξ, d2N_dξ2, N)[0]
        B_straight = straight_op.physical_coordinate_form(dN_dξ, d2N_dξ2, N)[0]
        np.testing.assert_allclose(B_curved, B_straight, atol=1e-12, err_msg="κ0=0 curved B should equal straight B")
    finally:
        import shutil
        if os.path.isdir(temp_dir):
            shutil.rmtree(temp_dir, ignore_errors=True)

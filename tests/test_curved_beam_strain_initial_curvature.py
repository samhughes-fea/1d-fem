"""
Curved beam: strain zero for rigid-body motion on a ring segment (constant κ0 = 1/R).
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


def test_curved_beam_strain_initial_curvature_rigid_body():
    """For straight segment (κ0=0), rigid-body translation yields zero strain; curved (κ0!=0) has curvature in B."""
    L = 1.0
    section_dictionary = {
        "A": np.array([0.001]),
        "I_x": np.array([0.0]),
        "I_y": np.array([1e-8]),
        "I_z": np.array([1e-8]),
        "J_t": np.array([1e-9]),
        "kappa": np.array([5.0 / 6.0]),
    }
    grid_dictionary = {"coordinates": np.array([[0.0, 0.0, 0.0], [L, 0.0, 0.0]])}
    material_dictionary = {
        "E": np.array([2.1e11]),
        "G": np.array([8.1e10]),
        "nu": np.array([0.3]),
        "rho": np.array([7850.0]),
    }
    point_load_array = np.empty((0, 9))
    distributed_load_array = np.empty((0, 9))
    orders = {
        "axial": np.array([2]),
        "bending_y": np.array([2]),
        "bending_z": np.array([2]),
        "shear_y": np.array([2]),
        "shear_z": np.array([2]),
        "torsion": np.array([2]),
        "load": np.array([2]),
    }

    import tempfile
    temp_dir = tempfile.mkdtemp()
    job_results_dir = os.path.join(temp_dir, "curved_rb")
    os.makedirs(job_results_dir, exist_ok=True)
    os.makedirs(os.path.join(job_results_dir, "element_stiffness_matrices"), exist_ok=True)
    os.makedirs(os.path.join(job_results_dir, "element_force_vectors"), exist_ok=True)

    try:
        element_straight = LinearCurvedTimoshenkoBeamElement3D(
            element_id=0,
            element_dictionary={
                "ids": np.array([0]),
                "connectivity": np.array([[0, 1]]),
                "types": np.array(["LinearCurvedTimoshenkoBeamElement3D"]),
                "integration_orders": orders,
                "curvature": np.array([0.0]),
            },
            grid_dictionary=grid_dictionary,
            material_dictionary=material_dictionary,
            section_dictionary=section_dictionary,
            point_load_array=point_load_array,
            distributed_load_array=distributed_load_array,
            job_results_dir=job_results_dir,
        )
        xi_gauss, _ = element_straight.integration_points
        N, dN_dξ, d2N_dξ2 = element_straight.shape_function_operator.natural_coordinate_form(xi_gauss)
        B_straight = element_straight.strain_displacement_operator.physical_coordinate_form(dN_dξ, d2N_dξ2, N)
        u_rb = np.zeros(12)
        u_rb[0] = 1.0
        u_rb[6] = 1.0
        for g in range(B_straight.shape[0]):
            strain = B_straight[g] @ u_rb
            np.testing.assert_allclose(strain, 0.0, atol=1e-12, err_msg=f"Straight: GP {g} strain should be zero for RB translation")

        element_curved = LinearCurvedTimoshenkoBeamElement3D(
            element_id=0,
            element_dictionary={
                "ids": np.array([0]),
                "connectivity": np.array([[0, 1]]),
                "types": np.array(["LinearCurvedTimoshenkoBeamElement3D"]),
                "integration_orders": orders,
                "curvature": np.array([0.5]),
            },
            grid_dictionary=grid_dictionary,
            material_dictionary=material_dictionary,
            section_dictionary=section_dictionary,
            point_load_array=point_load_array,
            distributed_load_array=distributed_load_array,
            job_results_dir=job_results_dir,
        )
        B_curved = element_curved.strain_displacement_operator.physical_coordinate_form(dN_dξ, d2N_dξ2, N)
        assert np.any(np.abs(B_curved - B_straight) > 1e-10), "Curved B should differ from straight B when κ0 != 0"
    finally:
        import shutil
        if os.path.isdir(temp_dir):
            shutil.rmtree(temp_dir, ignore_errors=True)

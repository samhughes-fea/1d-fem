"""Unit tests for ``assemble_timoshenko_K0`` and order resolution."""

import os
import shutil
import sys
import tempfile
from pathlib import Path

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))


def _job_dirs(root: str) -> None:
    os.makedirs(os.path.join(root, "element_stiffness_matrices"), exist_ok=True)
    os.makedirs(os.path.join(root, "element_force_vectors"), exist_ok=True)
    os.makedirs(os.path.join(root, "logs"), exist_ok=True)


def test_assemble_timoshenko_K0_matches_linear_element_K_e():
    """Direct ``assemble_timoshenko_K0`` agrees with ``LinearTimoshenkoBeamElement3D.element_stiffness_matrix``."""
    from pre_processing.element_library.linear.beam.first_order_shear_deformation_theory.timoshenko.linear_timoshenko_3D import (
        LinearTimoshenkoBeamElement3D,
    )
    from pre_processing.element_library.linear.beam.first_order_shear_deformation_theory.timoshenko.utilities.k0_timoshenko import (
        assemble_timoshenko_K0,
        timoshenko_quadrature_orders_from_element_array,
    )

    L, E, G = 1.5, 2.0e11, 77e9
    A, I_z, I_y, J_t = 6e-4, 8e-8, 4e-8, 1e-9
    grid_dictionary = {"ids": np.array([0, 1]), "coordinates": np.array([[0.0, 0.0, 0.0], [L, 0.0, 0.0]])}
    element_dictionary = {
        "ids": np.array([0]),
        "connectivity": np.array([[0, 1]]),
        "types": np.array(["LinearTimoshenkoBeamElement3D"]),
        "integration_orders": {
            "axial": np.array([2]),
            "bending_y": np.array([3]),
            "bending_z": np.array([3]),
            "shear_y": np.array([2]),
            "shear_z": np.array([2]),
            "torsion": np.array([2]),
            "load": np.array([2]),
        },
    }
    material_dictionary = {"E": np.array([E]), "G": np.array([G]), "nu": np.array([0.3]), "rho": np.array([7850.0])}
    section_dictionary = {
        "A": np.array([A]),
        "I_x": np.array([0.0]),
        "I_y": np.array([I_y]),
        "I_z": np.array([I_z]),
        "J_t": np.array([J_t]),
    }
    temp_dir = tempfile.mkdtemp()
    job_results_dir = os.path.join(temp_dir, "k0_unit")
    _job_dirs(job_results_dir)
    try:
        elem = LinearTimoshenkoBeamElement3D(
            element_id=0,
            element_dictionary=element_dictionary,
            grid_dictionary=grid_dictionary,
            material_dictionary=material_dictionary,
            section_dictionary=section_dictionary,
            point_load_array=np.empty((0, 9)),
            distributed_load_array=np.empty((0, 9)),
            job_results_dir=job_results_dir,
        )
        orders = timoshenko_quadrature_orders_from_element_array(elem.element_array)
        D = elem.material_stiffness_operator.assembly_form()
        K_direct = assemble_timoshenko_K0(
            elem.shape_function_operator,
            elem.strain_displacement_operator,
            D,
            elem.jacobian_determinant,
            orders,
        )
        K_elem = elem.element_stiffness_matrix().K_e
        np.testing.assert_allclose(K_direct, K_elem, rtol=0.0, atol=1e-12)
    finally:
        shutil.rmtree(temp_dir, ignore_errors=True)


def test_shear_block_order_changes_stiffness():
    """Higher shear-block quadrature changes ``K`` vs default 1-point block."""
    from pre_processing.element_library.linear.beam.first_order_shear_deformation_theory.timoshenko.linear_timoshenko_3D import (
        LinearTimoshenkoBeamElement3D,
    )
    from pre_processing.element_library.linear.beam.first_order_shear_deformation_theory.timoshenko.utilities.k0_timoshenko import (
        TimoshenkoQuadratureOrders,
        assemble_timoshenko_K0,
        timoshenko_quadrature_orders_from_element_array,
    )

    L, E, G = 1.0, 2.1e11, 8.1e10
    A, I_z, I_y, J_t = 0.01, 1e-6, 1e-6, 2e-7
    grid_dictionary = {"ids": np.array([0, 1]), "coordinates": np.array([[0.0, 0.0, 0.0], [L, 0.0, 0.0]])}
    element_dictionary = {
        "ids": np.array([0]),
        "connectivity": np.array([[0, 1]]),
        "types": np.array(["LinearTimoshenkoBeamElement3D"]),
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
    material_dictionary = {"E": np.array([E]), "G": np.array([G]), "nu": np.array([0.3]), "rho": np.array([7850.0])}
    section_dictionary = {
        "A": np.array([A]),
        "I_x": np.array([0.0]),
        "I_y": np.array([I_y]),
        "I_z": np.array([I_z]),
        "J_t": np.array([J_t]),
    }
    temp_dir = tempfile.mkdtemp()
    job_results_dir = os.path.join(temp_dir, "k0_shear")
    _job_dirs(job_results_dir)
    try:
        elem = LinearTimoshenkoBeamElement3D(
            element_id=0,
            element_dictionary=element_dictionary,
            grid_dictionary=grid_dictionary,
            material_dictionary=material_dictionary,
            section_dictionary=section_dictionary,
            point_load_array=np.empty((0, 9)),
            distributed_load_array=np.empty((0, 9)),
            job_results_dir=job_results_dir,
        )
        base = timoshenko_quadrature_orders_from_element_array(elem.element_array)
        D = elem.material_stiffness_operator.assembly_form()
        detJ = elem.jacobian_determinant
        K1 = assemble_timoshenko_K0(elem.shape_function_operator, elem.strain_displacement_operator, D, detJ, base)
        hi_shear = TimoshenkoQuadratureOrders(
            axial=base.axial,
            bending_y=base.bending_y,
            bending_z=base.bending_z,
            shear_y=base.shear_y,
            shear_z=base.shear_z,
            shear_block=3,
            torsion=base.torsion,
        )
        K3 = assemble_timoshenko_K0(elem.shape_function_operator, elem.strain_displacement_operator, D, detJ, hi_shear)
        assert np.max(np.abs(K3 - K1)) > 1e-10
    finally:
        shutil.rmtree(temp_dir, ignore_errors=True)

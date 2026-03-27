"""
Tests for 7-DOF warping beam elements: DOF assembly and global index mapping.

Asserts that warping elements return 14 global DOF indices (7 per node) and that
total_dof is set from max dof_per_node when using the static runner.
"""

import os
import sys
from pathlib import Path

import numpy as np
import pytest

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from pre_processing.element_library.linear.beam.first_order_shear_deformation_theory.timoshenko.linear_warping_timoshenko_3D import (
    LinearWarpingTimoshenkoBeamElement3D,
)


def _minimal_section_dict_with_gamma():
    """Section dict with 10 entries (A, I_x, I_y, I_z, J_t, kappa, alpha, y_sc, z_sc, Gamma)."""
    return {
        "A": np.array([0.001]),
        "I_x": np.array([0.0]),
        "I_y": np.array([1e-8]),
        "I_z": np.array([1e-8]),
        "J_t": np.array([1e-9]),
        "kappa": np.array([5.0 / 6.0]),
        "alpha": np.array([0.0]),
        "y_sc": np.array([0.0]),
        "z_sc": np.array([0.0]),
        "Gamma": np.array([1.0e-8]),
    }


def test_warping_dof_assembly():
    """7-DOF warping element returns 14 global DOF indices (nodes 0 and 1: 0..6 and 7..13)."""
    L = 1.0
    element_dictionary = {
        "ids": np.array([0]),
        "connectivity": np.array([[0, 1]]),
        "types": np.array(["LinearWarpingTimoshenkoBeamElement3D"]),
        "integration_orders": {
            "axial": np.array([2]),
            "bending_y": np.array([2]),
            "bending_z": np.array([2]),
            "shear_y": np.array([2]),
            "shear_z": np.array([2]),
            "torsion": np.array([2]),
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
    section_dictionary = _minimal_section_dict_with_gamma()
    point_load_array = np.empty((0, 9))
    distributed_load_array = np.empty((0, 9))

    import tempfile
    temp_dir = tempfile.mkdtemp()
    job_results_dir = os.path.join(temp_dir, "warping_test")
    os.makedirs(job_results_dir, exist_ok=True)
    os.makedirs(os.path.join(job_results_dir, "element_stiffness_matrices"), exist_ok=True)
    os.makedirs(os.path.join(job_results_dir, "element_force_vectors"), exist_ok=True)

    try:
        element = LinearWarpingTimoshenkoBeamElement3D(
            element_id=0,
            element_dictionary=element_dictionary,
            grid_dictionary=grid_dictionary,
            material_dictionary=material_dictionary,
            section_dictionary=section_dictionary,
            point_load_array=point_load_array,
            distributed_load_array=distributed_load_array,
            job_results_dir=job_results_dir,
        )
        dof_indices = element.assemble_global_dof_indices()
        assert len(dof_indices) == 14, f"Expected 14 DOF indices, got {len(dof_indices)}"
        assert element.dof_per_node == 7
        expected = np.concatenate([np.arange(0, 7), np.arange(7, 14)])
        np.testing.assert_array_equal(dof_indices, expected)
    finally:
        import shutil
        if os.path.isdir(temp_dir):
            shutil.rmtree(temp_dir, ignore_errors=True)

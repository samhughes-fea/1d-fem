"""
Tests that section pipeline supplies Γ and warping element receives it via section_array.

General section integration (or parser with Gamma column) should provide Gamma;
element reads section_array[9] when section_array.size >= 10.
"""

import os
import sys
from pathlib import Path

import numpy as np
import pytest

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from pre_processing.element_library.linear.timoshenko.linear_warping_timoshenko_3D import (
    LinearWarpingTimoshenkoBeamElement3D,
)


def test_section_integration_Gamma_element_receives():
    """Warping element reads Γ from section_array when section_dictionary has Gamma."""
    L = 1.0
    Gamma_value = 2.5e-9
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
    section_dictionary = {
        "A": np.array([0.001]),
        "I_x": np.array([0.0]),
        "I_y": np.array([1e-8]),
        "I_z": np.array([1e-8]),
        "J_t": np.array([1e-9]),
        "kappa": np.array([5.0 / 6.0]),
        "alpha": np.array([0.0]),
        "y_sc": np.array([0.0]),
        "z_sc": np.array([0.0]),
        "Gamma": np.array([Gamma_value]),
    }
    point_load_array = np.empty((0, 9))
    distributed_load_array = np.empty((0, 9))

    import tempfile
    temp_dir = tempfile.mkdtemp()
    job_results_dir = os.path.join(temp_dir, "gamma_test")
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
        assert element.section_array.size >= 10, "section_array must include Gamma (index 9)"
        assert element.Gamma == Gamma_value, f"Element Gamma {element.Gamma} != {Gamma_value}"
        # D[6,6] = E*Gamma
        D = element._D_matrix_7x7()
        assert D[6, 6] == element.E * Gamma_value
    finally:
        import shutil
        if os.path.isdir(temp_dir):
            shutil.rmtree(temp_dir, ignore_errors=True)


def test_section_integration_Gamma_optional_zero():
    """When Gamma is not in section_dictionary, element uses Γ = 0 (no warping stiffness)."""
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
    # No Gamma key: section_array has 9 entries (y_sc, z_sc only)
    section_dictionary = {
        "A": np.array([0.001]),
        "I_x": np.array([0.0]),
        "I_y": np.array([1e-8]),
        "I_z": np.array([1e-8]),
        "J_t": np.array([1e-9]),
        "kappa": np.array([5.0 / 6.0]),
        "alpha": np.array([0.0]),
        "y_sc": np.array([0.0]),
        "z_sc": np.array([0.0]),
    }
    point_load_array = np.empty((0, 9))
    distributed_load_array = np.empty((0, 9))

    import tempfile
    temp_dir = tempfile.mkdtemp()
    job_results_dir = os.path.join(temp_dir, "gamma_zero_test")
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
        # section_array has 9 entries (no Gamma), so element.Gamma defaults to 0
        assert element.Gamma == 0.0, "When Gamma not in section, element should use 0"
        D = element._D_matrix_7x7()
        assert D[6, 6] == 0.0, "D[6,6] should be 0 when Γ = 0"
    finally:
        import shutil
        if os.path.isdir(temp_dir):
            shutil.rmtree(temp_dir, ignore_errors=True)

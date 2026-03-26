"""
Unit tests for Truss and Bar 3D elements (axial, transverse, torsion).
"""

import os
import sys
import tempfile
from pathlib import Path

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from pre_processing.element_library.linear.truss.linear_truss_3D import LinearTrussElement3D
from pre_processing.element_library.linear.bar.linear_bar_3D import LinearBarElement3D


def _make_job_dir():
    d = tempfile.mkdtemp()
    job = os.path.join(d, "test_results")
    os.makedirs(job, exist_ok=True)
    os.makedirs(os.path.join(job, "element_stiffness_matrices"), exist_ok=True)
    os.makedirs(os.path.join(job, "element_force_vectors"), exist_ok=True)
    return job


def test_truss_element_stiffness_shape_and_symmetry():
    """Truss K_e is 12×12 and symmetric; has axial, transverse, torsion contributions."""
    E, G = 2.1e11, 8.1e10
    A, J_t = 0.001, 1e-8
    L = 1.0
    node_coords = np.array([[0.0, 0.0, 0.0], [L, 0.0, 0.0]], dtype=np.float64)

    element_dictionary = {
        "ids": np.array([0]),
        "connectivity": np.array([[0, 1]]),
        "types": np.array(["TrussElement3D"]),
        "integration_orders": {
            "axial": np.array([1]),
            "bending_y": np.array([0]),
            "bending_z": np.array([0]),
            "shear_y": np.array([0]),
            "shear_z": np.array([0]),
            "torsion": np.array([1]),
            "load": np.array([1]),
        },
    }
    grid_dictionary = {"coordinates": node_coords}
    section_dictionary = {
        "A": np.array([A]),
        "I_x": np.array([0.0]),
        "I_y": np.array([0.0]),
        "I_z": np.array([0.0]),
        "J_t": np.array([J_t]),
    }
    material_dictionary = {
        "E": np.array([E]),
        "G": np.array([G]),
        "nu": np.array([0.3]),
        "rho": np.array([7850.0]),
    }
    job_results_dir = _make_job_dir()

    el = LinearTrussElement3D(
        element_id=0,
        element_dictionary=element_dictionary,
        grid_dictionary=grid_dictionary,
        section_dictionary=section_dictionary,
        material_dictionary=material_dictionary,
        point_load_array=np.empty((0, 9)),
        distributed_load_array=np.empty((0, 9)),
        job_results_dir=job_results_dir,
    )
    obj = el.element_stiffness_matrix()
    K_e = obj.K_e

    assert K_e.shape == (12, 12), "Truss K_e must be 12×12"
    np.testing.assert_allclose(K_e, K_e.T, err_msg="K_e must be symmetric")
    # Axial along x: DOF 0 and 6 coupled; torsion DOF 3 and 9
    assert K_e[0, 0] > 0 and K_e[6, 6] > 0, "Axial diagonal entries positive"
    assert K_e[3, 3] > 0 and K_e[9, 9] > 0, "Torsion diagonal entries positive"
    # Expected axial stiffness contribution at (0,0): (EA/L) * cx^2 = EA/L for aligned element
    expected_axial = E * A / L
    np.testing.assert_allclose(K_e[0, 0], expected_axial, rtol=1e-5)


def test_bar_element_stiffness_shape_and_blocks():
    """Bar K_e is 12×12; axial and torsion blocks only (no transverse)."""
    E, G = 2.1e11, 8.1e10
    A, J_t = 0.001, 1e-8
    L = 1.0
    node_coords = np.array([[0.0, 0.0, 0.0], [L, 0.0, 0.0]], dtype=np.float64)

    element_dictionary = {
        "ids": np.array([0]),
        "connectivity": np.array([[0, 1]]),
        "types": np.array(["LinearBarElement3D"]),
        "integration_orders": {
            "axial": np.array([1]),
            "bending_y": np.array([0]),
            "bending_z": np.array([0]),
            "shear_y": np.array([0]),
            "shear_z": np.array([0]),
            "torsion": np.array([1]),
            "load": np.array([1]),
        },
    }
    grid_dictionary = {"coordinates": node_coords}
    section_dictionary = {
        "A": np.array([A]),
        "I_x": np.array([0.0]),
        "I_y": np.array([0.0]),
        "I_z": np.array([0.0]),
        "J_t": np.array([J_t]),
    }
    material_dictionary = {
        "E": np.array([E]),
        "G": np.array([G]),
        "nu": np.array([0.3]),
        "rho": np.array([7850.0]),
    }
    job_results_dir = _make_job_dir()

    el = LinearBarElement3D(
        element_id=0,
        element_dictionary=element_dictionary,
        grid_dictionary=grid_dictionary,
        section_dictionary=section_dictionary,
        material_dictionary=material_dictionary,
        point_load_array=np.empty((0, 9)),
        distributed_load_array=np.empty((0, 9)),
        job_results_dir=job_results_dir,
    )
    obj = el.element_stiffness_matrix()
    K_e = obj.K_e

    assert K_e.shape == (12, 12), "Bar K_e must be 12×12"
    np.testing.assert_allclose(K_e, K_e.T, err_msg="K_e must be symmetric")
    expected_axial = E * A / L
    expected_torsion = G * J_t / L
    np.testing.assert_allclose(K_e[0, 0], expected_axial, rtol=1e-5)
    np.testing.assert_allclose(K_e[3, 3], expected_torsion, rtol=1e-5)
    np.testing.assert_allclose(K_e[9, 9], expected_torsion, rtol=1e-5)


def test_truss_force_vector():
    """Truss element_force_vector returns F_e shape (12,)."""
    L = 1.0
    node_coords = np.array([[0.0, 0.0, 0.0], [L, 0.0, 0.0]], dtype=np.float64)
    element_dictionary = {
        "ids": np.array([0]),
        "connectivity": np.array([[0, 1]]),
        "types": np.array(["TrussElement3D"]),
        "integration_orders": {k: np.array([1]) for k in ["axial", "bending_y", "bending_z", "shear_y", "shear_z", "torsion", "load"]},
    }
    grid_dictionary = {"coordinates": node_coords}
    section_dictionary = {"A": np.array([0.001]), "I_x": np.array([0.0]), "I_y": np.array([0.0]), "I_z": np.array([0.0]), "J_t": np.array([1e-8])}
    material_dictionary = {"E": np.array([2.1e11]), "G": np.array([8.1e10]), "nu": np.array([0.3]), "rho": np.array([7850.0])}
    job_results_dir = _make_job_dir()

    el = LinearTrussElement3D(
        element_id=0,
        element_dictionary=element_dictionary,
        grid_dictionary=grid_dictionary,
        section_dictionary=section_dictionary,
        material_dictionary=material_dictionary,
        point_load_array=np.empty((0, 9)),
        distributed_load_array=np.empty((0, 9)),
        job_results_dir=job_results_dir,
    )
    force_obj = el.element_force_vector()
    assert force_obj.F_e.shape == (12,)
    assert force_obj.element_type == "Truss-3D"


def test_bar_force_vector():
    """Bar element_force_vector returns F_e shape (12,)."""
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
    job_results_dir = _make_job_dir()

    el = LinearBarElement3D(
        element_id=0,
        element_dictionary=element_dictionary,
        grid_dictionary=grid_dictionary,
        section_dictionary=section_dictionary,
        material_dictionary=material_dictionary,
        point_load_array=np.empty((0, 9)),
        distributed_load_array=np.empty((0, 9)),
        job_results_dir=job_results_dir,
    )
    force_obj = el.element_force_vector()
    assert force_obj.F_e.shape == (12,)
    assert force_obj.element_type == "Bar-3D"

"""
Tests for warping beam element stiffness: singular/restraint behaviour.

Cantilever with warping free at tip: K_e has rigid-body modes; with warping restrained
at root we get a non-singular block. Asserts shape (14,14) and basic symmetry.
"""

import os
import sys
from pathlib import Path

import numpy as np
import pytest

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from pre_processing.element_library.linear.beam.first_order_shear_deformation_theory.timoshenko.linear_timoshenko_3D import (
    LinearTimoshenkoBeamElement3D,
)


def _minimal_section_dict_with_gamma(gamma=1.0e-8):
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
        "Gamma": np.array([gamma]),
    }


def test_warping_stiffness_matrix_restraint():
    """Warping Timoshenko K_e is (14,14), symmetric; with Γ > 0 the warping block contributes."""
    L = 1.0
    element_dictionary = {
        "ids": np.array([0]),
        "connectivity": np.array([[0, 1]]),
        "types": np.array(["LinearTimoshenkoBeamElement3D"]),
        "warping": np.array([1], dtype=np.int8),
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
    section_dictionary = _minimal_section_dict_with_gamma(gamma=1.0e-8)
    point_load_array = np.empty((0, 9))
    distributed_load_array = np.empty((0, 9))

    import tempfile
    temp_dir = tempfile.mkdtemp()
    job_results_dir = os.path.join(temp_dir, "warping_stiff_test")
    os.makedirs(job_results_dir, exist_ok=True)
    os.makedirs(os.path.join(job_results_dir, "element_stiffness_matrices"), exist_ok=True)
    os.makedirs(os.path.join(job_results_dir, "element_force_vectors"), exist_ok=True)

    try:
        element = LinearTimoshenkoBeamElement3D(
            element_id=0,
            element_dictionary=element_dictionary,
            grid_dictionary=grid_dictionary,
            material_dictionary=material_dictionary,
            section_dictionary=section_dictionary,
            point_load_array=point_load_array,
            distributed_load_array=distributed_load_array,
            job_results_dir=job_results_dir,
        )
        obj = element.element_stiffness_matrix()
        Ke = obj.K_e
        assert Ke.shape == (14, 14), f"K_e shape {Ke.shape} != (14, 14)"
        np.testing.assert_allclose(Ke, Ke.T, rtol=0, atol=1e-12, err_msg="K_e must be symmetric")
        # Warping DOFs 12, 13: block [12:14, 12:14] should have non-zero diagonal when Γ > 0
        K_warp = Ke[12:14, 12:14]
        assert np.linalg.norm(K_warp) > 1e-20, "Warping block should contribute when Γ > 0"
        # Element free (no BCs): singular (at least 6 rigid-body modes for 3D beam)
        ev = np.linalg.eigvalsh(Ke)
        n_zero = np.sum(np.abs(ev) < 1e-6)
        assert n_zero >= 6, f"Expect at least 6 near-zero eigenvalues for free element, got {n_zero}"
    finally:
        import shutil
        if os.path.isdir(temp_dir):
            shutil.rmtree(temp_dir, ignore_errors=True)

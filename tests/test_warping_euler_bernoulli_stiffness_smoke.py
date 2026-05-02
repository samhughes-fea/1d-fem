"""
Smoke test: warping Euler–Bernoulli element stiffness shape and symmetry.

Uses ``LinearEulerBernoulliBeamElement3D`` with explicit ``warping`` in the element dictionary
(replaces the removed ``LinearWarpingEulerBernoulliBeamElement3D`` alias).
"""

import os
import sys
import tempfile
from pathlib import Path

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from pre_processing.element_library.linear.beam.zero_order_shear_deformation_theory.euler_bernoulli.linear_euler_bernoulli_3D import (
    LinearEulerBernoulliBeamElement3D,
)


def _minimal_section_with_gamma(gamma: float = 1.0e-8):
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


def test_warping_euler_bernoulli_stiffness_smoke():
    """K_e is (14,14), symmetric; warping block non-zero when Γ > 0."""
    L = 1.0
    element_dictionary = {
        "ids": np.array([0]),
        "connectivity": np.array([[0, 1]]),
        "types": np.array(["LinearEulerBernoulliBeamElement3D"]),
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
    section_dictionary = _minimal_section_with_gamma(gamma=1.0e-8)
    point_load_array = np.empty((0, 9))
    distributed_load_array = np.empty((0, 9))

    temp_dir = tempfile.mkdtemp()
    job_results_dir = os.path.join(temp_dir, "warping_eb_smoke")
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
            point_load_array=point_load_array,
            distributed_load_array=distributed_load_array,
            job_results_dir=job_results_dir,
        )
        obj = element.element_stiffness_matrix()
        Ke = obj.K_e
        assert Ke.shape == (14, 14)
        np.testing.assert_allclose(Ke, Ke.T, rtol=1e-12, atol=1e-8)
        assert np.linalg.norm(Ke[12:14, 12:14]) > 1e-20
    finally:
        import shutil

        if os.path.isdir(temp_dir):
            shutil.rmtree(temp_dir, ignore_errors=True)

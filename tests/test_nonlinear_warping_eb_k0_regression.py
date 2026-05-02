"""Nonlinear warping EB: K_T(U=0) matches linear warping-EB K_e (14×14)."""

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


def _shared_dicts(L: float = 1.0):
    """Grid + element with [warping]=1; section has Γ at index 9 (full optional block)."""
    grid_dictionary = {
        "ids": np.array([0, 1]),
        "coordinates": np.array([[0.0, 0.0, 0.0], [L, 0.0, 0.0]]),
    }
    element_dictionary = {
        "ids": np.array([0]),
        "connectivity": np.array([[0, 1]]),
        "types": np.array(["LinearEulerBernoulliBeamElement3D"]),
        "warping": np.array([1], dtype=np.int8),
        "integration_orders": {
            "axial": np.array([3]),
            "bending_y": np.array([3]),
            "bending_z": np.array([3]),
            "shear_y": np.array([0]),
            "shear_z": np.array([0]),
            "torsion": np.array([3]),
            "load": np.array([2]),
        },
    }
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
        "kappa": np.array([1.0]),
        "alpha": np.array([0.0]),
        "y_sc": np.array([0.0]),
        "z_sc": np.array([0.0]),
        "Gamma": np.array([1.0e-8]),
    }
    return grid_dictionary, element_dictionary, material_dictionary, section_dictionary


def test_nonlinear_warping_eb_K0_matches_linear_warping_eb():
    from pre_processing.element_library.linear.beam.zero_order_shear_deformation_theory.euler_bernoulli.linear_euler_bernoulli_3D import (
        LinearEulerBernoulliBeamElement3D,
    )
    from pre_processing.element_library.nonlinear.euler_bernoulli.nonlinear_euler_bernoulli_3D import (
        NonlinearEulerBernoulliBeamElement3D,
    )

    tmp = tempfile.mkdtemp()
    job = os.path.join(tmp, "job")
    _job_dirs(job)
    try:
        g, e, m, s = _shared_dicts()
        e_nl = {**e, "types": np.array(["NonlinearEulerBernoulliBeamElement3D"])}
        pl = np.empty((0, 9))
        dl = np.empty((0, 9))

        lin = LinearEulerBernoulliBeamElement3D(
            element_id=0,
            element_dictionary=e,
            grid_dictionary=g,
            material_dictionary=m,
            section_dictionary=s,
            point_load_array=pl,
            distributed_load_array=dl,
            job_results_dir=job,
        )
        nlin = NonlinearEulerBernoulliBeamElement3D(
            element_id=0,
            element_dictionary=e_nl,
            grid_dictionary=g,
            material_dictionary=m,
            section_dictionary=s,
            point_load_array=pl,
            distributed_load_array=dl,
            job_results_dir=job,
        )
        K_lin = lin.element_stiffness_matrix().K_e
        U0 = np.zeros(14, dtype=np.float64)
        K_n = nlin.tangent_stiffness_matrix(U0)
        np.testing.assert_allclose(K_n, K_lin, rtol=1e-10, atol=1e-8)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

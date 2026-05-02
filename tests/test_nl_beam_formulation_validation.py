"""Roadmap validation: tangent symmetry and equilibrium sanity for NL EB / Timoshenko beams."""

import os
import shutil
import sys
import tempfile
from pathlib import Path

import numpy as np
import pytest

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))


def _job_dirs(root: str) -> None:
    os.makedirs(os.path.join(root, "element_stiffness_matrices"), exist_ok=True)
    os.makedirs(os.path.join(root, "element_force_vectors"), exist_ok=True)
    os.makedirs(os.path.join(root, "logs"), exist_ok=True)


def _cantilever_dicts(L=1.0, E=2.1e11, G=8.1e10):
    grid_dictionary = {
        "ids": np.array([0, 1]),
        "coordinates": np.array([[0.0, 0.0, 0.0], [L, 0.0, 0.0]]),
    }
    element_dictionary = {
        "ids": np.array([0]),
        "connectivity": np.array([[0, 1]]),
        "types": np.array(["NonlinearTimoshenkoBeamElement3D"]),
        "integration_orders": {
            "axial": np.array([3]),
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
        "A": np.array([0.001]),
        "I_x": np.array([0.0]),
        "I_y": np.array([1e-8]),
        "I_z": np.array([1e-8]),
        "J_t": np.array([1e-9]),
    }
    return grid_dictionary, element_dictionary, material_dictionary, section_dictionary


@pytest.mark.parametrize("seed", [0, 1, 42])
def test_nonlinear_timoshenko_tangent_symmetric_random_displacement(seed: int) -> None:
    """``K_T`` is numerically symmetric at random ``U_e``."""
    from pre_processing.element_library.nonlinear.timoshenko.nonlinear_timoshenko_3D import (
        NonlinearTimoshenkoBeamElement3D,
    )

    grid_dictionary, element_dictionary, material_dictionary, section_dictionary = _cantilever_dicts()
    point_load_array = np.empty((0, 9))
    distributed_load_array = np.empty((0, 9))

    temp_dir = tempfile.mkdtemp()
    job_results_dir = os.path.join(temp_dir, "nl_val")
    _job_dirs(job_results_dir)
    try:
        elem = NonlinearTimoshenkoBeamElement3D(
            element_id=0,
            element_dictionary=element_dictionary,
            grid_dictionary=grid_dictionary,
            material_dictionary=material_dictionary,
            section_dictionary=section_dictionary,
            point_load_array=point_load_array,
            distributed_load_array=distributed_load_array,
            job_results_dir=job_results_dir,
        )
        rng = np.random.default_rng(seed)
        U = rng.standard_normal(12) * 0.01
        K = elem.tangent_stiffness_matrix(U)
        np.testing.assert_allclose(K, K.T, rtol=1e-12, atol=1e-12)
    finally:
        shutil.rmtree(temp_dir, ignore_errors=True)


def test_nonlinear_euler_bernoulli_tangent_symmetric_random_displacement() -> None:
    """NL EB ``K_T`` symmetric at small random displacement."""
    from pre_processing.element_library.nonlinear.euler_bernoulli.nonlinear_euler_bernoulli_3D import (
        NonlinearEulerBernoulliBeamElement3D,
    )

    L = 1.0
    grid_dictionary = {
        "ids": np.array([0, 1]),
        "coordinates": np.array([[0.0, 0.0, 0.0], [L, 0.0, 0.0]]),
    }
    element_dictionary = {
        "ids": np.array([0]),
        "connectivity": np.array([[0, 1]]),
        "types": np.array(["NonlinearEulerBernoulliBeamElement3D"]),
        "integration_orders": {
            "axial": np.array([3]),
            "bending_y": np.array([3]),
            "bending_z": np.array([3]),
            "shear_y": np.array([2]),
            "shear_z": np.array([2]),
            "torsion": np.array([2]),
            "load": np.array([2]),
        },
    }
    material_dictionary = {"E": np.array([2.1e11]), "G": np.array([8.1e10]), "nu": np.array([0.3]), "rho": np.array([7850.0])}
    section_dictionary = {
        "A": np.array([0.001]),
        "I_x": np.array([0.0]),
        "I_y": np.array([1e-8]),
        "I_z": np.array([1e-8]),
        "J_t": np.array([1e-9]),
    }
    temp_dir = tempfile.mkdtemp()
    job_results_dir = os.path.join(temp_dir, "nl_eb_val")
    _job_dirs(job_results_dir)
    try:
        elem = NonlinearEulerBernoulliBeamElement3D(
            element_id=0,
            element_dictionary=element_dictionary,
            grid_dictionary=grid_dictionary,
            material_dictionary=material_dictionary,
            section_dictionary=section_dictionary,
            point_load_array=np.empty((0, 9)),
            distributed_load_array=np.empty((0, 9)),
            job_results_dir=job_results_dir,
        )
        rng = np.random.default_rng(7)
        U = rng.standard_normal(12) * 0.005
        K = elem.tangent_stiffness_matrix(U)
        np.testing.assert_allclose(K, K.T, rtol=1e-11, atol=1e-11)
    finally:
        shutil.rmtree(temp_dir, ignore_errors=True)


def test_nonlinear_euler_bernoulli_internal_force_uses_B_tot_nonzero_when_nl_active() -> None:
    """With nonzero displacement, NL EB internal force differs from ``B_lin.T @ S`` when ``B_nl`` is nonzero."""
    from pre_processing.element_library.nonlinear.euler_bernoulli.nonlinear_euler_bernoulli_3D import (
        NonlinearEulerBernoulliBeamElement3D,
    )

    L = 2.0
    grid_dictionary = {
        "ids": np.array([0, 1]),
        "coordinates": np.array([[0.0, 0.0, 0.0], [L, 0.0, 0.0]]),
    }
    element_dictionary = {
        "ids": np.array([0]),
        "connectivity": np.array([[0, 1]]),
        "types": np.array(["NonlinearEulerBernoulliBeamElement3D"]),
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
    material_dictionary = {"E": np.array([2.1e11]), "G": np.array([8.1e10]), "nu": np.array([0.3]), "rho": np.array([7850.0])}
    section_dictionary = {
        "A": np.array([0.01]),
        "I_x": np.array([0.0]),
        "I_y": np.array([1e-6]),
        "I_z": np.array([1e-6]),
        "J_t": np.array([1e-7]),
    }
    temp_dir = tempfile.mkdtemp()
    job_results_dir = os.path.join(temp_dir, "nl_eb_fint")
    _job_dirs(job_results_dir)
    try:
        elem = NonlinearEulerBernoulliBeamElement3D(
            element_id=0,
            element_dictionary=element_dictionary,
            grid_dictionary=grid_dictionary,
            material_dictionary=material_dictionary,
            section_dictionary=section_dictionary,
            point_load_array=np.empty((0, 9)),
            distributed_load_array=np.empty((0, 9)),
            job_results_dir=job_results_dir,
        )
        D = elem.material_stiffness_operator.assembly_form()
        U = np.zeros(12, dtype=np.float64)
        U[0] = 0.02
        U[1] = 0.001
        F_tot = elem.internal_force_vector(U)
        # Recompute legacy-style F_lin only at one GP for magnitude check (not identical API)
        xi, w = elem.integration_points
        detJ = elem.jacobian_determinant
        dξ_dx = 2.0 / elem.L
        d2ξ_dx2 = 4.0 / (elem.L**2)
        F_lin_style = np.zeros(12, dtype=np.float64)
        for xi_g, w_g in zip(xi, w):
            N, dN_dξ, d2N_dξ2 = elem.shape_function_operator.natural_coordinate_form(np.array([xi_g]))
            dN_dx = dN_dξ.copy() * dξ_dx
            d2N_dx2 = d2N_dξ2.copy() * d2ξ_dx2
            B_lin = elem.green_lagrange_strain_operator.linearized_strain_displacement(dN_dx[0], d2N_dx2[0])
            E_lin = elem.green_lagrange_strain_operator.strain_linear_part(dN_dx[0], d2N_dx2[0], U)
            E_nl = elem.green_lagrange_strain_operator.strain_nonlinear_part(dN_dx[0], d2N_dx2[0], U)
            E = E_lin + E_nl
            S = D @ E
            F_lin_style += (B_lin.T @ S) * w_g * detJ
        assert np.linalg.norm(F_tot - F_lin_style) > 1e-6 * max(1.0, np.linalg.norm(F_tot))
    finally:
        shutil.rmtree(temp_dir, ignore_errors=True)

"""Regression: NL TL beam material stiffness matches linear reference when warping (14 DOF) is on."""

from __future__ import annotations

import os
import shutil
import sys
import tempfile
from pathlib import Path

import numpy as np
import pytest

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))


def _job_root(parent: str) -> str:
    os.makedirs(os.path.join(parent, "element_stiffness_matrices"), exist_ok=True)
    os.makedirs(os.path.join(parent, "element_force_vectors"), exist_ok=True)
    os.makedirs(os.path.join(parent, "logs"), exist_ok=True)
    return parent


def _section_gamma(gamma: float) -> dict:
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


@pytest.fixture
def warping_mesh_common():
    L = 1.0
    grid_dictionary = {"coordinates": np.array([[0.0, 0.0, 0.0], [L, 0.0, 0.0]], dtype=np.float64)}
    integration_orders = {
        "axial": np.array([2]),
        "bending_y": np.array([2]),
        "bending_z": np.array([2]),
        "shear_y": np.array([2]),
        "shear_z": np.array([2]),
        "torsion": np.array([2]),
        "load": np.array([2]),
    }
    material_dictionary = {
        "E": np.array([2.1e11]),
        "G": np.array([8.1e10]),
        "nu": np.array([0.3]),
        "rho": np.array([7850.0]),
    }
    section_dictionary = _section_gamma(1.0e-7)
    point_load_array = np.empty((0, 9))
    distributed_load_array = np.empty((0, 9))
    return (
        grid_dictionary,
        integration_orders,
        material_dictionary,
        section_dictionary,
        point_load_array,
        distributed_load_array,
    )


def test_nonlinear_timoshenko_K0_matches_linear_timoshenko_ke_warping(warping_mesh_common):
    from pre_processing.element_library.linear.beam.first_order_shear_deformation_theory.timoshenko.linear_timoshenko_3D import (
        LinearTimoshenkoBeamElement3D,
    )
    from pre_processing.element_library.nonlinear.timoshenko.nonlinear_timoshenko_3D import (
        NonlinearTimoshenkoBeamElement3D,
    )

    grid_dictionary, integration_orders, material_dictionary, section_dictionary, pl, dl = warping_mesh_common
    tmp = tempfile.mkdtemp()
    job_results_dir = _job_root(os.path.join(tmp, "job"))
    try:
        element_dictionary_lin = {
            "ids": np.array([0]),
            "connectivity": np.array([[0, 1]]),
            "types": np.array(["LinearTimoshenkoBeamElement3D"]),
            "warping": np.array([1], dtype=np.int8),
            "integration_orders": integration_orders,
        }
        element_dictionary_nl = {
            **element_dictionary_lin,
            "types": np.array(["NonlinearTimoshenkoBeamElement3D"]),
        }
        kwargs = dict(
            element_id=0,
            grid_dictionary=grid_dictionary,
            material_dictionary=material_dictionary,
            section_dictionary=section_dictionary,
            point_load_array=pl,
            distributed_load_array=dl,
            job_results_dir=job_results_dir,
        )
        lin = LinearTimoshenkoBeamElement3D(element_dictionary=element_dictionary_lin, **kwargs)
        nl = NonlinearTimoshenkoBeamElement3D(element_dictionary=element_dictionary_nl, **kwargs)
        Ke = lin.element_stiffness_matrix().K_e
        K0 = nl._get_K_0()
        np.testing.assert_allclose(K0, Ke, rtol=1e-11, atol=1e-7)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def test_nonlinear_euler_bernoulli_tangent_at_zero_matches_linear_ke_warping(warping_mesh_common):
    from pre_processing.element_library.linear.beam.zero_order_shear_deformation_theory.euler_bernoulli.linear_euler_bernoulli_3D import (
        LinearEulerBernoulliBeamElement3D,
    )
    from pre_processing.element_library.nonlinear.euler_bernoulli.nonlinear_euler_bernoulli_3D import (
        NonlinearEulerBernoulliBeamElement3D,
    )

    grid_dictionary, integration_orders, material_dictionary, section_dictionary, pl, dl = warping_mesh_common
    tmp = tempfile.mkdtemp()
    job_results_dir = _job_root(os.path.join(tmp, "job"))
    try:
        element_dictionary_lin = {
            "ids": np.array([0]),
            "connectivity": np.array([[0, 1]]),
            "types": np.array(["LinearEulerBernoulliBeamElement3D"]),
            "warping": np.array([1], dtype=np.int8),
            "integration_orders": integration_orders,
        }
        element_dictionary_nl = {
            **element_dictionary_lin,
            "types": np.array(["NonlinearEulerBernoulliBeamElement3D"]),
        }
        kwargs = dict(
            element_id=0,
            grid_dictionary=grid_dictionary,
            material_dictionary=material_dictionary,
            section_dictionary=section_dictionary,
            point_load_array=pl,
            distributed_load_array=dl,
            job_results_dir=job_results_dir,
        )
        lin = LinearEulerBernoulliBeamElement3D(element_dictionary=element_dictionary_lin, **kwargs)
        nl = NonlinearEulerBernoulliBeamElement3D(element_dictionary=element_dictionary_nl, **kwargs)
        Ke = lin.element_stiffness_matrix().K_e
        u0 = np.zeros(14, dtype=np.float64)
        Kt = nl.tangent_stiffness_matrix(u0)
        np.testing.assert_allclose(Kt, Ke, rtol=1e-10, atol=1e-6)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def test_beam_warping_policy_matches_piecewise_helpers():
    from pre_processing.element_library.beam_warping import (
        beam_warping_policy,
        effective_warping_gamma,
        element_warping_stiffness_on,
        mesh_uses_warping_dof,
    )

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
    gamma_s = 2.5e-7
    pol = beam_warping_policy(element_dictionary, 0, str(element_dictionary["types"][0]), gamma_s)
    assert pol.mesh_allocates_chi_dof == mesh_uses_warping_dof(element_dictionary)
    assert pol.warping_stiffness_on == element_warping_stiffness_on(
        element_dictionary, 0, str(element_dictionary["types"][0])
    )
    assert pol.gamma_effective == effective_warping_gamma(gamma_s, pol.warping_stiffness_on)

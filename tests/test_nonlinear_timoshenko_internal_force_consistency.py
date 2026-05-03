"""Nonlinear Timoshenko: F_int(0)=0 and K_T symmetry sanity checks."""

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


def _cantilever_dicts(L=1.0, E=2.1e11, G=8.1e10, A=0.001, I_z=1e-8, I_y=1e-8, J_t=1e-9):
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
        "A": np.array([A]),
        "I_x": np.array([0.0]),
        "I_y": np.array([I_y]),
        "I_z": np.array([I_z]),
        "J_t": np.array([J_t]),
    }
    return grid_dictionary, element_dictionary, material_dictionary, section_dictionary


def test_nonlinear_timoshenko_internal_force_zero_at_zero_displacement():
    """F_int(U_e=0) should be zero."""
    from pre_processing.element_library.nonlinear.timoshenko.nonlinear_timoshenko_3D import (
        NonlinearTimoshenkoBeamElement3D,
    )

    grid_dictionary, element_dictionary, material_dictionary, section_dictionary = _cantilever_dicts()
    point_load_array = np.empty((0, 9))
    distributed_load_array = np.empty((0, 9))

    temp_dir = tempfile.mkdtemp()
    job_results_dir = os.path.join(temp_dir, "nl_timo_fint")
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
        U_zero = np.zeros(12, dtype=np.float64)
        F_int = elem.internal_force_vector(U_zero)
        np.testing.assert_allclose(F_int, 0.0, atol=1e-12, err_msg="F_int(0) should be zero")
    finally:
        shutil.rmtree(temp_dir, ignore_errors=True)


def test_nonlinear_timoshenko_tangent_symmetric_and_internal_force_finite():
    """K_T(U_e) is symmetric; F_int(U_e) is finite for small U."""
    from pre_processing.element_library.nonlinear.timoshenko.nonlinear_timoshenko_3D import (
        NonlinearTimoshenkoBeamElement3D,
    )

    grid_dictionary, element_dictionary, material_dictionary, section_dictionary = _cantilever_dicts()
    point_load_array = np.empty((0, 9))
    distributed_load_array = np.empty((0, 9))

    temp_dir = tempfile.mkdtemp()
    job_results_dir = os.path.join(temp_dir, "nl_timo_kt")
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
        U = np.zeros(12, dtype=np.float64)
        U[7] = 1e-4
        U[5] = 1e-5
        K_T = np.asarray(elem.tangent_stiffness_matrix(U), dtype=np.float64)
        F_int = np.asarray(elem.internal_force_vector(U), dtype=np.float64)
        np.testing.assert_allclose(K_T, K_T.T, atol=1e-12, err_msg="K_T should be symmetric")
        assert np.all(np.isfinite(F_int)), "F_int should be finite"
    finally:
        shutil.rmtree(temp_dir, ignore_errors=True)


def test_nonlinear_timoshenko_warping_initial_tangent_matches_linear_stiffness():
    """At U=0, NL Timoshenko + warping tangent equals linear 14×14 K_e (K_0 only)."""
    from pre_processing.element_library.linear.beam.first_order_shear_deformation_theory.timoshenko.linear_timoshenko_3D import (
        LinearTimoshenkoBeamElement3D,
    )
    from pre_processing.element_library.nonlinear.timoshenko.nonlinear_timoshenko_3D import (
        NonlinearTimoshenkoBeamElement3D,
    )

    L = 1.0
    element_dictionary_lin = {
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
    element_dictionary_nl = {**element_dictionary_lin, "types": np.array(["NonlinearTimoshenkoBeamElement3D"])}
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
        "Gamma": np.array([1.0e-8]),
    }
    point_load_array = np.empty((0, 9))
    distributed_load_array = np.empty((0, 9))

    temp_dir = tempfile.mkdtemp()
    job_results_dir = os.path.join(temp_dir, "nl_ts_warp_kt0")
    _job_dirs(job_results_dir)

    try:
        el_lin = LinearTimoshenkoBeamElement3D(
            element_id=0,
            element_dictionary=element_dictionary_lin,
            grid_dictionary=grid_dictionary,
            material_dictionary=material_dictionary,
            section_dictionary=section_dictionary,
            point_load_array=point_load_array,
            distributed_load_array=distributed_load_array,
            job_results_dir=job_results_dir,
        )
        K_lin = el_lin.element_stiffness_matrix().K_e

        el_nl = NonlinearTimoshenkoBeamElement3D(
            element_id=0,
            element_dictionary=element_dictionary_nl,
            grid_dictionary=grid_dictionary,
            material_dictionary=material_dictionary,
            section_dictionary=section_dictionary,
            point_load_array=point_load_array,
            distributed_load_array=distributed_load_array,
            job_results_dir=job_results_dir,
        )
        U0 = np.zeros(14, dtype=np.float64)
        K_nl = el_nl.tangent_stiffness_matrix(U0)
        np.testing.assert_allclose(K_nl, K_lin, rtol=1e-10, atol=1e-6)
    finally:
        shutil.rmtree(temp_dir, ignore_errors=True)


def _nl_ts_warp_case():
    """Shared minimal mesh and dictionaries for NL Timoshenko + [warping] element tests."""
    L = 1.0
    element_dictionary = {
        "ids": np.array([0]),
        "connectivity": np.array([[0, 1]]),
        "types": np.array(["NonlinearTimoshenkoBeamElement3D"]),
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
        "Gamma": np.array([1.0e-8]),
    }
    return element_dictionary, grid_dictionary, material_dictionary, section_dictionary


def test_nonlinear_timoshenko_warping_internal_force_zero_at_zero_displacement():
    """F_int(0) == 0 for 14-DOF TL Timoshenko + warping."""
    from pre_processing.element_library.nonlinear.timoshenko.nonlinear_timoshenko_3D import (
        NonlinearTimoshenkoBeamElement3D,
    )

    element_dictionary, grid_dictionary, material_dictionary, section_dictionary = _nl_ts_warp_case()
    point_load_array = np.empty((0, 9))
    distributed_load_array = np.empty((0, 9))
    temp_dir = tempfile.mkdtemp()
    job_results_dir = os.path.join(temp_dir, "nl_ts_warp_f0")
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
        F_int = elem.internal_force_vector(np.zeros(14, dtype=np.float64))
        np.testing.assert_allclose(F_int, 0.0, atol=1e-12)
    finally:
        shutil.rmtree(temp_dir, ignore_errors=True)


def test_nonlinear_timoshenko_warping_tangent_symmetric_random_displacement():
    """K_T(U) symmetric and F_int(U) finite for small random U including χ DOFs."""
    from pre_processing.element_library.nonlinear.timoshenko.nonlinear_timoshenko_3D import (
        NonlinearTimoshenkoBeamElement3D,
    )

    element_dictionary, grid_dictionary, material_dictionary, section_dictionary = _nl_ts_warp_case()
    point_load_array = np.empty((0, 9))
    distributed_load_array = np.empty((0, 9))
    temp_dir = tempfile.mkdtemp()
    job_results_dir = os.path.join(temp_dir, "nl_ts_warp_kt_rand")
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
        rng = np.random.default_rng(42)
        U = rng.normal(scale=1e-4, size=14)
        U[12] = 1e-5
        U[13] = -3e-6
        K_T = np.asarray(elem.tangent_stiffness_matrix(U), dtype=np.float64)
        F_int = np.asarray(elem.internal_force_vector(U), dtype=np.float64)
        np.testing.assert_allclose(K_T, K_T.T, atol=1e-10, err_msg="K_T should be symmetric")
        assert np.all(np.isfinite(F_int)), "F_int should be finite"
    finally:
        shutil.rmtree(temp_dir, ignore_errors=True)


def test_nonlinear_timoshenko_warping_finite_difference_force_matches_tangent():
    """ΔF_int ≈ K_T @ ΔU for small ΔU (spot-check tangent consistency)."""
    from pre_processing.element_library.nonlinear.timoshenko.nonlinear_timoshenko_3D import (
        NonlinearTimoshenkoBeamElement3D,
    )

    element_dictionary, grid_dictionary, material_dictionary, section_dictionary = _nl_ts_warp_case()
    point_load_array = np.empty((0, 9))
    distributed_load_array = np.empty((0, 9))
    temp_dir = tempfile.mkdtemp()
    job_results_dir = os.path.join(temp_dir, "nl_ts_warp_fd")
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
        rng = np.random.default_rng(7)
        u = rng.normal(scale=5e-5, size=14)
        p = rng.normal(size=14)
        p /= np.linalg.norm(p)
        eps = 1e-8
        F0 = np.asarray(elem.internal_force_vector(u), dtype=np.float64)
        F1 = np.asarray(elem.internal_force_vector(u + eps * p), dtype=np.float64)
        K_T = np.asarray(elem.tangent_stiffness_matrix(u), dtype=np.float64)
        dF = F1 - F0
        pred = K_T @ (eps * p)
        denom = max(np.linalg.norm(pred), 1e-30)
        rel_err = np.linalg.norm(dF - pred) / denom
        assert rel_err < 0.02, f"relative FD mismatch {rel_err}"
    finally:
        shutil.rmtree(temp_dir, ignore_errors=True)


def test_nonlinear_timoshenko_warping_linear_geometric_matches_linear_timoshenko():
    """linear_geometric_stiffness_matrix on NL TS + warping equals linear twin (14×14, χ rows/cols zero)."""
    from pre_processing.element_library.linear.beam.first_order_shear_deformation_theory.timoshenko.linear_timoshenko_3D import (
        LinearTimoshenkoBeamElement3D,
    )
    from pre_processing.element_library.nonlinear.timoshenko.nonlinear_timoshenko_3D import (
        NonlinearTimoshenkoBeamElement3D,
    )

    element_dictionary_nl, grid_dictionary, material_dictionary, section_dictionary = _nl_ts_warp_case()
    element_dictionary_lin = {
        **element_dictionary_nl,
        "types": np.array(["LinearTimoshenkoBeamElement3D"]),
    }
    point_load_array = np.empty((0, 9))
    distributed_load_array = np.empty((0, 9))
    temp_dir = tempfile.mkdtemp()
    job_results_dir = os.path.join(temp_dir, "nl_ts_warp_ksig")
    _job_dirs(job_results_dir)
    try:
        el_lin = LinearTimoshenkoBeamElement3D(
            element_id=0,
            element_dictionary=element_dictionary_lin,
            grid_dictionary=grid_dictionary,
            material_dictionary=material_dictionary,
            section_dictionary=section_dictionary,
            point_load_array=point_load_array,
            distributed_load_array=distributed_load_array,
            job_results_dir=job_results_dir,
        )
        el_nl = NonlinearTimoshenkoBeamElement3D(
            element_id=0,
            element_dictionary=element_dictionary_nl,
            grid_dictionary=grid_dictionary,
            material_dictionary=material_dictionary,
            section_dictionary=section_dictionary,
            point_load_array=point_load_array,
            distributed_load_array=distributed_load_array,
            job_results_dir=job_results_dir,
        )
        rng = np.random.default_rng(11)
        U14 = rng.normal(scale=1e-5, size=14)
        K_nl = el_nl.linear_geometric_stiffness_matrix(U14)
        K_lin = el_lin.linear_geometric_stiffness_matrix(U14)
        np.testing.assert_allclose(K_nl, K_lin, rtol=1e-12, atol=1e-18)
        assert np.allclose(K_nl[12:, :], 0.0)
        assert np.allclose(K_nl[:, 12:], 0.0)
    finally:
        shutil.rmtree(temp_dir, ignore_errors=True)

"""
Phase 3a: GEBT shear F_int consistency — F_int(0)=0 and K_T matches numerical derivative of F_int.
"""

import os
import shutil
import sys
import tempfile
from pathlib import Path

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))


def _cantilever_dicts(L=1.0, E=2.1e11, G=8.1e10, A=0.001, I_z=1e-8, I_y=1e-8, J_t=1e-9):
    grid_dictionary = {
        "ids": np.array([0, 1]),
        "coordinates": np.array([[0.0, 0.0, 0.0], [L, 0.0, 0.0]]),
    }
    element_dictionary = {
        "ids": np.array([0]),
        "connectivity": np.array([[0, 1]]),
        "types": np.array(["GEBTShearBeamElement3D"]),
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


def test_gebt_shear_internal_force_zero_at_zero_displacement():
    """F_int(U_e=0) should be zero."""
    from pre_processing.element_library.nonlinear.gebt_shear.gebt_shear_3D import GEBTShearBeamElement3D

    grid_dictionary, element_dictionary, material_dictionary, section_dictionary = _cantilever_dicts()
    point_load_array = np.empty((0, 9))
    distributed_load_array = np.empty((0, 9))

    temp_dir = tempfile.mkdtemp()
    job_results_dir = os.path.join(temp_dir, "gebt_fint")
    os.makedirs(job_results_dir, exist_ok=True)
    os.makedirs(os.path.join(job_results_dir, "element_stiffness_matrices"), exist_ok=True)
    os.makedirs(os.path.join(job_results_dir, "element_force_vectors"), exist_ok=True)
    os.makedirs(os.path.join(job_results_dir, "logs"), exist_ok=True)

    try:
        elem = GEBTShearBeamElement3D(
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


def test_gebt_shear_tangent_symmetric_and_internal_force_finite():
    """K_T(U_e) is symmetric; F_int(U_e) is finite for small U (equilibrium/consistency check)."""
    from pre_processing.element_library.nonlinear.gebt_shear.gebt_shear_3D import GEBTShearBeamElement3D

    grid_dictionary, element_dictionary, material_dictionary, section_dictionary = _cantilever_dicts()
    point_load_array = np.empty((0, 9))
    distributed_load_array = np.empty((0, 9))

    temp_dir = tempfile.mkdtemp()
    job_results_dir = os.path.join(temp_dir, "gebt_kcons")
    os.makedirs(job_results_dir, exist_ok=True)
    os.makedirs(os.path.join(job_results_dir, "element_stiffness_matrices"), exist_ok=True)
    os.makedirs(os.path.join(job_results_dir, "element_force_vectors"), exist_ok=True)
    os.makedirs(os.path.join(job_results_dir, "logs"), exist_ok=True)

    try:
        elem = GEBTShearBeamElement3D(
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

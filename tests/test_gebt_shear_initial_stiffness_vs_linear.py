"""
Phase 3a: At U_e=0, GEBT shear tangent stiffness K_T equals linear Timoshenko K_e.
"""

import os
import shutil
import sys
import tempfile
from pathlib import Path

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))


def _cantilever_dicts(L=2.0, E=2.1e11, G=8.1e10, A=0.00131, I_z=2.08769e-06, I_y=3.234e-07, J_t=2.60673e-08):
    grid_dictionary = {
        "ids": np.array([0, 1]),
        "coordinates": np.array([[0.0, 0.0, 0.0], [L, 0.0, 0.0]]),
    }
    element_dictionary = {
        "ids": np.array([0]),
        "connectivity": np.array([[0, 1]]),
        "types": np.array(["LinearTimoshenkoBeamElement3D"]),
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
    material_dictionary = {
        "E": np.array([E]),
        "G": np.array([G]),
        "nu": np.array([0.3]),
        "rho": np.array([7850.0]),
    }
    section_dictionary = {
        "A": np.array([A]),
        "I_x": np.array([0.0]),
        "I_y": np.array([I_y]),
        "I_z": np.array([I_z]),
        "J_t": np.array([J_t]),
    }
    return grid_dictionary, element_dictionary, material_dictionary, section_dictionary


def test_gebt_shear_initial_stiffness_vs_linear():
    """At U_e=0, GEBT shear tangent stiffness K_T equals linear Timoshenko K_e (same L, section, material)."""
    from pre_processing.element_library.linear.beam.first_order_shear_deformation_theory.timoshenko.linear_timoshenko_3D import LinearTimoshenkoBeamElement3D
    from pre_processing.element_library.nonlinear.gebt_shear.gebt_shear_3D import GEBTShearBeamElement3D

    grid_dictionary, element_dictionary, material_dictionary, section_dictionary = _cantilever_dicts()
    point_load_array = np.empty((0, 9))
    distributed_load_array = np.empty((0, 9))

    temp_dir = tempfile.mkdtemp()
    job_results_dir = os.path.join(temp_dir, "gebt_test")
    os.makedirs(job_results_dir, exist_ok=True)
    os.makedirs(os.path.join(job_results_dir, "element_stiffness_matrices"), exist_ok=True)
    os.makedirs(os.path.join(job_results_dir, "element_force_vectors"), exist_ok=True)
    os.makedirs(os.path.join(job_results_dir, "logs"), exist_ok=True)

    try:
        linear_elem = LinearTimoshenkoBeamElement3D(
            element_id=0,
            element_dictionary=element_dictionary,
            grid_dictionary=grid_dictionary,
            material_dictionary=material_dictionary,
            section_dictionary=section_dictionary,
            point_load_array=point_load_array,
            distributed_load_array=distributed_load_array,
            job_results_dir=job_results_dir,
        )
        gebt_elem = GEBTShearBeamElement3D(
            element_id=0,
            element_dictionary=element_dictionary,
            grid_dictionary=grid_dictionary,
            material_dictionary=material_dictionary,
            section_dictionary=section_dictionary,
            point_load_array=point_load_array,
            distributed_load_array=distributed_load_array,
            job_results_dir=job_results_dir,
        )

        K_linear = linear_elem.element_stiffness_matrix().K_e
        if hasattr(K_linear, "toarray"):
            K_linear = K_linear.toarray()
        K_linear = np.asarray(K_linear, dtype=np.float64)

        U_zero = np.zeros(12, dtype=np.float64)
        K_gebt = gebt_elem.tangent_stiffness_matrix(U_zero)
        K_gebt = np.asarray(K_gebt, dtype=np.float64)

        np.testing.assert_allclose(
            K_gebt, K_linear, rtol=1e-9, atol=1e-12,
            err_msg="GEBT shear initial tangent should equal linear Timoshenko stiffness",
        )
    finally:
        shutil.rmtree(temp_dir, ignore_errors=True)

"""
Verify nonlinear Timoshenko: at U=0 tangent stiffness equals linear Timoshenko K_e.
"""

import os
import shutil
import sys
import tempfile
from pathlib import Path

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))


def verify_nonlinear_timoshenko_initial_stiffness():
    """Compare nonlinear tangent at U=0 vs linear Timoshenko K_e."""
    from pre_processing.element_library.linear.beam.first_order_shear_deformation_theory.timoshenko.linear_timoshenko_3D import (
        LinearTimoshenkoBeamElement3D,
    )
    from pre_processing.element_library.nonlinear.timoshenko.nonlinear_timoshenko_3D import (
        NonlinearTimoshenkoBeamElement3D,
    )

    L, E, G = 2.0, 2.1e11, 8.1e10
    A, I_z, I_y, J_t = 0.00131, 2.08769e-06, 3.234e-07, 2.60673e-08
    grid_dictionary = {"ids": np.array([0, 1]), "coordinates": np.array([[0.0, 0.0, 0.0], [L, 0.0, 0.0]])}
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
    material_dictionary = {"E": np.array([E]), "G": np.array([G]), "nu": np.array([0.3]), "rho": np.array([7850.0])}
    section_dictionary = {
        "A": np.array([A]),
        "I_x": np.array([0.0]),
        "I_y": np.array([I_y]),
        "I_z": np.array([I_z]),
        "J_t": np.array([J_t]),
    }
    point_load_array = np.empty((0, 9))
    distributed_load_array = np.empty((0, 9))

    temp_dir = tempfile.mkdtemp()
    job_results_dir = os.path.join(temp_dir, "verify_nl_timo")
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
        nl_elem = NonlinearTimoshenkoBeamElement3D(
            element_id=0,
            element_dictionary=element_dictionary,
            grid_dictionary=grid_dictionary,
            material_dictionary=material_dictionary,
            section_dictionary=section_dictionary,
            point_load_array=point_load_array,
            distributed_load_array=distributed_load_array,
            job_results_dir=job_results_dir,
        )
        K_linear = np.asarray(linear_elem.element_stiffness_matrix().K_e, dtype=np.float64)
        K_nl = np.asarray(nl_elem.tangent_stiffness_matrix(np.zeros(12, dtype=np.float64)), dtype=np.float64)
        max_diff = np.max(np.abs(K_nl - K_linear))
        if max_diff > 1e-9:
            print("FAIL: initial tangent differs from linear Timoshenko; max |diff| =", max_diff)
            return 1
        print("PASS: nonlinear Timoshenko initial stiffness matches linear Timoshenko K_e")
        return 0
    finally:
        shutil.rmtree(temp_dir, ignore_errors=True)


if __name__ == "__main__":
    sys.exit(verify_nonlinear_timoshenko_initial_stiffness())

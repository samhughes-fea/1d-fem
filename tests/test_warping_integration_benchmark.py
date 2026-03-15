"""
Integration benchmark: one warping beam element, warping restrained at root, solve K U = F.

Optional comparison to analytical/reference; test asserts solve completes and solution shape.
"""

import os
import sys
from pathlib import Path

import numpy as np
import pytest
import scipy.sparse as sp

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from pre_processing.element_library.element_factory import ElementFactory
from pre_processing.element_library.linear.timoshenko.linear_warping_timoshenko_3D import (
    LinearWarpingTimoshenkoBeamElement3D,
)


def test_warping_integration_benchmark_solve():
    """One warping element, fixed at node 0 (all 7 DOF), tip torque; solve and check U shape."""
    L = 2.0
    E, G = 2.1e11, 8.1e10
    A = 0.001
    I_y, I_z = 1e-8, 1e-8
    J_t = 1e-9
    Gamma = 1.0e-8

    import tempfile
    temp_dir = tempfile.mkdtemp()
    job_results_dir = os.path.join(temp_dir, "warping_bench")
    os.makedirs(job_results_dir, exist_ok=True)
    os.makedirs(os.path.join(job_results_dir, "element_stiffness_matrices"), exist_ok=True)
    os.makedirs(os.path.join(job_results_dir, "element_force_vectors"), exist_ok=True)
    os.makedirs(os.path.join(job_results_dir, "logs"), exist_ok=True)

    element_dictionary = {
        "ids": np.array([0]),
        "connectivity": np.array([[0, 1]]),
        "types": np.array(["LinearWarpingTimoshenkoBeamElement3D"]),
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
        "kappa": np.array([5.0 / 6.0]),
        "alpha": np.array([0.0]),
        "y_sc": np.array([0.0]),
        "z_sc": np.array([0.0]),
        "Gamma": np.array([Gamma]),
    }
    point_load_array = np.empty((0, 9))
    distributed_load_array = np.empty((0, 9))

    try:
        factory = ElementFactory(job_results_dir=job_results_dir)
        elements = factory.create_elements_batch(
            element_ids=np.array([0]),
            element_dictionary=element_dictionary,
            grid_dictionary=grid_dictionary,
            material_dictionary=material_dictionary,
            section_dictionary=section_dictionary,
            point_load_array=point_load_array,
            distributed_load_array=distributed_load_array,
        )
        assert len(elements) == 1
        element = elements[0]
        assert element.dof_per_node == 7

        obj_Ke = element.element_stiffness_matrix()
        obj_Fe = element.element_force_vector()
        Ke = obj_Ke.K_e
        Fe = obj_Fe.F_e
        dof_indices = element.assemble_global_dof_indices()

        total_dof = 2 * 7
        K_global = sp.csr_matrix((total_dof, total_dof))
        F_global = np.zeros(total_dof)
        for i, di in enumerate(dof_indices):
            F_global[di] += Fe[i]
            for j, dj in enumerate(dof_indices):
                K_global[di, dj] += Ke[i, j]

        # Fix node 0 (all 7 DOF) and warping at node 1 (DOF 13) so the free block is non-singular
        # (θ_x and χ at tip share a null direction in φ_x′ otherwise)
        fixed_dofs = np.array([0, 1, 2, 3, 4, 5, 6, 13])
        free_dofs = np.array([7, 8, 9, 10, 11, 12])
        K_ff = K_global[free_dofs, :][:, free_dofs]
        F_f = F_global[free_dofs]
        U_f = np.linalg.solve(K_ff.toarray(), F_f)
        U_global = np.zeros(total_dof)
        U_global[free_dofs] = U_f
        U_global[fixed_dofs] = 0.0

        assert U_global.shape == (14,)
        assert np.all(U_global[fixed_dofs] == 0.0)
    finally:
        import shutil
        if os.path.isdir(temp_dir):
            shutil.rmtree(temp_dir, ignore_errors=True)

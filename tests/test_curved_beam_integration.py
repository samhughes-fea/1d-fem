"""
Integration test: semicircular arch (or curved segment), fixed at one end; solve and check solution shape.
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
from pre_processing.element_library.linear.beam.first_order_shear_deformation_theory.curved_beam.linear_curved_timoshenko_3D import (
    LinearCurvedTimoshenkoBeamElement3D,
)


def test_curved_beam_integration_solve():
    """One curved element (κ0 > 0), fixed at node 0; assemble and solve; assert U shape and fixed DOFs zero."""
    L = 2.0
    R = 5.0
    kappa0 = 1.0 / R
    E, G = 2.1e11, 8.1e10
    A, I_y, I_z, J_t = 0.001, 1e-8, 1e-8, 1e-9

    import tempfile
    temp_dir = tempfile.mkdtemp()
    job_results_dir = os.path.join(temp_dir, "curved_integration")
    os.makedirs(job_results_dir, exist_ok=True)
    os.makedirs(os.path.join(job_results_dir, "element_stiffness_matrices"), exist_ok=True)
    os.makedirs(os.path.join(job_results_dir, "element_force_vectors"), exist_ok=True)
    os.makedirs(os.path.join(job_results_dir, "logs"), exist_ok=True)

    element_dictionary = {
        "ids": np.array([0]),
        "connectivity": np.array([[0, 1]]),
        "types": np.array(["LinearCurvedTimoshenkoBeamElement3D"]),
        "integration_orders": {
            "axial": np.array([3]),
            "bending_y": np.array([3]),
            "bending_z": np.array([3]),
            "shear_y": np.array([2]),
            "shear_z": np.array([2]),
            "torsion": np.array([2]),
            "load": np.array([2]),
        },
        "curvature": np.array([kappa0]),
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
        assert element.curvature == kappa0

        obj_Ke = element.element_stiffness_matrix()
        obj_Fe = element.element_force_vector()
        Ke = obj_Ke.K_e
        Fe = obj_Fe.F_e
        dof_indices = element.assemble_global_dof_indices()

        total_dof = 2 * 6
        K_global = sp.lil_matrix((total_dof, total_dof))
        F_global = np.zeros(total_dof)
        for i, di in enumerate(dof_indices):
            F_global[di] += Fe[i]
            for j, dj in enumerate(dof_indices):
                K_global[di, dj] += Ke[i, j]
        K_global = K_global.tocsr()

        fixed_dofs = np.arange(0, 6)
        free_dofs = np.arange(6, 12)
        K_ff = K_global[free_dofs, :][:, free_dofs]
        F_f = F_global[free_dofs]
        U_f = np.linalg.solve(K_ff.toarray(), F_f)
        U_global = np.zeros(total_dof)
        U_global[free_dofs] = U_f

        assert U_global.shape == (12,)
        np.testing.assert_allclose(U_global[fixed_dofs], 0.0, atol=1e-12)
    finally:
        import shutil
        if os.path.isdir(temp_dir):
            shutil.rmtree(temp_dir, ignore_errors=True)

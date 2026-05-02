"""Cantilever with nonlinear Timoshenko: Newton-style solve smoke test."""

import os
import shutil
import sys
import tempfile
from pathlib import Path

import numpy as np
import scipy.sparse as sp

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))


def _cantilever_dicts(L=2.0, E=2.1e11, G=8.1e10, A=0.00131, I_z=2.08769e-06, I_y=3.234e-07, J_t=2.60673e-08, P=-100.0):
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
    point_load_array = np.array([[L, 0.0, 0.0, 0.0, P, 0.0, 0.0, 0.0, 0.0]])
    distributed_load_array = np.empty((0, 9))
    return (
        grid_dictionary,
        element_dictionary,
        material_dictionary,
        section_dictionary,
        point_load_array,
        distributed_load_array,
    )


def test_nonlinear_timoshenko_integration_cantilever_solve():
    """Assemble K_T and F_int at U=0, one Newton step; fixed DOFs remain zero."""
    from pre_processing.element_library.element_factory import ElementFactory
    from pre_processing.element_library.nonlinear.timoshenko.nonlinear_timoshenko_3D import (
        NonlinearTimoshenkoBeamElement3D,
    )

    (
        grid_dictionary,
        element_dictionary,
        material_dictionary,
        section_dictionary,
        point_load_array,
        distributed_load_array,
    ) = _cantilever_dicts()

    temp_dir = tempfile.mkdtemp()
    job_results_dir = os.path.join(temp_dir, "nl_timo_integration")
    os.makedirs(job_results_dir, exist_ok=True)
    os.makedirs(os.path.join(job_results_dir, "element_stiffness_matrices"), exist_ok=True)
    os.makedirs(os.path.join(job_results_dir, "element_force_vectors"), exist_ok=True)
    os.makedirs(os.path.join(job_results_dir, "logs"), exist_ok=True)

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
        elem = elements[0]
        assert isinstance(elem, NonlinearTimoshenkoBeamElement3D)

        U_e = np.zeros(12, dtype=np.float64)
        K_T = elem.tangent_stiffness_matrix(U_e)
        F_int = elem.internal_force_vector(U_e)
        F_ext = elem.element_force_vector().F_e
        R = np.asarray(F_ext, dtype=np.float64).ravel() - np.asarray(F_int, dtype=np.float64).ravel()

        dof = elem.assemble_global_dof_indices()
        total_dof = 12
        K_global = sp.lil_matrix((total_dof, total_dof))
        for i, di in enumerate(dof):
            for j, dj in enumerate(dof):
                K_global[di, dj] += K_T[i, j]
        K_global = K_global.tocsr()
        fixed = np.arange(0, 6)
        free = np.arange(6, 12)
        K_ff = K_global[free, :][:, free]
        R_f = R[free]
        dU_f = np.linalg.solve(K_ff.toarray(), R_f)
        U_global = np.zeros(total_dof)
        U_global[free] = dU_f

        assert U_global.shape == (12,)
        np.testing.assert_allclose(U_global[fixed], 0.0, atol=1e-12)
        assert np.all(np.isfinite(U_global))
    finally:
        shutil.rmtree(temp_dir, ignore_errors=True)

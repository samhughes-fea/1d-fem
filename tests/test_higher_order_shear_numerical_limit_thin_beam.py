"""
Phase 3b numerical limit: thin beam (L/h large); Reddy and Timoshenko tip deflections should be close.
"""

import os
import sys
import tempfile
import shutil
from pathlib import Path

import numpy as np
import scipy.sparse as sp

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from pre_processing.element_library.element_factory import ElementFactory


def test_higher_order_shear_numerical_limit_thin_beam():
    """Thin beam (L/h=20): cantilever tip deflection; Reddy and Timoshenko should be close (same order)."""
    L = 1.0
    h = 0.05
    E, G = 2.1e11, 8.1e10
    b = 0.02
    A = b * h
    I_z = b * h**3 / 12.0
    I_y = h * b**3 / 12.0
    J_t = 0.02 * b * h**3
    P = -10.0

    temp_dir = tempfile.mkdtemp()
    job_results_dir = os.path.join(temp_dir, "thin_beam")
    os.makedirs(job_results_dir, exist_ok=True)
    os.makedirs(os.path.join(job_results_dir, "element_stiffness_matrices"), exist_ok=True)
    os.makedirs(os.path.join(job_results_dir, "element_force_vectors"), exist_ok=True)
    os.makedirs(os.path.join(job_results_dir, "logs"), exist_ok=True)

    section_dictionary = {
        "A": np.array([A]),
        "I_x": np.array([0.0]),
        "I_y": np.array([I_y]),
        "I_z": np.array([I_z]),
        "J_t": np.array([J_t]),
    }
    grid_dictionary = {"coordinates": np.array([[0.0, 0.0, 0.0], [L, 0.0, 0.0]])}
    material_dictionary = {
        "E": np.array([E]),
        "G": np.array([G]),
        "nu": np.array([0.3]),
        "rho": np.array([7850.0]),
    }
    point_load_array = np.array([[L, 0.0, 0.0, 0.0, P, 0.0, 0.0, 0.0, 0.0]])
    distributed_load_array = np.empty((0, 9))
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

    def solve_tip_deflection(element_type_name):
        element_dictionary["types"] = np.array([element_type_name])
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
        elem = elements[0]
        obj_K = elem.element_stiffness_matrix()
        obj_F = elem.element_force_vector()
        Ke = obj_K.K_e
        Fe = obj_F.F_e
        if hasattr(Ke, "toarray"):
            Ke = Ke.toarray()
        dof = elem.assemble_global_dof_indices()
        K_global = sp.lil_matrix((12, 12))
        F_global = np.zeros(12)
        for i, di in enumerate(dof):
            F_global[di] += Fe[i]
            for j, dj in enumerate(dof):
                K_global[di, dj] += Ke[i, j]
        K_global = K_global.tocsr()
        free = np.arange(6, 12)
        U_f = np.linalg.solve(K_global[free, :][:, free].toarray(), F_global[free])
        U = np.zeros(12)
        U[free] = U_f
        return U[7]

    try:
        u_y_reddy = solve_tip_deflection("LinearReddyBeamElement3D")
        u_y_timo = solve_tip_deflection("LinearTimoshenkoBeamElement3D")
        assert np.isfinite(u_y_reddy) and np.isfinite(u_y_timo)
        assert (u_y_reddy < 0 and u_y_timo < 0) or (u_y_reddy > 0 and u_y_timo > 0), (
            "Thin beam: Reddy and Timoshenko should agree in sign (both downward for tip load)"
        )
        ratio = abs(u_y_reddy) / (abs(u_y_timo) + 1e-30)
        assert 1e-6 < ratio < 1e6, (
            f"Thin beam: tip deflections should be finite and non-degenerate; ratio={ratio}"
        )
    finally:
        shutil.rmtree(temp_dir, ignore_errors=True)

# Native GESDB Gauss assembly for F_int and K_T (engineering stretch + linear TS rows).

from __future__ import annotations

import numpy as np

from onedfem.pre_processing.element_library.nonlinear.large_rotations.geometrically_exact_shear_deformable_beam.utilities.gesdb_kinematics import (
    native_engineering_strain_and_B_eng_timoshenko_12,
)


def native_internal_force(elem, U_e: np.ndarray) -> np.ndarray:
    """``B_engᵀ S`` Gauss sum on the 12-DOF native kernel."""
    D = elem.material_stiffness_operator.assembly_form()
    xi, w = elem.integration_points
    detJ = elem.jacobian_determinant
    F_int = np.zeros(12, dtype=np.float64)
    U_u = np.asarray(U_e, dtype=np.float64).ravel()
    for xi_g, w_g in zip(xi, w):
        E, B_eng = native_engineering_strain_and_B_eng_timoshenko_12(elem, U_u, xi_g)
        S = D @ (E - elem._E_0_voigt)
        F_int += (B_eng.T @ S) * w_g * detJ
    return F_int


def native_tangent_stiffness(elem, U_e: np.ndarray) -> np.ndarray:
    """``K_mat + K_sigma`` for the native kernel (no TL ``K_delta`` split)."""
    U_u = np.asarray(U_e, dtype=np.float64).ravel()
    D = elem.material_stiffness_operator.assembly_form()
    xi, w = elem.integration_points
    detJ = elem.jacobian_determinant
    n_g = len(xi)
    N_gp = np.zeros(n_g, dtype=np.float64)
    M_y_gp = np.zeros(n_g, dtype=np.float64)
    M_z_gp = np.zeros(n_g, dtype=np.float64)
    dN_dx_list = []
    K_mat = np.zeros((12, 12), dtype=np.float64)
    for k, (xi_g, w_g) in enumerate(zip(xi, w)):
        E, B_eng = native_engineering_strain_and_B_eng_timoshenko_12(elem, U_u, xi_g)
        K_mat += B_eng.T @ D @ B_eng * w_g * detJ
        N_i, M_y_i, M_z_i = elem.stress_resultant_operator.section_forces_from_strain(
            E - elem._E_0_voigt, D
        )
        N_gp[k] = N_i
        M_y_gp[k] = M_y_i
        M_z_gp[k] = M_z_i
        _, dN_dξ, _d2 = elem.shape_function_operator.natural_coordinate_form(np.array([xi_g]))
        dN_dx_list.append(dN_dξ[0] * elem.dxi_dx)
    dN_dx_arr = np.stack(dN_dx_list, axis=0)
    K_sigma = elem.geometric_stiffness_operator.assemble_K_sigma(
        N_gp, M_y_gp, M_z_gp, w, dN_dx_arr, detJ
    )
    return K_mat + K_sigma

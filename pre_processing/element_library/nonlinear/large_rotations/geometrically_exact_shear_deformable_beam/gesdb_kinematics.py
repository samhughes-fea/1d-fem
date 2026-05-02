# pre_processing/element_library/nonlinear/large_rotations/geometrically_exact_shear_deformable_beam/gesdb_kinematics.py
"""
GESDB strain kernels and mapping to the repository Voigt layout.

The reference weak form is summarised in ``docs/element_library/gesdb_weak_form.md``.
"""

from __future__ import annotations

from typing import Any, Tuple

import numpy as np

from pre_processing.element_library.nonlinear.timoshenko.tl_green_lagrange_voigt import (
    chord_frame_green_lagrange_voigt_timoshenko_12,
)


def gesdb_director_voigt_strain_timoshenko_12(
    elem: Any, U_e: np.ndarray, xi_g: float
) -> np.ndarray:
    """
    Director-based Voigt strains at ``xi_g`` for the 2-node Timoshenko GESDB element (12 DOF).

    For the current shape-function choice and fixed chord frame, the work-conjugate strain vector
    implemented in production matches the chord-frame Green–Lagrange reduction used by
    :class:`~pre_processing.element_library.nonlinear.timoshenko.nonlinear_timoshenko_3D.NonlinearTimoshenkoBeamElement3D`
    (see ``gesdb_weak_form.md``: locked reference until a distinct SVQ ``B`` is wired).

    Parameters
    ----------
    elem, U_e, xi_g
        Same contract as :func:`~pre_processing.element_library.nonlinear.timoshenko.tl_green_lagrange_voigt.chord_frame_green_lagrange_voigt_timoshenko_12`.

    Returns
    -------
    np.ndarray
        Six-component Voigt strain at the Gauss point.
    """
    return chord_frame_green_lagrange_voigt_timoshenko_12(elem, U_e, xi_g)


def native_engineering_strain_and_B_eng_timoshenko_12(
    elem: Any, U_e: np.ndarray, xi_g: float
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Native GESDB pathway: engineering axial stretch plus linear Timoshenko rows 1–5.

    Axial strain uses ``||r2−r1|| / L_ref − 1`` from deformed chord ``r_i = X_i + u_i``.
    Rows ``κ_y``, ``κ_z``, ``γ_xy``, ``γ_xz``, ``φ_x`` match ``B_lin @ U`` at ``xi_g``.

    The resulting ``B_eng`` is the Jacobian ``∂E/∂U`` so that
    ``F_int = Σ w_g |J| B_engᵀ S`` and ``K_mat ≈ Σ w_g |J| B_engᵀ D B_eng`` are consistent with
    ``σ = D (E − E₀)`` for this strain vector (see ``gesdb_weak_form.md``, native kernel section).

    Parameters
    ----------
    elem
        Nonlinear Timoshenko-like element with ``node_coords``, ``L``, shape and strain operators.
    U_e
        Element displacement, shape ``(12,)``.
    xi_g
        Natural coordinate.

    Returns
    -------
    E : np.ndarray, shape (6,)
        Voigt strain vector.
    B_eng : np.ndarray, shape (6, 12)
        Strain-displacement Jacobian for the native vector.
    """
    U12 = np.asarray(U_e, dtype=np.float64).ravel()[:12]
    xi_g = float(xi_g)
    X1 = elem.node_coords[0]
    X2 = elem.node_coords[1]
    r1 = X1 + U12[0:3]
    r2 = X2 + U12[6:9]
    d_vec = r2 - r1
    Ld = float(np.linalg.norm(d_vec))
    Lref = float(elem.L)
    if Lref <= 0.0:
        raise ValueError("Reference length L must be positive.")
    eps_ax = Ld / Lref - 1.0

    N, dN_dξ, d2N_dξ2 = elem.shape_function_operator.natural_coordinate_form(np.array([xi_g]))
    B_lin = elem.strain_displacement_operator.physical_coordinate_form(dN_dξ, d2N_dξ2, N)[0]
    e_lin = (B_lin @ U12).ravel()
    E = np.array(
        [eps_ax, e_lin[1], e_lin[2], e_lin[3], e_lin[4], e_lin[5]],
        dtype=np.float64,
    )

    B_eng = np.zeros((6, 12), dtype=np.float64)
    if Ld > 1e-16:
        t_hat = d_vec / Ld
        B_eng[0, 0:3] = -t_hat / Lref
        B_eng[0, 6:9] = t_hat / Lref
    B_eng[1:6, :] = B_lin[1:6, :]
    return E, B_eng


def native_engineering_voigt_strain_timoshenko_12(
    elem: Any, U_e: np.ndarray, xi_g: float
) -> np.ndarray:
    """Voigt strain only (tests / strain hook)."""
    E, _ = native_engineering_strain_and_B_eng_timoshenko_12(elem, U_e, xi_g)
    return E

# pre_processing/element_library/nonlinear/timoshenko/tl_green_lagrange_voigt.py
"""
Chord-frame Green–Lagrange Voigt strain for 12-DOF Total Lagrangian Timoshenko.

Used by :class:`NonlinearTimoshenkoBeamElement3D` and referenced from the GESDB strain hook
(:doc:`docs/element_library/gesdb_weak_form.md`). Voigt ordering follows
``FORMULATION_DOCSTRING_STANDARDS.md``.
"""

from __future__ import annotations

from typing import Any

import numpy as np


def chord_frame_green_lagrange_voigt_timoshenko_12(
    elem: Any, U_e: np.ndarray, xi_g: float
) -> np.ndarray:
    """
    Total Lagrangian Voigt strain ``E = E_lin + E_nl`` at natural coordinate ``xi_g``.

    Parameters
    ----------
    elem
        Nonlinear Timoshenko element exposing ``shape_function_operator``,
        ``strain_displacement_operator``, ``green_lagrange_strain_operator``, ``dxi_dx``.
    U_e
        Element displacement, shape ``(12,)`` (first 12 entries used if longer).
    xi_g
        Gauss point in ``[-1, 1]``.

    Returns
    -------
    np.ndarray
        Strain vector, shape ``(6,)``.
    """
    U12 = np.asarray(U_e, dtype=np.float64).ravel()
    if U12.size < 12:
        raise ValueError("Expected at least 12 displacement DOFs for Timoshenko TL strain.")
    U12 = U12[:12]
    xi_g = float(xi_g)
    N, dN_dξ, d2N_dξ2 = elem.shape_function_operator.natural_coordinate_form(np.array([xi_g]))
    B = elem.strain_displacement_operator.physical_coordinate_form(dN_dξ, d2N_dξ2, N)[0]
    E_lin = (B @ U12).ravel()
    dN_dx = dN_dξ[0] * elem.dxi_dx
    d2N_dx2 = d2N_dξ2[0] * (elem.dxi_dx**2)
    E_nl = elem.green_lagrange_strain_operator.strain_nonlinear_part(
        dN_dx, d2N_dx2, U12, N[0]
    )
    return np.asarray(E_lin + E_nl, dtype=np.float64).ravel()

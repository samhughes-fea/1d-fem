# pre_processing/element_library/utilities/total_lagrangian_beam.py
"""
Legacy shim for Total Lagrangian beam helpers.

**Prefer** importing ``GreenLagrangeStrainOperator`` from:

- ``pre_processing.element_library.nonlinear.euler_bernoulli.utilities.green_lagrange_strain`` (Euler–Bernoulli NL), or
- ``pre_processing.element_library.nonlinear.timoshenko.utilities.green_lagrange_strain`` (Timoshenko NL).

Those modules carry the implemented strain splits (\ ``E = E_lin + E_nl``\ ), ``B_nl``, and docstrings aligned with
``nonlinear_*_3D`` elements. This file keeps a minimal frozen dataclass API for older imports; it is **not** the
canonical formulation source.

Operators: ``GreenLagrangeStrainOperator`` (partial), ``StressResultantOperator``, ``GeometricStiffnessOperator``.
Tangent patterns for elements: NL EB uses ``K_mat + K_sigma`` with ``B_tot``; NL Timoshenko uses ``K_0 + K_delta + K_sigma``.
"""

from dataclasses import dataclass
from typing import Tuple

import numpy as np

from pre_processing.element_library.nonlinear.euler_bernoulli.utilities.geometric_stiffness import (
    GeometricStiffnessOperator,
)


@dataclass(frozen=True)
class GreenLagrangeStrainOperator:
    """
    Minimal Green–Lagrange strain helper (legacy).

    For nonlinear beam elements, use the dedicated EB or Timoshenko ``GreenLagrangeStrainOperator`` under
    ``nonlinear/.../utilities/`` — they implement ``strain_nonlinear_part(dN_dx, d2N_dx2, u_e)`` and ``B_nl``.
    """

    element_length: float
    include_shear: bool = False

    def __post_init__(self) -> None:
        if self.element_length <= 0:
            raise ValueError(f"element_length must be positive, got {self.element_length}")

    def strain_linear_part(
        self,
        dN_dx: np.ndarray,
        d2N_dx2: np.ndarray,
        u_e: np.ndarray,
    ) -> np.ndarray:
        """``E_lin = B_lin @ u_e`` using ``linearized_strain_displacement``."""
        B_lin = self.linearized_strain_displacement(dN_dx, d2N_dx2)
        return (B_lin @ u_e).ravel()

    def strain_nonlinear_part(
        self,
        dN_dx: np.ndarray,
        u_e: np.ndarray,
    ) -> np.ndarray:
        """Axial nonlinear strain only (½‖∇u‖² along x); no ``d2N`` / curvature NL — use EB/Timoshenko utilities instead."""
        du_dx = dN_dx[0, 0] * u_e[0] + dN_dx[6, 0] * u_e[6]
        dv_dx = dN_dx[1, 1] * u_e[1] + dN_dx[7, 1] * u_e[7]
        dw_dx = dN_dx[2, 2] * u_e[2] + dN_dx[8, 2] * u_e[8]
        e_nl = np.zeros(6, dtype=np.float64)
        e_nl[0] = 0.5 * (du_dx**2 + dv_dx**2 + dw_dx**2)
        return e_nl

    def linearized_strain_displacement(
        self,
        dN_dx: np.ndarray,
        d2N_dx2: np.ndarray,
        N: np.ndarray | None = None,
    ) -> np.ndarray:
        """``B_lin`` (6, 12); EB-style curvature rows when shear off."""
        B_lin = np.zeros((6, 12), dtype=np.float64)
        B_lin[0, 0] = dN_dx[0, 0]
        B_lin[0, 6] = dN_dx[6, 0]
        B_lin[1, 2] = d2N_dx2[2, 2]
        B_lin[1, 8] = d2N_dx2[8, 2]
        B_lin[1, 4] = d2N_dx2[4, 4]
        B_lin[1, 10] = d2N_dx2[10, 4]
        B_lin[2, 1] = d2N_dx2[1, 1]
        B_lin[2, 7] = d2N_dx2[7, 1]
        B_lin[2, 5] = d2N_dx2[5, 5]
        B_lin[2, 11] = d2N_dx2[11, 5]
        if self.include_shear:
            B_lin[3, 1] = dN_dx[1, 1]
            B_lin[3, 7] = dN_dx[7, 1]
            B_lin[3, 5] = -N[5, 5] if N is not None else -1.0
            B_lin[3, 11] = -N[11, 5] if N is not None else -1.0
            B_lin[4, 2] = dN_dx[2, 2]
            B_lin[4, 8] = dN_dx[8, 2]
            B_lin[4, 4] = -N[4, 4] if N is not None else -1.0
            B_lin[4, 10] = -N[10, 4] if N is not None else -1.0
        B_lin[5, 3] = dN_dx[3, 3]
        B_lin[5, 9] = dN_dx[9, 3]
        return B_lin


@dataclass(frozen=True)
class StressResultantOperator:
    """
    Section forces from ``S = D @ E`` at a point (same as EB/Timoshenko utilities).
    """

    def section_forces_from_strain(
        self,
        E: np.ndarray,
        D: np.ndarray,
    ) -> Tuple[float, float, float]:
        """Return ``N, M_y, M_z`` from first three stress components."""
        S = D @ E
        N = float(S[0])
        M_y = float(S[1])
        M_z = float(S[2])
        return N, M_y, M_z

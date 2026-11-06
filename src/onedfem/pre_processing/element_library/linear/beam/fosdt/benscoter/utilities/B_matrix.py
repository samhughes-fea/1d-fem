# pre_processing/element_library/linear/beam/fosdt/benscoter/utilities/B_matrix.py
"""Strain-displacement **B** (8, 14) — infinitesimal Benscoter kinematics.

**Discrete:** ``eps = B @ U_e`` with
``eps = [eps_x, kappa_y, kappa_z, gamma_xy, gamma_xz, phi_x, d, psi']``,
``d = psi - theta_x'``, ``psi' = d psi / dx``.

Rows 0–5 reuse the parent Timoshenko ``StrainDisplacementOperator`` (6×12 block).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


_N_STD = 12
_N_STRAIN = 8


@dataclass(frozen=True)
class StrainDisplacementOperator:
    """Benscoter **B** at Gauss stations; composes Timoshenko 6×12 kinematics."""

    element_length: float
    n_dof: int
    psi_cols: tuple
    parent: object

    def physical_coordinate_form(
        self,
        dN_dξ: np.ndarray,
        d2N_dξ2: np.ndarray,
        N: np.ndarray,
        xi=None,
    ) -> np.ndarray:
        n_gauss = dN_dξ.shape[0]
        B6 = self.parent.physical_coordinate_form(dN_dξ, d2N_dξ2, N)
        xi = np.zeros(n_gauss) if xi is None else np.atleast_1d(np.asarray(xi, dtype=np.float64))
        p1, p2 = self.psi_cols
        B = np.zeros((n_gauss, _N_STRAIN, self.n_dof), dtype=np.float64)
        B[:, :6, :_N_STD] = B6
        L = float(self.element_length)
        for g in range(n_gauss):
            B[g, 6, :] = -B[g, 5, :]
            B[g, 6, p1] += 0.5 * (1.0 - xi[g])
            B[g, 6, p2] += 0.5 * (1.0 + xi[g])
            B[g, 7, p1] = -1.0 / L
            B[g, 7, p2] = 1.0 / L
        return B

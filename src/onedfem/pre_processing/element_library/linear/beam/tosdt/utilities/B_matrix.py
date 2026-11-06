# pre_processing/element_library/linear/beam/tosdt/utilities/B_matrix.py
"""Strain-displacement **B** (8, 16) — infinitesimal third-order shear kinematics.

**Discrete:** ``e8 = B @ U_e`` with Hermite ``v``, ``w`` and linear rotations.
``chi = gamma' - R psi6``. Element signs: ``S = diag(1,-1,1,1,1,1,1,1)`` on rows.

**Stiffness:** Reddy ``K = int B8.T D B8 dx``; Levinson ``K = int B6.T D[:6,:] B8 dx``.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from onedfem.pre_processing.element_library.linear.beam.tosdt.utilities.shape_functions import (
    N_DOF,
    _V_COLS,
    _W_COLS,
    hermite_x,
)

_SIGN8 = np.array([1.0, -1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0])


@dataclass(frozen=True)
class StrainDisplacementOperator:
    element_length: float
    rate: np.ndarray

    def at(self, xi: float) -> np.ndarray:
        L = float(self.element_length)
        L1, L2 = 0.5 * (1.0 - xi), 0.5 * (1.0 + xi)
        dL = np.array([-1.0, 1.0]) / L
        _, dH, d2H = hermite_x(float(xi), L)
        B = np.zeros((8, N_DOF))
        B[0, [0, 6]] = dL
        B[1, [4, 10]] = dL
        B[2, [5, 11]] = dL
        B[3, _V_COLS] = dH
        B[3, [5, 11]] = [-L1, -L2]
        B[4, _W_COLS] = dH
        B[4, [4, 10]] = [L1, L2]
        B[5, [3, 9]] = dL
        B[6, _V_COLS] = d2H
        B[6, [5, 11]] = -dL
        B[7, _W_COLS] = d2H
        B[7, [4, 10]] = dL
        B *= _SIGN8[:, None]
        B[6:] -= np.asarray(self.rate, dtype=np.float64) @ B[:6]
        return B

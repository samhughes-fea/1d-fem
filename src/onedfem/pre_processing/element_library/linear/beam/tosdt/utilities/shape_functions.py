# pre_processing/element_library/linear/beam/tosdt/utilities/shape_functions.py
"""Hermite field interpolation for 16-DOF third-order (Reddy/Levinson) beams.

``N_q`` maps ``U_e`` (16,) to ``[u, v, w, theta_x, theta_y, theta_z]``.
``v`` and ``w`` are cubic Hermite; the remaining fields are linear Lagrange.
Stiffness **B** is built in ``B_matrix.py`` from these derivatives — not from the
Timoshenko registry operator.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

N_DOF = 16
_V_COLS = [1, 12, 7, 14]
_W_COLS = [2, 13, 8, 15]


def hermite_x(xi: float, L: float):
    """Cubic Hermite on ``x = (1 + xi) L / 2``: values, d/dx, d2/dx2 for ``(f1, f1', f2, f2')``."""
    x = 0.5 * (1.0 + xi)
    H = np.array([1 - 3 * x**2 + 2 * x**3, L * (x - 2 * x**2 + x**3), 3 * x**2 - 2 * x**3, L * (-x**2 + x**3)])
    dH = np.array([-6 * x + 6 * x**2, L * (1 - 4 * x + 3 * x**2), 6 * x - 6 * x**2, L * (-2 * x + 3 * x**2)]) / L
    d2H = np.array([-6 + 12 * x, L * (-4 + 6 * x), 6 - 12 * x, L * (-2 + 6 * x)]) / L**2
    return H, dH, d2H


@dataclass(frozen=True)
class ShapeFunctionOperator:
    """Field interpolation **N_q** (6, 16) at a natural coordinate."""

    element_length: float

    def field_interpolation(self, xi: float, element_signs: bool = True) -> np.ndarray:
        L1, L2 = 0.5 * (1.0 - xi), 0.5 * (1.0 + xi)
        H, _, _ = hermite_x(float(xi), self.element_length)
        Nq = np.zeros((6, N_DOF))
        for r in (0, 3, 4, 5):
            Nq[r, [r, r + 6]] = [L1, L2]
        Nq[1, _V_COLS] = H
        Nq[2, _W_COLS] = H
        if element_signs:
            Nq[4] *= -1.0
        return Nq

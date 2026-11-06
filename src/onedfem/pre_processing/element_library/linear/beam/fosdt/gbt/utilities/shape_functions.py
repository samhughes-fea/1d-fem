# pre_processing/element_library/linear/beam/fosdt/gbt/utilities/shape_functions.py
"""Hermite interpolation for GBT distortional mode amplitudes along the span.

``a_k(x)`` uses cubic Hermite on ``xi in [-1, 1]`` with nodal ``(a, a')``.
This is spanwise 1D enrichment, not a section-plane shape function.
"""

from __future__ import annotations

import numpy as np


def gbt_modes_of_row(element_dictionary: dict, idx: int, key: str = "gbt_modes") -> int:
    modes = element_dictionary.get(key)
    return 0 if modes is None else int(np.asarray(modes).ravel()[idx])


def hermite_x(xi: float, L: float):
    """Hermite shapes for ``(a_1, a_1', a_2, a_2')`` on ``x = (1 + xi) L / 2``: values, d/dx, d2/dx2."""
    x = 0.5 * (1.0 + xi)
    H = np.array([1 - 3 * x**2 + 2 * x**3, L * (x - 2 * x**2 + x**3), 3 * x**2 - 2 * x**3, L * (-x**2 + x**3)])
    dH = np.array([-6 * x + 6 * x**2, L * (1 - 4 * x + 3 * x**2), 6 * x - 6 * x**2, L * (-2 * x + 3 * x**2)]) / L
    d2H = np.array([-6 + 12 * x, L * (-4 + 6 * x), 6 - 12 * x, L * (-2 + 6 * x)]) / L**2
    return H, dH, d2H

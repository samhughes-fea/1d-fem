# post_processing/graphical_visualisers/beam_field_sampling.py
"""Displacement fields along a straight beam from its nodal solution, through the element shape functions.

A curve of deflection or twist along a beam is the finite element field
``[u, v, w, theta_x, theta_y, theta_z](xi) = N(xi) U_e`` of each element, not a line drawn
between nodal values. ``beam_curve`` samples that field inside every element with the
element's own interpolation (``beam_load_interpolation.beam_field_interpolation``: cubic
Hermite bending for Euler-Bernoulli elements, linear Lagrange for Timoshenko, Benscoter and
GBT elements, as ``framed_element._bending_kind`` assigns) and returns the samples with a
``NaN`` between elements, so a plot shows each element's interpolant on its own.

Public API
----------
interpolation_kind, beam_curve.
"""

from __future__ import annotations

from typing import Sequence, Tuple

import numpy as np
import numpy.typing as npt

from pre_processing.element_library.beam_load_interpolation import beam_field_interpolation

__all__ = ["interpolation_kind", "beam_curve"]


def interpolation_kind(element_type: str) -> str:
    """``"hermite"`` for Euler-Bernoulli element types, else ``"linear"`` (as ``_bending_kind``)."""
    return "hermite" if "EulerBernoulli" in element_type else "linear"


def beam_curve(x_nodes: Sequence[float], U: npt.ArrayLike, component: int, element_type: str,
               n_per_element: int = 11) -> Tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]:
    """``(x, value)`` of one field component along a straight beam on the ``x`` axis.

    ``U`` (n_nodes, >= 6) holds the right-handed nodal DOFs ``u, v, w, theta_x, theta_y,
    theta_z`` in its first six columns (extra DOFs such as warping are not used); node ``i``
    and ``i + 1`` bound element ``i``. ``component`` indexes ``u, v, w, theta_x, theta_y,
    theta_z``. Elements are separated by ``NaN``.
    """
    x = np.asarray(x_nodes, dtype=np.float64)
    U = np.asarray(U, dtype=np.float64)[:, :6]
    kind = interpolation_kind(element_type)
    xi = np.linspace(-1.0, 1.0, n_per_element)
    xs, vs = [], []
    for e in range(x.size - 1):
        L = float(x[e + 1] - x[e])
        Ue = np.concatenate([U[e], U[e + 1]])
        xs.append(x[e] + 0.5 * (1.0 + xi) * L)
        vs.append(np.array([beam_field_interpolation(t, L, kind)[component] @ Ue for t in xi]))
        xs.append([np.nan])
        vs.append([np.nan])
    return np.concatenate(xs[:-1]), np.concatenate(vs[:-1])

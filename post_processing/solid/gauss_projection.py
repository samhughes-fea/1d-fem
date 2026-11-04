# post_processing/solid/gauss_projection.py
"""Fields over 2D section elements from their Gauss-point values.

A ``Q8`` section element is integrated on 3 x 3 Gauss points (``solid.elements.quadrature``),
so its stresses are first known there (``SolidSectionAnalysis.point_stresses``, element by
element in quadrature order). To show a stress over the element, those nine values are
projected with the nine-node Lagrange shape functions written in the Gauss points' own
coordinates, ``(xi, eta) / sqrt(3/5)``: the interpolant passes through every Gauss value and
extrapolates to the element edges and nodes as the standard Gauss-to-node projection does.
Nothing is smoothed or fitted, and elements are not averaged with their neighbours.

Plots then draw the projected field per element: ``element_patches`` splits each element into
sub-cells through its own isoparametric map and colours each by the projection at the
sub-cell centre; ``element_line`` samples it along a line of constant ``xi``.

Public API
----------
GAUSS_1D, gauss_lagrange_q8, project_q8, element_patches, element_line.
"""

from __future__ import annotations

from typing import List, Sequence, Tuple

import numpy as np
import numpy.typing as npt

from pre_processing.section_integration.solid.elements import quadrature, shape

__all__ = ["GAUSS_1D", "gauss_lagrange_q8", "project_q8", "element_patches", "element_line"]

GAUSS_1D = np.sqrt(0.6) * np.array([-1.0, 0.0, 1.0])


def _lagrange3(t: float) -> npt.NDArray[np.float64]:
    """1D quadratic Lagrange shapes through ``-1, 0, 1``."""
    return np.array([0.5 * t * (t - 1.0), 1.0 - t * t, 0.5 * t * (t + 1.0)])


def gauss_lagrange_q8(xi: float, eta: float) -> npt.NDArray[np.float64]:
    """Weights (9,) of the nine Gauss values of a ``Q8`` at ``(xi, eta)``, in quadrature order.

    ``quadrature("Q8")`` orders the points with ``xi`` outer and ``eta`` inner, so point
    ``3 i + j`` sits at ``(GAUSS_1D[i], GAUSS_1D[j])``.
    """
    a = _lagrange3(float(xi) / np.sqrt(0.6))
    b = _lagrange3(float(eta) / np.sqrt(0.6))
    return np.outer(a, b).ravel()


def project_q8(gauss_values: npt.ArrayLike, xi: float, eta: float) -> npt.NDArray[np.float64]:
    """The projected field at ``(xi, eta)`` from the element's nine Gauss values (9,) or (9, k)."""
    v = np.asarray(gauss_values, dtype=np.float64)
    if v.shape[0] != 9:
        raise ValueError("a Q8 element has nine Gauss values")
    return np.tensordot(gauss_lagrange_q8(xi, eta), v, axes=(0, 0))


def _check_order() -> None:
    pts, _ = quadrature("Q8")
    for k, (x, e) in enumerate(pts):
        if not (np.isclose(x, GAUSS_1D[k // 3]) and np.isclose(e, GAUSS_1D[k % 3])):
            raise RuntimeError("Q8 quadrature order changed; update gauss_lagrange_q8")


_check_order()


def element_patches(nodes: npt.ArrayLike, connectivity: Sequence[int], gauss_values: npt.ArrayLike,
                    n_sub: int = 4) -> Tuple[List[npt.NDArray[np.float64]], npt.NDArray[np.float64]]:
    """Sub-cell polygons of one ``Q8`` and the projected value at each sub-cell centre.

    The sub-cell corners come from the element's own map ``x(xi, eta) = N(xi, eta) X``, so the
    patches tile the curved element exactly as the element sees it.
    """
    X = np.asarray(nodes, dtype=np.float64)[list(connectivity)]
    t = np.linspace(-1.0, 1.0, n_sub + 1)
    c = 0.5 * (t[:-1] + t[1:])

    def at(x, e):
        return shape("Q8", x, e)[0] @ X

    polys, vals = [], []
    for i in range(n_sub):
        for j in range(n_sub):
            polys.append(np.array([at(t[i], t[j]), at(t[i + 1], t[j]), at(t[i + 1], t[j + 1]), at(t[i], t[j + 1])]))
            vals.append(project_q8(gauss_values, c[i], c[j]))
    return polys, np.asarray(vals)


def element_line(gauss_values: npt.ArrayLike, xi: float = 0.0, n: int = 21
                 ) -> Tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]:
    """``(eta, value)`` of the projected field along ``xi`` = const, ``eta`` from -1 to 1."""
    eta = np.linspace(-1.0, 1.0, n)
    return eta, np.array([project_q8(gauss_values, xi, e) for e in eta])

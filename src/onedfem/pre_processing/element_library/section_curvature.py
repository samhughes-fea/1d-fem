# pre_processing/element_library/section_curvature.py
"""Winkler curvature correction of a beam section.

On a curved axis the fibres inside the bend are shorter than those outside. A fibre
at ``(y, z)`` has length ``g = 1 - z k_y - y k_z`` per length of the axis, so the
axial and bending block of the section stiffness becomes ``∫ E a aᵀ / g dA`` with
``a = (1, -z, -y)``. ``with_curvature`` applies that correction to a provider
``SectionAnalysis`` (the class itself lives in the section provider; this is the
1d-fem side of it, next to ``SectionAnalysis.with_pretwist``).
"""

from __future__ import annotations

import dataclasses

import numpy as np


def with_curvature(section, k_y: float, k_z: float):
    """``section`` on a beam curved at ``(k_y, k_z)`` (1/m, element axes).

    The axial and bending block of every stiffness (the 6 × 6, classical, Benscoter,
    third-order and GBT ones) gains::

        D[i, j] += ∫ E a_i a_j (1/g - 1) dA,    a = (1, -z, -y),
        g = 1 - z k_y - y k_z,    i, j in (eps_x, kappa_y, kappa_z)

    This is the Winkler-Bach curved beam. To first order it couples extension and
    bending by ``-k EI`` and is otherwise a ``(k h)²`` change. The integrals use the
    section's ``axial_points``. An isotropic ``section.txt`` rectangle is stored on
    3 × 3 points (exact for the pretwist polynomials); those are rebuilt here on
    5 × 5 points, which is what ``1/g`` needs. Shear, torsion, warping and mass are
    not corrected. A curvature that puts the centre of curvature inside the section
    is refused. Zero curvature returns ``section`` itself.
    """
    if not k_y and not k_z:
        return section
    points = _points(section)
    y, z, wE = np.asarray(points, dtype=np.float64).T
    g = 1.0 - z * float(k_y) - y * float(k_z)
    if np.any(g <= 0.0):
        raise ValueError(
            f"curvature ({k_y:g}, {k_z:g}) 1/m puts the centre of curvature inside "
            "the section (1 - z k_y - y k_z <= 0 at a section point)"
        )
    a = np.stack([np.ones_like(y), -z, -y])
    dD = (a * (wE * (1.0 / g - 1.0))) @ a.T
    dD = 0.5 * (dD + dD.T)

    def add(K):
        if K is None:
            return None
        K = np.array(K, dtype=np.float64, copy=True)
        K[:3, :3] += dD
        return K

    return dataclasses.replace(
        section,
        stiffness=add(section.stiffness),
        classical_stiffness=add(section.classical_stiffness),
        benscoter_stiffness=add(section.benscoter_stiffness),
        third_order_stiffness=add(section.third_order_stiffness),
        gbt_stiffness=add(section.gbt_stiffness),
    )


def _points(section) -> np.ndarray:
    pts = section.axial_points
    if pts is None:
        raise ValueError(
            "the curvature correction needs the section's axial_points "
            "(a section with area and inertias)"
        )
    pts = np.asarray(pts, dtype=np.float64)
    if getattr(section, "source", None) == "isotropic" and pts.shape == (9, 3):
        return _rectangle_5x5(pts)
    return pts


def _rectangle_5x5(points: np.ndarray) -> np.ndarray:
    """5 × 5 Gauss points of the uniform rectangle implied by a 3 × 3 set.

    The weights are ``E dA``. The 3 × 3 rule is exact for ``∫ y²``, ``∫ z²`` and
    ``∫ 1``, so the extents and ``EA`` recovered here are the rectangle's own.
    """
    y, z, w = points.T
    wsum = float(w.sum())
    if wsum == 0.0:
        return points
    b = float(np.sqrt(max(12.0 * float(w @ (y * y)) / wsum, 0.0)))
    h = float(np.sqrt(max(12.0 * float(w @ (z * z)) / wsum, 0.0)))
    if b == 0.0 or h == 0.0:
        return points
    xi, gw = np.polynomial.legendre.leggauss(5)
    return np.array(
        [[0.5 * b * yi, 0.5 * h * zi, wsum * 0.25 * wy * wz]
         for yi, wy in zip(xi, gw) for zi, wz in zip(xi, gw)]
    )

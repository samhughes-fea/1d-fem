# pre_processing/element_library/element_frame.py
"""Local frames of beam elements: straight elements in any direction, and curved, twisted ones.

A beam element's own axes are ``x`` along the beam and ``y``, ``z`` in its section, the
axes its stiffness, section and stresses are written in. Before this module every linear
beam assumed those were the global axes, so a beam had to lie along global ``x``.

Straight elements
-----------------
``x`` runs from node 1 to node 2. With zero section twist, ``y`` is horizontal,
``y = Z cross x / |Z cross x|`` (so ``z`` points up, ``z . Z >= 0``), and for a vertical
element ``y = Y``. The section twist ``beta`` (the ``[twist]`` column of ``element.txt``,
in degrees) turns ``y`` and ``z`` about ``x`` by ``beta``, right handed. An element along
``+X`` with zero twist has its axes equal to the global axes, so existing jobs are
unchanged.

Curved and twisted elements
---------------------------
``precurvature.txt`` gives an element a constant reference curvature
``k = [k_x0, k_y0, k_z0]`` in its own axes: twist rate and the two bending curvatures.
The element is then the arc of that curvature between its nodes, its frame along the
arc ``R(s) = R_mid exp(s [k]x)`` for arc length ``s`` measured from the middle. ``R_mid``
is chosen so that the arc's chord runs from node 1 to node 2 and the twist ``beta`` sets
the section's turn about the chord, as for a straight element. The arc length follows
from the chord length.

Public API
----------
rotation_exp, straight_frame, arc_geometry, ArcGeometry, node_dof_rotation.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import numpy as np
import numpy.typing as npt

_Z = np.array([0.0, 0.0, 1.0])
_Y = np.array([0.0, 1.0, 0.0])


def _skew(a: npt.ArrayLike) -> np.ndarray:
    x, y, z = np.asarray(a, dtype=np.float64)
    return np.array([[0.0, -z, y], [z, 0.0, -x], [-y, x, 0.0]])


def rotation_exp(a: npt.ArrayLike) -> np.ndarray:
    """Rotation matrix of the rotation vector ``a`` (Rodrigues)."""
    a = np.asarray(a, dtype=np.float64)
    t = float(np.linalg.norm(a))
    K = _skew(a)
    if t < 1e-8:
        return np.eye(3) + K + 0.5 * K @ K
    return np.eye(3) + np.sin(t) / t * K + (1.0 - np.cos(t)) / t**2 * K @ K


def straight_frame(p0: npt.ArrayLike, p1: npt.ArrayLike, twist_deg: float = 0.0) -> np.ndarray:
    """Columns: the element's ``x``, ``y``, ``z`` in global axes (see the module notes)."""
    d = np.asarray(p1, dtype=np.float64) - np.asarray(p0, dtype=np.float64)
    L = float(np.linalg.norm(d))
    if L <= 0.0:
        raise ValueError("element has zero length")
    ex = d / L
    y = np.cross(_Z, ex)
    if np.linalg.norm(y) < 1e-9:          # vertical element
        y = _Y - (_Y @ ex) * ex
    ey = y / np.linalg.norm(y)
    ez = np.cross(ex, ey)
    R = np.column_stack([ex, ey, ez])
    if twist_deg:
        R = R @ rotation_exp([np.deg2rad(twist_deg), 0.0, 0.0])
    return R


def _chord_integral(k: np.ndarray, length: float, n: int = 64) -> np.ndarray:
    """``int_{-l/2}^{l/2} exp(s [k]x) e_x ds`` in the frame at the middle of the arc."""
    x, w = np.polynomial.legendre.leggauss(n)
    s = 0.5 * length * x
    return 0.5 * length * sum(wi * rotation_exp(si * k)[:, 0] for si, wi in zip(s, w))


def _align(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """The smallest rotation taking unit vector ``a`` to unit vector ``b``."""
    v = np.cross(a, b)
    s, c = float(np.linalg.norm(v)), float(a @ b)
    if s < 1e-14:
        return np.eye(3)
    return rotation_exp(v / s * np.arctan2(s, c))


@dataclass(frozen=True)
class ArcGeometry:
    """An element arc: its frames and length, and the straight sub-chords that approximate it."""

    R_mid: np.ndarray            # (3, 3) frame at the middle of the arc
    curvature: np.ndarray        # (3,) k in the element's own axes
    length: float                # arc length
    chord: float                 # chord length

    def frame_at(self, s: float) -> np.ndarray:
        """Frame at arc length ``s`` from the middle."""
        return self.R_mid @ rotation_exp(s * self.curvature)

    @property
    def R0(self) -> np.ndarray:
        return self.frame_at(-0.5 * self.length)

    @property
    def R1(self) -> np.ndarray:
        return self.frame_at(0.5 * self.length)

    def position(self, s: float, p_mid: np.ndarray) -> np.ndarray:
        """Point at arc length ``s`` from the middle, given the middle point."""
        if abs(s) < 1e-15:
            return p_mid
        x, w = np.polynomial.legendre.leggauss(32)
        t = 0.5 * s * (x + 1.0)
        return p_mid + self.R_mid @ (0.5 * s * sum(wi * rotation_exp(ti * self.curvature)[:, 0]
                                                   for ti, wi in zip(t, w)))

    def sub_frames(self, m: int) -> list:
        """Frames of ``m`` straight sub-elements along the arc, each with ``x`` along its own chord."""
        h = self.length / m
        c = _chord_integral(self.curvature, h)
        S = _align(c / np.linalg.norm(c), np.array([1.0, 0.0, 0.0]))
        return [self.frame_at(-0.5 * self.length + (i + 0.5) * h) @ S.T for i in range(m)]

    def sub_chord(self, m: int) -> float:
        return float(np.linalg.norm(_chord_integral(self.curvature, self.length / m)))


def arc_geometry(p0: npt.ArrayLike, p1: npt.ArrayLike, curvature: npt.ArrayLike,
                 twist_deg: float = 0.0) -> ArcGeometry:
    """The arc of constant ``curvature`` (element axes) from ``p0`` to ``p1``."""
    k = np.asarray(curvature, dtype=np.float64).reshape(3)
    p0, p1 = np.asarray(p0, dtype=np.float64), np.asarray(p1, dtype=np.float64)
    chord = float(np.linalg.norm(p1 - p0))
    R_chord = straight_frame(p0, p1, twist_deg)
    length = chord
    for _ in range(100):                      # chord(l) = chord: fixed point on the length
        c = _chord_integral(k, length)
        new = length * chord / float(np.linalg.norm(c))
        if abs(new - length) <= 1e-14 * chord:
            length = new
            break
        length = new
    c = _chord_integral(k, length)
    # R_mid c / |c| must be the chord direction, with the section turned about the chord as R_chord says
    S = _align(c / np.linalg.norm(c), np.array([1.0, 0.0, 0.0]))
    return ArcGeometry(R_chord @ S, k, float(length), chord)


def node_dof_rotation(R0: np.ndarray, R1: np.ndarray, n_dof: int, blocks: Optional[list] = None,
                      pairs: Optional[list] = None) -> np.ndarray:
    """``T`` with ``u_local = T @ u_global`` for an element whose nodes carry frames ``R0``, ``R1``.

    ``blocks`` lists ``(start, node)``: three DOFs from ``start`` are a vector at that node
    (default: the translations and rotations of the standard 12, ``[0, 3]`` at node 0 and
    ``[6, 9]`` at node 1). ``pairs`` lists ``(start, node, N)``: two DOFs from ``start``
    are the cross-section (``y``, ``z``) components of a vector held at the node in the
    frame ``N`` (the node's own frame), mapped to the element's by
    ``[e_y, e_z]^T [n_y, n_z]``: a plane rotation when the frames share their ``x``.
    Every other DOF (warping, mode amplitudes) is a scalar and maps to itself.
    """
    T = np.eye(n_dof)
    for start, node in (blocks or [(0, 0), (3, 0), (6, 1), (9, 1)]):
        T[start:start + 3, start:start + 3] = (R0 if node == 0 else R1).T
    for start, node, N in (pairs or ()):
        R = R0 if node == 0 else R1
        T[start:start + 2, start:start + 2] = R[:, 1:3].T @ np.asarray(N)[:, 1:3]
    return T

# pre_processing/element_library/nonlinear/large_rotations/corotational/corotational_section.py
"""Corotational beam elements on any SectionAnalysis, with third- and fourth-order local energy.

A linear element (Euler-Bernoulli, Timoshenko or Benscoter, with its coupled section
stiffness) is the local core. The corotational frame takes out the rigid motion of each
element, and the element is left with small local displacements ``p``:

- ``p`` holds the stretch ``l - l0`` of the chord and the two nodal rotations relative to
  the frame, ``theta_i = log(R_r^T R_i)``. Extra DOFs (the Benscoter warping amplitude) are
  scalars of the section and pass through unchanged.
- The frame ``R_r`` has its first axis along the current chord. Its second axis is the
  mean of the nodal second axes made orthogonal to the chord (Battini and Pacoste 2002).

Nodal rotations are total rotation vectors, ``R_i = exp(theta_i) E_0``, so the
additive Newton update of the runner is consistent.

Local energy, in the core element's DOFs::

    Pi(p) = 1/2 p^T K p + 1/2 p^T K_sigma(p) p + Pi_4(p)

``K`` is the core's linear stiffness and ``K_sigma`` its second-order geometric stiffness
(``docs/element_library/second_order.md``), which is linear in ``p``. The cubic term is
the work of the stresses on the quadratic strains of the rigid in-plane section motion.
It carries the shear-centre and tension-centre offsets, the Wagner term and the bimoment.
``Pi_4 = 1/2 int sum_k k (eps_2)^2`` with
``eps_2 = 1/2 [(v' - z phi')^2 + (w' + y phi')^2]`` over the section's axial points
completes the axial strain of large twist (Wagner torsional stiffening).

Global force ``F = B^T f(p)`` with ``B = dp/dU``. The tangent is
``B^T H(p) B + sum_i f_i d2p_i/dU2``. ``B`` comes from the complex step on the closed-form
map ``p(U)``, exact to round-off, so the residual is exact. ``d2p/dU2`` comes from central
differences, and ``f`` and ``H`` are exact.

Public API
----------
CorotationalEulerBernoulliBeamElement3D, CorotationalTimoshenkoBeamElement3D,
CorotationalBenscoterBeamElement3D, rotation_exp, rotation_log.
"""

from __future__ import annotations

from typing import Tuple

import numpy as np

from pre_processing.element_library.linear.beam.first_order_shear_deformation_theory.benscoter.linear_benscoter_3D import (
    LinearBenscoterBeamElement3D,
)
from pre_processing.element_library.linear.beam.first_order_shear_deformation_theory.timoshenko.linear_timoshenko_3D import (
    LinearTimoshenkoBeamElement3D,
)
from pre_processing.element_library.linear.beam.zero_order_shear_deformation_theory.euler_bernoulli.linear_euler_bernoulli_3D import (
    LinearEulerBernoulliBeamElement3D,
)

_GAUSS_X, _GAUSS_W = np.polynomial.legendre.leggauss(4)
#: central-difference step for d2p/dU2 (rotations; translations scale with l0)
_H2 = 1e-4


def _skew(a: np.ndarray) -> np.ndarray:
    return np.array([[0.0, -a[2], a[1]], [a[2], 0.0, -a[0]], [-a[1], a[0], 0.0]])


def rotation_exp(theta: np.ndarray) -> np.ndarray:
    """Rotation matrix of the rotation vector ``theta`` (Rodrigues)."""
    t = float(np.linalg.norm(theta))
    S = _skew(theta)
    if t < 1e-8:
        return np.eye(3) + S + 0.5 * S @ S
    return np.eye(3) + np.sin(t) / t * S + (1.0 - np.cos(t)) / (t * t) * S @ S


def rotation_log(R: np.ndarray) -> np.ndarray:
    """Rotation vector of ``R`` (through the quaternion, robust up to pi)."""
    tr = np.trace(R)
    w2 = np.array([1.0 + tr, 1.0 + 2 * R[0, 0] - tr, 1.0 + 2 * R[1, 1] - tr, 1.0 + 2 * R[2, 2] - tr])
    k = int(np.argmax(w2))
    q = np.zeros(4)
    q[k] = 0.5 * np.sqrt(w2[k])
    if k == 0:
        q[1] = (R[2, 1] - R[1, 2]) / (4 * q[0])
        q[2] = (R[0, 2] - R[2, 0]) / (4 * q[0])
        q[3] = (R[1, 0] - R[0, 1]) / (4 * q[0])
    elif k == 1:
        q[0] = (R[2, 1] - R[1, 2]) / (4 * q[1])
        q[2] = (R[0, 1] + R[1, 0]) / (4 * q[1])
        q[3] = (R[0, 2] + R[2, 0]) / (4 * q[1])
    elif k == 2:
        q[0] = (R[0, 2] - R[2, 0]) / (4 * q[2])
        q[1] = (R[0, 1] + R[1, 0]) / (4 * q[2])
        q[3] = (R[1, 2] + R[2, 1]) / (4 * q[2])
    else:
        q[0] = (R[1, 0] - R[0, 1]) / (4 * q[3])
        q[1] = (R[0, 2] + R[2, 0]) / (4 * q[3])
        q[2] = (R[1, 2] + R[2, 1]) / (4 * q[3])
    if q[0] < 0:
        q = -q
    v = q[1:]
    s = float(np.linalg.norm(v))
    if s < 1e-12:
        return 2.0 * v
    return 2.0 * np.arctan2(s, q[0]) / s * v


def _norm(a: np.ndarray) -> np.ndarray:
    """``sqrt(a . a)`` along the last axis, analytic (complex-step safe)."""
    return np.sqrt(np.sum(a * a, axis=-1))


def _exp_batch(th: np.ndarray) -> np.ndarray:
    """Rotation matrices (m, 3, 3) of rotation vectors (m, 3); complex-step safe."""
    t2 = np.sum(th * th, axis=1)
    small = np.abs(t2) < 1e-12
    t = np.sqrt(np.where(small, 1.0, t2))
    a = np.where(small, 1.0 - t2 / 6.0, np.sin(t) / t)
    b = np.where(small, 0.5 - t2 / 24.0, (1.0 - np.cos(t)) / (t * t))
    S = np.zeros((th.shape[0], 3, 3), dtype=th.dtype)
    S[:, 0, 1], S[:, 0, 2], S[:, 1, 2] = -th[:, 2], th[:, 1], -th[:, 0]
    S[:, 1, 0], S[:, 2, 0], S[:, 2, 1] = th[:, 2], -th[:, 1], th[:, 0]
    return np.eye(3) + a[:, None, None] * S + b[:, None, None] * (S @ S)


def _log_batch(R: np.ndarray) -> np.ndarray:
    """Rotation vectors (m, 3) of rotations (m, 3, 3) below pi/2 in angle (the local rotations).

    Uses ``t = arctan(s / c)``, analytic for the complex step.
    """
    v = 0.5 * np.stack([R[:, 2, 1] - R[:, 1, 2], R[:, 0, 2] - R[:, 2, 0], R[:, 1, 0] - R[:, 0, 1]], axis=1)
    s2 = np.sum(v * v, axis=1)
    c = 0.5 * (np.trace(R, axis1=1, axis2=2) - 1.0)
    small = np.abs(s2) < 1e-24
    s = np.sqrt(np.where(small, 1.0, s2))
    f = np.where(small, 1.0 / c + s2 / 6.0, np.arctan(s / c) / s)
    return f[:, None] * v


def _reference_frame(d: np.ndarray) -> np.ndarray:
    """Undeformed element frame: ``e1`` along the chord, ``e2`` towards global ``y`` (or ``z``)."""
    e1 = d / np.linalg.norm(d)
    ref = np.array([0.0, 1.0, 0.0]) if abs(e1[1]) < 0.9 else np.array([0.0, 0.0, 1.0])
    e2 = ref - (ref @ e1) * e1
    e2 /= np.linalg.norm(e2)
    return np.column_stack([e1, e2, np.cross(e1, e2)])


class CorotationalMixin:
    """Corotational large-displacement behaviour on top of a linear element core."""

    supports_element_frame = False   # rotated or curved elements: not yet

    supports_section_analysis = True

    # ------------------------------------------------------------------ #
    def _cr_setup(self) -> None:
        if getattr(self, "_cr_ready", False):
            return
        if np.any(self._E_0_voigt):
            raise NotImplementedError(f"Element {self.element_id}: precurvature is not supported by the corotational elements yet")
        X = np.asarray(self.node_coords, dtype=np.float64)
        self._X1, self._X2 = X[0].copy(), X[1].copy()
        self._E0 = _reference_frame(self._X2 - self._X1)
        self._l0 = float(np.linalg.norm(self._X2 - self._X1))
        n = self._n_dof
        logger, self.logger_operator = self.logger_operator, None
        try:
            self._K_core = np.asarray(self.element_stiffness_matrix().K_e, dtype=np.float64)
        finally:
            self.logger_operator = logger
        eye = np.eye(n)
        self._C = np.array([self.linear_geometric_stiffness_matrix(eye[j]) for j in range(n)])
        pts = self.section_analysis.axial_points
        self._axial_pts = None if pts is None else np.asarray(pts, dtype=np.float64)
        self._G_gauss = [self._geometric_state(float(x), np.zeros(n))[0] for x in _GAUSS_X]
        self._cr_ready = True

    # ------------------------------------------------------------------ #
    # Kinematics
    # ------------------------------------------------------------------ #
    def _local_dofs(self, U_e: np.ndarray) -> np.ndarray:
        """Local deformational DOFs ``p`` (core element order) of the global state ``U_e``."""
        return self._local_dofs_batch(np.asarray(U_e, dtype=np.float64)[None, :])[0]

    def _local_dofs_batch(self, U: np.ndarray) -> np.ndarray:
        """``p`` (m, n) for global states ``U`` (m, n)."""
        d = (self._X2 + U[:, 6:9]) - (self._X1 + U[:, 0:3])
        l = _norm(d)
        if np.any(np.abs(l) <= 1e-12 * self._l0):
            raise ValueError(f"Element {self.element_id}: corotational chord collapsed")
        R1 = _exp_batch(U[:, 3:6]) @ self._E0
        R2 = _exp_batch(U[:, 9:12]) @ self._E0
        e1 = d / l[:, None]
        q = 0.5 * (R1[:, :, 1] + R2[:, :, 1])
        e3 = np.cross(e1, q)
        e3 = e3 / _norm(e3)[:, None]
        Rr = np.stack([e1, np.cross(e3, e1), e3], axis=2)
        RrT = np.transpose(Rr, (0, 2, 1))
        p = np.zeros(U.shape, dtype=U.dtype)
        p[:, 3:6] = _log_batch(RrT @ R1)
        p[:, 6] = l - self._l0
        p[:, 9:12] = _log_batch(RrT @ R2)
        p[:, 12:] = U[:, 12:]
        return p

    def _steps(self, h: float) -> np.ndarray:
        s = np.full(12, h)
        s[[0, 1, 2, 6, 7, 8]] *= self._l0
        return s

    def _local_jacobian(self, U_e: np.ndarray) -> np.ndarray:
        """``B = dp/dU`` by the complex step (exact to round-off)."""
        n = self._n_dof
        h = 1e-30
        Z = np.tile(U_e.astype(np.complex128), (12, 1))
        Z[np.arange(12), np.arange(12)] += 1j * h
        B = np.zeros((n, n))
        B[:, :12] = (self._local_dofs_batch(Z).imag / h).T
        B[12:, 12:] = np.eye(n - 12)
        return B

    def _local_hessian_contraction(self, U_e: np.ndarray, f: np.ndarray) -> np.ndarray:
        """``sum_i f_i d2p_i/dU2`` on the twelve standard DOFs (``p`` is linear in the extras)."""
        n = self._n_dof
        hs = self._steps(_H2)
        j, k = np.triu_indices(12)
        m = j.size
        D = np.zeros((4, m, n))
        for c, (sj, sk) in enumerate(((1, 1), (1, -1), (-1, 1), (-1, -1))):
            D[c, np.arange(m), j] += sj * hs[j]
            D[c, np.arange(m), k] += sk * hs[k]
        vals = (self._local_dofs_batch((U_e[None, :] + D.reshape(-1, n))) @ f).reshape(4, m)
        h2 = (vals[0] - vals[1] - vals[2] + vals[3]) / (4 * hs[j] * hs[k])
        Kg = np.zeros((n, n))
        Kg[j, k] = h2
        Kg[k, j] = h2
        return Kg

    # ------------------------------------------------------------------ #
    # Local energy
    # ------------------------------------------------------------------ #
    def _local_force_and_tangent(self, p: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        C = self._C
        Ks = np.tensordot(p, C, axes=(0, 0))       # K_sigma(p) = sum_j p_j C_j
        M = C @ p                                   # M[j, k] = (C_j p)_k
        f = self._K_core @ p + Ks @ p + 0.5 * np.einsum("jkl,k,l->j", C, p, p)
        H = self._K_core + Ks + M + M.T
        pts = self._axial_pts
        if pts is not None:
            y, z, k = pts[:, 0], pts[:, 1], pts[:, 2]
            detJ = 0.5 * self._l0
            for G, wg in zip(self._G_gauss, _GAUSS_W):
                v1, w1, f1, _ = G @ p
                a = v1 - z * f1
                b = w1 + y * f1
                e2 = 0.5 * (a * a + b * b)
                Ja = np.stack([a, b, -z * a + y * b], axis=1)    # d eps_2 / d[v', w', phi']
                grad = (k * e2) @ Ja
                Hg = Ja.T @ (k[:, None] * Ja)
                ke = k * e2
                Hg += np.array([[ke.sum(), 0.0, -(ke * z).sum()],
                                [0.0, ke.sum(), (ke * y).sum()],
                                [-(ke * z).sum(), (ke * y).sum(), (ke * (y * y + z * z)).sum()]])
                G3 = G[:3]
                f = f + wg * detJ * G3.T @ grad
                H = H + wg * detJ * G3.T @ Hg @ G3
        return f, H

    # ------------------------------------------------------------------ #
    # Nonlinear element interface
    # ------------------------------------------------------------------ #
    def local_deformation(self, U_e: np.ndarray) -> np.ndarray:
        """Local deformational DOFs ``p`` (for results and tests)."""
        self._cr_setup()
        return self._local_dofs(np.asarray(U_e, dtype=np.float64).ravel())

    def internal_force_vector(self, U_e: np.ndarray) -> np.ndarray:
        self._cr_setup()
        U = np.asarray(U_e, dtype=np.float64).ravel()
        p = self._local_dofs(U)
        f, _ = self._local_force_and_tangent(p)
        return self._local_jacobian(U).T @ f

    def tangent_stiffness_matrix(self, U_e: np.ndarray) -> np.ndarray:
        self._cr_setup()
        U = np.asarray(U_e, dtype=np.float64).ravel()
        p = self._local_dofs(U)
        f, H = self._local_force_and_tangent(p)
        B = self._local_jacobian(U)
        K = B.T @ H @ B + self._local_hessian_contraction(U, f)
        return 0.5 * (K + K.T) if np.allclose(self._K_core, self._K_core.T) else K


class CorotationalEulerBernoulliBeamElement3D(CorotationalMixin, LinearEulerBernoulliBeamElement3D):
    """Corotational Euler-Bernoulli beam (shear-rigid core, classical section stiffness)."""

    element_type_name = "Corotational-EB-3D"


class CorotationalTimoshenkoBeamElement3D(CorotationalMixin, LinearTimoshenkoBeamElement3D):
    """Corotational Timoshenko beam on the coupled Saint-Venant section stiffness."""

    element_type_name = "Corotational-Timoshenko-3D"


class CorotationalBenscoterBeamElement3D(CorotationalMixin, LinearBenscoterBeamElement3D):
    """Corotational Benscoter beam: restrained warping and large twist on thin-walled sections."""

    element_type_name = "Corotational-Benscoter-3D"

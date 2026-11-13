# pre_processing/element_library/nonlinear/large_rotations/corotational/corotational_section.py
"""Corotational beam elements on any SectionAnalysis, with third- and fourth-order local energy.

A linear element (Euler-Bernoulli, Timoshenko or Benscoter, with its coupled section
stiffness) is the local core. The corotational frame takes out the rigid motion of each
element, and the element is left with small local displacements ``p``:

- ``p`` holds the stretch ``l - l0`` of the chord and the two nodal rotations relative to
  the frame, ``theta_i = log(R_r^T R_i)``. Extra DOFs that are scalars of the section (the
  Benscoter warping amplitude, the GBT warping and distortional mode amplitudes) pass
  through unchanged. The third-order beams (Reddy, Levinson) hold the section shear
  strains ``[gamma_xy, gamma_xz]`` at the nodes, in the node's own rotating frame, and
  their local slopes are ``v' = gamma_xy + theta_z`` and ``w' = gamma_xz - theta_y`` with
  the local nodal rotations (``_local_extras``).
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
CorotationalBenscoterBeamElement3D, CorotationalGBTBeamElement3D,
CorotationalReddyBeamElement3D, CorotationalLevinsonBeamElement3D, rotation_exp, rotation_log.
"""

from __future__ import annotations

from typing import Tuple

import numpy as np

from onedfem.pre_processing.element_library.linear.beam.fosdt.benscoter.linear_benscoter_3D import (
    LinearBenscoterBeamElement3D,
)
from onedfem.pre_processing.element_library.linear.beam.fosdt.gbt.linear_gbt_3D import (
    LinearGBTBeamElement3D,
)
from onedfem.pre_processing.element_library.linear.beam.second_order import GeometricStiffnessMixin, _linear_rows
from onedfem.pre_processing.element_library.linear.beam.tosdt.levinson.linear_levinson_3D import (
    LinearLevinsonBeamElement3D,
)
from onedfem.pre_processing.element_library.linear.beam.tosdt.reddy.linear_reddy_3D import (
    LinearReddyBeamElement3D,
)
from onedfem.pre_processing.element_library.linear.beam.fosdt.timoshenko.linear_timoshenko_3D import (
    LinearTimoshenkoBeamElement3D,
)
from onedfem.pre_processing.element_library.linear.beam.zosdt.euler_bernoulli.linear_euler_bernoulli_3D import (
    LinearEulerBernoulliBeamElement3D,
)
from onedfem.pre_processing.element_library.nonlinear.large_rotations.corotational.utilities.frame import (
    reference_frame as _reference_frame,
)
from onedfem.pre_processing.element_library.nonlinear.large_rotations.corotational.utilities.local_energy import (
    cubic_quartic_force_and_tangent,
)
from onedfem.pre_processing.element_library.nonlinear.large_rotations.corotational.utilities.local_map import (
    battini_p_batch,
)
from onedfem.pre_processing.element_library.nonlinear.large_rotations.corotational.utilities.rotation import (
    rotation_exp,
    rotation_log,
)

_GAUSS_X, _GAUSS_W = np.polynomial.legendre.leggauss(4)
#: central-difference step for d2p/dU2 (rotations; translations scale with l0)
_H2 = 1e-4


class CorotationalMixin:
    """Corotational large-displacement behaviour on top of a linear element core."""

    supports_element_frame = True

    supports_section_analysis = True

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._apply_cr_reference_curvature()

    def _sync_pretwisted_section_operators(self) -> None:
        """Rebuild cores that snapshot ``D`` in ``__init__`` before pretwist is applied.

        Timoshenko reads ``section_analysis.stiffness`` at assembly time. Benscoter, GBT
        and tosdt copy ``D`` during ``super().__init__``, which runs before this mixin
        sets ``_pretwist_rate`` from ``k_x0``.
        """
        sa = self.section_analysis
        if getattr(self, "benscoter_material_stiffness_operator", None) is not None:
            if sa.benscoter_stiffness is not None:
                from onedfem.pre_processing.element_library.linear.beam.fosdt.benscoter.utilities.D_matrix import (
                    MaterialStiffnessOperator as BenscoterD,
                )

                self.benscoter_material_stiffness_operator = BenscoterD(
                    D=np.array(sa.benscoter_stiffness, dtype=np.float64)
                )
        if getattr(self, "gbt_material_stiffness_operator", None) is not None:
            from onedfem.pre_processing.element_library.linear.beam.fosdt.gbt.utilities.D_matrix import (
                MaterialStiffnessOperator as GBTD,
            )

            K, M = sa.gbt_matrices(self.n_modes, self.n_shear)
            self.gbt_material_stiffness_operator = GBTD(D=K, M=M)
            self._D_gbt, self._M_gbt = K, M
        if getattr(self, "_D8", None) is not None and sa.third_order_stiffness is not None:
            from onedfem.pre_processing.element_library.linear.beam.tosdt.utilities.D_matrix import (
                MaterialStiffnessOperator as ThirdOrderD,
            )

            self._D8 = np.asarray(sa.third_order_stiffness, dtype=np.float64)
            self.third_order_material_stiffness_operator = ThirdOrderD(D=self._D8)

    def _apply_cr_reference_curvature(self) -> None:
        """Pretwist ``k_x0`` couples ``D``; ``k_y0``, ``k_z0`` are constitutive ``E_0`` on the chord.

        CR stays geometrically straight (kinks carry discrete curvature). A fictitious
        arc from ``precurvature.txt`` is collapsed: pretwist-only keeps the mid-length
        section frame; bending precurvature uses the nodal chord frame.
        """
        E = np.array(self._E_0_voigt, dtype=np.float64, copy=True)
        if E.size < 6:
            E = np.zeros(6, dtype=np.float64)
        E[5] = 0.0
        arc = getattr(self, "_frame_arc", None)
        if arc is not None:
            k = np.asarray(arc.curvature, dtype=np.float64).reshape(3)
            kx, ky, kz = (float(v) for v in k)
            self._frame_arc = None
            self._frame_ctor = None
            if abs(ky) > 0.0 or abs(kz) > 0.0:
                from onedfem.pre_processing.element_library.element_frame import _align, _chord_integral

                c = _chord_integral(k, arc.length)
                S = _align(c / np.linalg.norm(c), np.array([1.0, 0.0, 0.0]))
                R_chord = np.asarray(arc.R_mid, dtype=np.float64) @ S.T
                self._frame_R0 = self._frame_R1 = R_chord
            else:
                self._frame_R0 = self._frame_R1 = arc.R_mid
            if abs(kx) > 0.0:
                self._pretwist_rate = kx
                self.__dict__.pop("_pretwisted_section_analysis", None)
                self._sync_pretwisted_section_operators()
            E[1], E[2] = ky, kz
        self._E_0_voigt = E

    # ------------------------------------------------------------------ #
    def _cr_setup(self) -> None:
        if getattr(self, "_cr_ready", False):
            return
        X = np.asarray(self.node_coords, dtype=np.float64)
        if getattr(self, "_frame_mesh", False):
            # framed mesh: the kernels work in the element's axes; the frame gives the geometry
            l0 = float(X[1, 0])
            self._E0 = np.asarray(self._frame_R0, dtype=np.float64).copy()
            self._X1 = np.asarray(self._frame_p0, dtype=np.float64).copy()
            self._X2 = self._X1 + l0 * self._E0[:, 0]
        else:
            self._X1, self._X2 = X[0].copy(), X[1].copy()
            self._E0 = _reference_frame(self._X2 - self._X1)
        self._l0 = float(np.linalg.norm(self._X2 - self._X1))
        n = self._n_dof
        logger, self.logger_operator = self.logger_operator, None
        depth = self.__dict__.get("_frame_depth", 0)
        self._frame_depth = 1        # the core's own (local) matrices, not turned to global
        try:
            self._K_core = np.asarray(self.element_stiffness_matrix().K_e, dtype=np.float64)
            E_saved = np.array(self._E_0_voigt, dtype=np.float64, copy=True)
            self._E_0_voigt = np.zeros(6, dtype=np.float64)
            try:
                eye = np.eye(n)
                self._C = np.array([self.linear_geometric_stiffness_matrix(eye[j]) for j in range(n)])
            finally:
                self._E_0_voigt = E_saved
            z = np.zeros(n, dtype=np.float64)
            self._K_sigma0 = np.asarray(self.linear_geometric_stiffness_matrix(z), dtype=np.float64)
            self._f_E0 = np.asarray(self._precurvature_equivalent_load(), dtype=np.float64)
        finally:
            self.logger_operator = logger
            self._frame_depth = depth
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
        return battini_p_batch(
            U,
            self._X1,
            self._X2,
            self._E0,
            self._l0,
            self._local_extras,
            self.element_id,
        )

    def _local_extras(self, U: np.ndarray, p: np.ndarray) -> np.ndarray:
        """Core-element extra DOFs from the global ones (scalars of the section: unchanged)."""
        return U[:, 12:]

    def _steps(self, h: float) -> np.ndarray:
        s = np.full(12, h)
        s[[0, 1, 2, 6, 7, 8]] *= self._l0
        return s

    def _local_jacobian(self, U_e: np.ndarray) -> np.ndarray:
        """``B = dp/dU`` by the complex step (exact to round-off)."""
        n = self._n_dof
        h = 1e-30
        Z = np.tile(U_e.astype(np.complex128), (n, 1))
        Z[np.arange(n), np.arange(n)] += 1j * h
        return (self._local_dofs_batch(Z).imag / h).T

    def _local_hessian_contraction(self, U_e: np.ndarray, f: np.ndarray) -> np.ndarray:
        """``sum_i f_i d2p_i/dU2`` on the twelve standard DOFs (``p`` is linear in the extras,
        with coefficients that do not depend on the standard DOFs)."""
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
        return cubic_quartic_force_and_tangent(
            p,
            self._K_core,
            self._C,
            self._axial_pts,
            self._G_gauss,
            _GAUSS_W,
            self._l0,
            f_E0=self._f_E0,
            K_sigma0=self._K_sigma0,
        )

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

    reduced_isotropic_section = True
    element_type_name = "Corotational-EB-3D"


class CorotationalTimoshenkoBeamElement3D(CorotationalMixin, LinearTimoshenkoBeamElement3D):
    """Corotational Timoshenko beam on the coupled Saint-Venant section stiffness."""

    reduced_isotropic_section = True
    element_type_name = "Corotational-Timoshenko-3D"


class CorotationalBenscoterBeamElement3D(CorotationalMixin, LinearBenscoterBeamElement3D):
    """Corotational Benscoter beam: restrained warping and large twist on thin-walled sections."""

    element_type_name = "Corotational-Benscoter-3D"


class CorotationalGBTBeamElement3D(CorotationalMixin, LinearGBTBeamElement3D):
    """Corotational GBT beam: large rotations of the section, with its warping and
    distortional modes linear in the element's own frame.

    The second-order terms are those of the global section forces (``_geometric_state``):
    the corotational frame carries the large rotations, and the cubic and quartic terms
    the Wagner, offset and stretch effects of moderate local rotation. The distortional
    modes have no wall-level geometric stiffness, so local and distortional buckling are
    not captured here (``LinearGBTBeamElement3D.linear_geometric_stiffness_matrix`` has them).

    **Brazier flattening.** A mode moves the walls in the section plane, so the beam
    curvature strains the moved fibres: ``e_B = sum_k a_k g_k . kappa``
    (``SectionAnalysis.gbt_brazier_terms``). Its work against the wall stresses of the GBT
    strains, ``sum_k a_k kappa^T T_k E``, and its own energy
    ``1/2 sum_kl a_k a_l kappa^T Q_kl kappa`` are added to the local energy at the Gauss
    points: the compressed flange of a bent box bows towards the neutral axis and the
    bending stiffness falls (Brazier 1927). Both are of third and fourth order in the local
    DOFs, as the other corotational terms.
    """

    element_type_name = "Corotational-GBT-3D"

    linear_geometric_stiffness_matrix = GeometricStiffnessMixin.linear_geometric_stiffness_matrix

    def _cr_setup(self) -> None:
        if getattr(self, "_cr_ready", False):
            return
        super()._cr_setup()
        bz = self.section_analysis.gbt_brazier_terms(self.n_modes, self.n_shear) if self.n_modes else None
        self._brazier = None
        if bz is not None:
            m, p = self.n_modes, self.n_warping
            r = self._reduced_rows
            rows = []
            for x in _GAUSS_X:
                B = self._B_at(float(x))[0].copy()
                B[r] = self._B_at(0.0)[0][r]
                Na = self._field_interpolation(float(x))[6 + p:6 + p + m]
                rows.append(np.vstack([Na, B[1:3], B]))
            self._brazier = (bz[0], bz[1], rows)

    def _local_force_and_tangent(self, p: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        f, H = super()._local_force_and_tangent(p)
        if self._brazier is None:
            return f, H
        T, Q, rows = self._brazier
        m = T.shape[0]
        nE = T.shape[2]
        detJ = 0.5 * self._l0
        for Z, wg in zip(rows, _GAUSS_W):
            g, Hz = brazier_energy_derivatives(T, Q, Z @ p)
            f = f + wg * detJ * Z.T @ g
            H = H + wg * detJ * Z.T @ Hz @ Z
        return f, H

    def _geometric_state(self, xi: float, U_e: np.ndarray):
        """Timoshenko rows of ``[v', w', phi', phi]`` and the terms of the GBT section forces.

        The section forces are the first six rows of ``D (B U_e - E_0)``: the resultants on
        the plane-section strains, with the warping amplitudes' part. The shear strains are
        taken at the centre, as the reduced shear rule sees them.
        """
        U_e = np.asarray(U_e, dtype=np.float64).ravel()
        e = self._B_at(float(xi))[0] @ U_e
        r = self._reduced_rows
        e[r] = (self._B_at(0.0)[0] @ U_e)[r]
        e[:6] -= self._E_0_voigt
        F = self._D_gbt[:6] @ e
        return _linear_rows(self.L, float(xi), self._n_dof), self._geometric_terms_matrix() @ F


def brazier_energy_derivatives(T: np.ndarray, Q: np.ndarray, z: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Gradient and Hessian of ``W = a^T (kappa . T E) + 1/2 a a kappa kappa . Q`` in ``z = [a, kappa, E]``.

    ``T`` is ``(m, 2, n)`` and ``Q`` ``(m, m, 2, 2)`` (``SectionAnalysis.gbt_brazier_terms``).
    """
    m, n = T.shape[0], T.shape[2]
    a, k, E = z[:m], z[m:m + 2], z[m + 2:]
    ia, ik, ie = slice(0, m), slice(m, m + 2), slice(m + 2, m + 2 + n)
    TE = T @ E                                   # (m, 2)
    Tk = np.einsum("kce,c->ke", T, k)            # (m, n)
    g = np.zeros(m + 2 + n)
    H = np.zeros((m + 2 + n, m + 2 + n))
    g[ia] = TE @ k
    g[ik] = a @ TE
    g[ie] = a @ Tk
    H[ia, ik] = TE
    H[ia, ie] = Tk
    H[ik, ie] = np.einsum("k,kce->ce", a, T)
    # quartic part
    Qkk = np.einsum("klcd,c,d->kl", Q, k, k)     # (m, m)
    Qaa = np.einsum("klcd,k,l->cd", Q, a, a)     # (2, 2)
    Qak = np.einsum("klcd,l,d->kc", Q, a, k) + np.einsum("klcd,l,c->kd", Q, a, k)
    g[ia] += Qkk @ a
    g[ik] += Qaa @ k
    H[ia, ia] += Qkk
    H[ik, ik] += Qaa
    H[ia, ik] += Qak
    H = np.triu(H) + np.triu(H, 1).T
    return g, H


class _CorotationalThirdOrder(CorotationalMixin):
    """Third-order beams: nodal shear strains in the node's own frame give the local slopes."""

    def _cr_setup(self) -> None:
        if getattr(self, "_cr_ready", False):
            return
        super()._cr_setup()
        # shear strains are held in the node frame N (kinked or twisted meshes); the
        # element's section axes at the node are its own frame E0 there
        Q = []
        for N in (getattr(self, "_frame_N0", None), getattr(self, "_frame_N1", None)):
            Q.append(np.eye(2) if N is None else (self._E0.T @ np.asarray(N, dtype=np.float64))[1:, 1:])
        self._slope_Q = Q

    def _local_extras(self, U: np.ndarray, p: np.ndarray) -> np.ndarray:
        Q0, Q1 = self._slope_Q
        g0 = U[:, 12:14] @ Q0.T
        g1 = U[:, 14:16] @ Q1.T
        out = U[:, 12:].copy()
        out[:, 0] = g0[:, 0] + p[:, 5]
        out[:, 1] = g0[:, 1] - p[:, 4]
        out[:, 2] = g1[:, 0] + p[:, 11]
        out[:, 3] = g1[:, 1] - p[:, 10]
        return out


class CorotationalReddyBeamElement3D(_CorotationalThirdOrder, LinearReddyBeamElement3D):
    """Corotational Reddy beam: third-order shear (shear lag of the slaved warping) at large rotation."""

    element_type_name = "Corotational-Reddy-3D"


class CorotationalLevinsonBeamElement3D(_CorotationalThirdOrder, LinearLevinsonBeamElement3D):
    """Corotational Levinson beam (non-variational third-order shear; unsymmetric tangent)."""

    element_type_name = "Corotational-Levinson-3D"

# pre_processing/element_library/linear/beam/third_order_shear_deformation_theory/third_order_beam.py
"""
2-node 3D third-order shear beam: the shared kinematics of the Reddy and Levinson elements.

**Kinematics.** The section moves as a Timoshenko section (plane rotation plus the
Saint-Venant warping of its six strains). The shear part of the warping,
``Psi(s) gamma``, is slaved to the local shear strain, so where the shear strain
changes along the beam faster than in the Saint-Venant solution it stretches the
section::

    eps_xx += Psi_y(s) chi_y + Psi_z(s) chi_z,     chi = gamma' - R psi6

``Psi`` and ``R`` come from ``SectionAnalysis.third_order_stiffness`` and
``third_order_rate``. For a solid rectangle ``Psi`` is Reddy's cubic
``z - 4 z^3 / (3 h^2)`` and the section terms are the familiar ``z^4`` and ``z^6``
integrals (``EI / 84``); for thin walls it is the shear warping from the wall
analysis. With ``chi`` free (no warping restraint) and no distributed load the
element gives the Saint-Venant (Timoshenko) beam.

**DOFs.** ``chi`` needs ``gamma'``, hence ``v''`` and ``w''``, so ``v`` and ``w`` are
cubic Hermite and each node carries their slopes as two extra DOFs:
``[u, v, w, theta_x, theta_y, theta_z, v', w']`` (``Q1 = v'``, ``Q2 = w'``). ``u``,
``theta_x``, ``theta_y`` and ``theta_z`` are linear. Element order: the 12 standard
DOFs (node 1, node 2), then ``[v'_1, w'_1, v'_2, w'_2]``. A clamp that also fixes
``Q1`` and ``Q2`` makes ``gamma = 0`` there, which restrains the shear warping
(Reddy's ``w = w' = phi = 0``).

**Strains** (element signs, ``S = diag(1, -1, 1, 1, 1, 1, 1, 1)``):
``e8 = [eps_x, kappa_y, kappa_z, gamma_xy, gamma_xz, phi_x, chi_y, chi_z]`` with
``gamma_xy = v' - theta_z``, ``gamma_xz = w' + theta_y`` (right-handed ``theta_y``) and
``D = SectionAnalysis.third_order_stiffness`` (8, 8).

**Weak forms.** Reddy (variationally consistent): ``K_e = int B8^T D B8 dx``.
Levinson (vectorial equilibrium of the classical resultants):
``K_e = int B6^T D[:6, :] B8 dx``; the higher-order resultant has no equation of
its own. Both use a 4-point Gauss rule (exact); the Kirchhoff limit
``gamma = 0`` leaves continuous linear rotations, so there is no shear locking.

**Loads.** ``F_e = int N_q^T q dx`` and ``N_q(x_p)^T P`` with ``N_q`` the field
interpolation (Hermite for ``v`` and ``w``).

**Results cache.** The Gauss cache stores ``B6`` and the 6 x 6 section stiffness, so
the standard section-force output is the Timoshenko part of the resultants.
"""

from __future__ import annotations

import logging

import numpy as np

from pre_processing.element_library.linear.beam.first_order_shear_deformation_theory.benscoter.linear_gbt_3D import (
    _hermite_x,
)
from pre_processing.element_library.linear.beam.first_order_shear_deformation_theory.timoshenko.linear_timoshenko_3D import (
    LinearTimoshenkoBeamElement3D,
)
from pre_processing.element_library.linear.beam.first_order_shear_deformation_theory.timoshenko.utilities.interpolate_loads import (
    LoadInterpolationOperator,
)

logger = logging.getLogger(__name__)

N_DOF = 16
SLOPE_DOFS = 2
_V_COLS = [1, 12, 7, 14]   # v1, v'1, v2, v'2
_W_COLS = [2, 13, 8, 15]   # w1, w'1, w2, w'2
_SIGN8 = np.array([1.0, -1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0])
_N_GAUSS = 4


class ThirdOrderShearBeamElement3D(LinearTimoshenkoBeamElement3D):
    """Shared Reddy / Levinson element (see module docstring). Subclasses set ``variational``."""

    element_type_name = "ThirdOrder-3D"
    variational = True

    @classmethod
    def _nodal_dof_count(cls, element_dictionary: dict, idx: int, warp_mesh: bool) -> int:
        return 6 + SLOPE_DOFS

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        ed = kwargs["element_dictionary"]
        idx = int(np.where(np.asarray(ed["ids"]) == self.element_id)[0][0])
        if self._warp_policy.warping_stiffness_on:
            raise NotImplementedError(
                f"Element {self.element_id}: {type(self).__name__} does not take [warping] = 1; "
                "use LinearBenscoterBeamElement3D for restrained torsional warping."
            )
        sa = self.section_analysis
        if sa.third_order_stiffness is None:
            raise ValueError(f"Element {self.element_id}: the section has no third-order shear terms.")
        self._D8 = np.asarray(sa.third_order_stiffness, dtype=np.float64)
        self._R = np.asarray(sa.third_order_rate, dtype=np.float64)
        self._idx = idx

    # ------------------------------------------------------------------ #
    # Operators
    # ------------------------------------------------------------------ #
    def D_matrix(self) -> np.ndarray:
        return self._D8

    def _D6(self) -> np.ndarray:
        return self._D8[:6, :6]

    def _B_at(self, xi: float) -> np.ndarray:
        """Element-signed ``B8`` (8, 16) at ``xi``."""
        L = self.L
        L1, L2 = 0.5 * (1.0 - xi), 0.5 * (1.0 + xi)
        dL = np.array([-1.0, 1.0]) / L
        _, dH, d2H = _hermite_x(float(xi), L)
        B = np.zeros((8, N_DOF))
        B[0, [0, 6]] = dL                       # eps = u'
        B[1, [4, 10]] = dL                      # kappa_y (right-handed) = theta_y'
        B[2, [5, 11]] = dL                      # kappa_z = theta_z'
        B[3, _V_COLS] = dH                      # gamma_xy = v' - theta_z
        B[3, [5, 11]] = [-L1, -L2]
        B[4, _W_COLS] = dH                      # gamma_xz = w' + theta_y
        B[4, [4, 10]] = [L1, L2]
        B[5, [3, 9]] = dL                       # phi = theta_x'
        B[6, _V_COLS] = d2H                     # gamma_xy'
        B[6, [5, 11]] = -dL
        B[7, _W_COLS] = d2H                     # gamma_xz'
        B[7, [4, 10]] = dL
        B *= _SIGN8[:, None]                    # element signs (kappa_y)
        B[6:] -= self._R @ B[:6]                # chi = gamma' - R psi6
        return B

    def _field_interpolation(self, xi: float, element_signs: bool = True) -> np.ndarray:
        """``N_q`` (6, 16) for ``[u, v, w, theta_x, theta_y, theta_z]`` at ``xi``.

        With ``element_signs`` the ``theta_y`` row is negated, as the Timoshenko
        element's load interpolation (``M_y`` loads in the element sign); without it the
        fields are right-handed (mass).
        """
        L1, L2 = 0.5 * (1.0 - xi), 0.5 * (1.0 + xi)
        H, _, _ = _hermite_x(float(xi), self.L)
        Nq = np.zeros((6, N_DOF))
        for r in (0, 3, 4, 5):
            Nq[r, [r, r + 6]] = [L1, L2]
        Nq[1, _V_COLS] = H
        Nq[2, _W_COLS] = H
        if element_signs:
            Nq[4] *= -1.0
        return Nq

    # ------------------------------------------------------------------ #
    # Stiffness
    # ------------------------------------------------------------------ #
    def element_stiffness_matrix(self):
        from pre_processing.element_library.gauss_point_data import ElementObject, StiffnessGaussPointData

        self._assert_logging_ready()
        D = self._D8
        detJ = self.jacobian_determinant
        Ke = np.zeros((N_DOF, N_DOF))
        gauss_cache = []
        for xi_g, w_g in zip(*np.polynomial.legendre.leggauss(_N_GAUSS)):
            B = self._B_at(xi_g)
            T = B if self.variational else B[:6]
            Dt = D if self.variational else D[:6, :]
            Ke += T.T @ Dt @ B * w_g * detJ
            Nq = self._field_interpolation(xi_g)
            gauss_cache.append(StiffnessGaussPointData(
                xi=float(xi_g), weight=float(w_g), B_matrix=B[:6].copy(), D_matrix=D[:6, :6].copy(),
                jacobian=float(detJ), shape_functions=Nq.T.copy(), shape_derivatives=B[:6].T.copy(),
            ))
        if self.variational:
            Ke = 0.5 * (Ke + Ke.T)
        if self.logger_operator:
            self.logger_operator.log_text("stiffness", f"\n=== Element {self.element_id} {self.element_type_name} K_e {Ke.shape} ===")
            self.logger_operator.log_matrix("stiffness", Ke, {"name": "Element Stiffness Matrix"})
            self.logger_operator.flush("stiffness")
        return ElementObject(
            element_id=self.element_id,
            element_type=self.element_type_name,
            K_e=Ke,
            gauss_data=gauss_cache,
            integration_scheme="Gauss-Legendre",
        )

    def linear_geometric_stiffness_matrix(self, U_e: np.ndarray) -> np.ndarray:
        raise NotImplementedError(f"{type(self).__name__} is linear only (second order arrives in a later PR).")

    # ------------------------------------------------------------------ #
    # Loads
    # ------------------------------------------------------------------ #
    def _precurvature_equivalent_load(self) -> np.ndarray:
        if not np.any(self._E_0_voigt):
            return np.zeros(N_DOF)
        E0 = np.concatenate([self._E_0_voigt, -self._R @ self._E_0_voigt])
        f_e = np.zeros(N_DOF)
        for xi_g, w_g in zip(*np.polynomial.legendre.leggauss(_N_GAUSS)):
            B = self._B_at(xi_g)
            T = B if self.variational else B[:6]
            Dt = self._D8 if self.variational else self._D8[:6, :]
            f_e += T.T @ Dt @ E0 * w_g * self.jacobian_determinant
        return f_e

    def _wall_load_vector(self) -> np.ndarray:
        Fe = np.zeros(N_DOF)
        f = self._wall_load_resultant()
        if f is None:
            return Fe
        f6 = np.asarray(f, dtype=np.float64)[:6]  # right-handed, as the DOFs
        for xi_g, w_g in zip(*np.polynomial.legendre.leggauss(_N_GAUSS)):
            Fe += self._field_interpolation(xi_g, element_signs=False).T @ f6 * w_g * self.jacobian_determinant
        return Fe

    def element_force_vector(self):
        from pre_processing.element_library.gauss_point_data import ForceGaussPointData, ForceObject
        from pre_processing.element_library.point_load_utils import add_phased_increment, point_load_phase_rad

        self._assert_logging_ready()
        Fe = self._precurvature_equivalent_load() + self._wall_load_vector()
        gauss_cache = []
        detJ = self.jacobian_determinant
        if self.distributed_load_array.size > 0:
            xi_g, w_g = np.polynomial.legendre.leggauss(max(self.quadrature_order, _N_GAUSS))
            x_g = (xi_g + 1) * (self.L / 2) + self.x_start
            q_g = LoadInterpolationOperator(
                distributed_loads_array=self.distributed_load_array,
                boundary_mode="error",
                interpolation_order="cubic",
                n_gauss_points=len(xi_g),
            ).interpolate(x_g)
            for k, (xi, w) in enumerate(zip(xi_g, w_g)):
                Nq = self._field_interpolation(xi)
                Fe += Nq.T @ q_g[k] * w * detJ
                gauss_cache.append(ForceGaussPointData(
                    xi=float(xi), weight=float(w), shape_functions=Nq.T.copy(), jacobian=float(detJ),
                    distributed_load=q_g[k].copy(),
                ))
        point_cache = None
        if self.point_load_array.size > 0:
            for load in self.point_load_array:
                x_p = float(load[0])
                if not self._is_point_in_element(x_p):
                    continue
                xi_p = 2 * (x_p - self.x_start) / self.L - 1
                inc = self._field_interpolation(xi_p).T @ load[3:9].astype(np.float64)
                Fe = add_phased_increment(Fe, inc, point_load_phase_rad(load))
            point_cache = self.point_load_array.copy()
        if self.logger_operator:
            self.logger_operator.log_matrix("force", Fe.reshape(1, -1), {"name": "Final Force Vector"})
            self.logger_operator.flush("force")
        return ForceObject(
            element_id=self.element_id,
            element_type=self.element_type_name,
            F_e=Fe,
            gauss_data=gauss_cache,
            point_loads=point_cache,
        )

    # ------------------------------------------------------------------ #
    # Mass
    # ------------------------------------------------------------------ #
    def element_mass_matrix(self):
        from pre_processing.element_library.gauss_point_data import MassObject

        self._assert_logging_ready()
        Ms = self.section_analysis.mass_matrix()
        M_e = np.zeros((N_DOF, N_DOF))
        for xi_g, w_g in zip(*np.polynomial.legendre.leggauss(_N_GAUSS)):
            Nq = self._field_interpolation(xi_g, element_signs=False)
            M_e += Nq.T @ Ms @ Nq * w_g * self.jacobian_determinant
        return MassObject(element_id=self.element_id, element_type=self.element_type_name, M_e=0.5 * (M_e + M_e.T))

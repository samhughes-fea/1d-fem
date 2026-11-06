# Unshearable geometrically exact beam (Kirchhoff): stretch + curvature on the live chord.

"""12-DOF Euler–Bernoulli core; deformational ``p`` from the Battini chord map.

Local energy is quadratic ``½ pᵀ K_e p − pᵀ f_{E0}`` (no Green–Lagrange on the
fixed mesh, no CR cubic/quartic wall energy). Global ``F = Bᵀ f`` with ``B = dp/dU``.
"""

from __future__ import annotations

import numpy as np

from onedfem.pre_processing.element_library.linear.beam.zosdt.euler_bernoulli.linear_euler_bernoulli_3D import (
    LinearEulerBernoulliBeamElement3D,
)
from onedfem.pre_processing.element_library.nonlinear.large_rotations.corotational.utilities.frame import (
    reference_frame,
)
from onedfem.pre_processing.element_library.nonlinear.large_rotations.gebt_unshearable.utilities.kinematics import (
    hessian_contraction,
    jacobian_p,
    local_force,
    local_p,
)


class GEBTUnshearableBeamElement3D(LinearEulerBernoulliBeamElement3D):
    """Unshearable GEBT: Kirchhoff stretch and curvature, linear EB section law."""

    element_type_name = "GEBTUnshearable-3D"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        X = np.asarray(self.node_coords, dtype=np.float64)
        self._X1, self._X2 = X[0].copy(), X[1].copy()
        self._E0 = reference_frame(self._X2 - self._X1)
        self._l0 = float(np.linalg.norm(self._X2 - self._X1))
        self._K_core = None
        self._f_E0 = None

    def _core(self):
        if self._K_core is None:
            logger, self.logger_operator = self.logger_operator, None
            try:
                self._K_core = np.asarray(self.element_stiffness_matrix().K_e, dtype=np.float64)
                self._f_E0 = np.asarray(self._precurvature_equivalent_load(), dtype=np.float64)
            finally:
                self.logger_operator = logger
        return self._K_core, self._f_E0

    def _p(self, U_e: np.ndarray) -> np.ndarray:
        return local_p(U_e, self._X1, self._X2, self._E0, self._l0, self.element_id)

    def internal_force_vector(self, U_e: np.ndarray) -> np.ndarray:
        U = np.asarray(U_e, dtype=np.float64).ravel()[:12]
        K, f0 = self._core()
        p = self._p(U)
        f, _ = local_force(p, K, f0)
        B = jacobian_p(U, self._X1, self._X2, self._E0, self._l0, self.element_id)
        return B.T @ f

    def tangent_stiffness_matrix(self, U_e: np.ndarray) -> np.ndarray:
        U = np.asarray(U_e, dtype=np.float64).ravel()[:12]
        K, f0 = self._core()
        p = self._p(U)
        f, H = local_force(p, K, f0)
        B = jacobian_p(U, self._X1, self._X2, self._E0, self._l0, self.element_id)
        Kt = B.T @ H @ B + hessian_contraction(U, f, self._X1, self._X2, self._E0, self._l0, self.element_id)
        return 0.5 * (Kt + Kt.T)

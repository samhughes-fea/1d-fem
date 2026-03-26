# pre_processing/element_library/nonlinear/timoshenko/utilities/geometric_stiffness.py
"""
Geometric stiffness K_sigma for 2-node 3D Timoshenko beam (Total Lagrangian).
Axial N and bending moments M_y, M_z: classical 4x4 beam-column template per bending plane
(Cook / Bathe style P/(30L) matrix) plus axial u_x coupling from N.
Moment terms use M/L^2 scaling on the same template (dimensions consistent with N/L for translations).
"""
from dataclasses import dataclass
import numpy as np


def _plane_k4(L: float) -> np.ndarray:
    return np.array(
        [
            [36.0, 3.0 * L, -36.0, 3.0 * L],
            [3.0 * L, 4.0 * L * L, -3.0 * L, -L * L],
            [-36.0, -3.0 * L, 36.0, -3.0 * L],
            [3.0 * L, -L * L, -3.0 * L, 4.0 * L * L],
        ],
        dtype=np.float64,
    )


def _embed(K: np.ndarray, idx, K4: np.ndarray) -> None:
    for a, ia in enumerate(idx):
        for b, ib in enumerate(idx):
            K[ia, ib] += K4[a, b]


@dataclass(frozen=True)
class GeometricStiffnessOperator:
    element_length: float

    def __post_init__(self) -> None:
        if self.element_length <= 0:
            raise ValueError(f"element_length must be positive, got {self.element_length}")

    def assemble_K_sigma(
        self,
        N: float,
        M_y: float,
        M_z: float,
        xi: np.ndarray,
        weights: np.ndarray,
        dN_dx: np.ndarray,
        jacobian: float,
    ) -> np.ndarray:
        L = self.element_length
        K_sigma = np.zeros((12, 12), dtype=np.float64)
        # Axial u_x (0, 6) — first-derivative coupling from N
        for k, wk in enumerate(weights):
            dN = dN_dx[k]
            for i in (0, 6):
                for j in (0, 6):
                    K_sigma[i, j] += N * dN[i, 0] * dN[j, 0] * wk * jacobian
        # Classical 4x4 templates (constant resultant along element in TL iterate)
        K4 = _plane_k4(L)
        cN = N / (30.0 * L) if L > 0 else 0.0
        cMy = M_y / (30.0 * L * L) if L > 0 else 0.0
        cMz = M_z / (30.0 * L * L) if L > 0 else 0.0
        # Bending about z (x-y plane): u_y, theta_z -> DOF 1,5,7,11
        Kz = (cN + cMz) * K4
        _embed(K_sigma, [1, 5, 7, 11], Kz)
        # Bending about y (x-z plane): u_z, theta_y -> 2,4,8,10
        Ky = (cN + cMy) * K4
        _embed(K_sigma, [2, 4, 8, 10], Ky)
        K_sigma = 0.5 * (K_sigma + K_sigma.T)
        return K_sigma

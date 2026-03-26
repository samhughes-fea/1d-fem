# pre_processing/element_library/nonlinear/euler_bernoulli/utilities/geometric_stiffness.py
"""
Geometric stiffness **K_σ** (12, 12) for 2-node 3D TL beam — **weak-form Gauss sum** only.

At each Gauss point ``g`` with weight ``w_g`` and ``|J|``:

- **Axial (u_x DOFs 0, 6):** ``K_σ[i,j] += N_g * (∂N_i/∂x)(∂N_j/∂x) w_g |J|`` for ``i,j ∈ {0,6}`` (component ``u_x``).

- **Bending XY (DOFs 1, 5, 7, 11):** let ``h_xy = [∂N_1/∂x, ∂N_5/∂x, ∂N_7/∂x, ∂N_11/∂x]`` (columns ``u_y`` / ``θ_z``). Then
  ``K_σ += w_g |J| (N_g + M_z,g / L) h_xy h_xyᵀ`` embedded on those DOFs.

- **Bending XZ (DOFs 2, 4, 8, 10):** ``h_xz = [∂N_2/∂x, -∂N_4/∂x, ∂N_8/∂x, -∂N_10/∂x]`` (``θ_y`` sign convention matches linear EB/Timoshenko shape layout). Then
  ``K_σ += w_g |J| (N_g + M_y,g / L) h_xz h_xzᵀ``.

For Euler-Bernoulli Hermite shapes, this sum matches the classical ``P/(30L)`` beam-column block for constant ``N`` and ``M``.
For Timoshenko linear shapes, the same formula uses ``dN/dx`` at each point (nonlinear Timoshenko element).
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np


def _embed(K: np.ndarray, idx: list[int], K4: np.ndarray) -> None:
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
        N_gp: np.ndarray,
        M_y_gp: np.ndarray,
        M_z_gp: np.ndarray,
        weights: np.ndarray,
        dN_dx: np.ndarray,
        jacobian: float,
    ) -> np.ndarray:
        """
        Sum ``K_sigma += ... * w_k * jacobian`` over Gauss points (see module docstring for term layout).

        Parameters
        ----------
        N_gp, M_y_gp, M_z_gp
            Section resultants at each Gauss point, shape (n_gp,).
        weights
            Gauss weights, shape (n_gp,).
        dN_dx
            dN/dx per point, shape (n_gp, 12, 6).
        jacobian
            detJ = dx/dxi (constant L/2 for straight 2-node chord map).
        """
        L = self.element_length
        K_sigma = np.zeros((12, 12), dtype=np.float64)
        for k, wk in enumerate(weights):
            dN = dN_dx[k]
            Ng = float(N_gp[k])
            Myg = float(M_y_gp[k])
            Mzg = float(M_z_gp[k])
            fac = wk * jacobian
            for i in (0, 6):
                for j in (0, 6):
                    K_sigma[i, j] += Ng * dN[i, 0] * dN[j, 0] * fac
            h_xy = np.array(
                [dN[1, 1], dN[5, 5], dN[7, 1], dN[11, 5]],
                dtype=np.float64,
            )
            h_xz = np.array(
                [dN[2, 2], -dN[4, 4], dN[8, 2], -dN[10, 4]],
                dtype=np.float64,
            )
            s_xy = Ng + Mzg / L
            s_xz = Ng + Myg / L
            _embed(K_sigma, [1, 5, 7, 11], s_xy * fac * np.outer(h_xy, h_xy))
            _embed(K_sigma, [2, 4, 8, 10], s_xz * fac * np.outer(h_xz, h_xz))
        K_sigma = 0.5 * (K_sigma + K_sigma.T)
        return K_sigma

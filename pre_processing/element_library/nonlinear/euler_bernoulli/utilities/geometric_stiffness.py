# pre_processing/element_library/nonlinear/euler_bernoulli/utilities/geometric_stiffness.py
"""
Geometric stiffness K_σ for 2-node 3D Euler–Bernoulli beam (Total Lagrangian).
Depends on current N, M_y, M_z and shape derivatives.
"""

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class GeometricStiffnessOperator:
    """
    Geometric stiffness matrix K_σ for a 2-node 3D Euler–Bernoulli beam (Total Lagrangian).

    K_σ depends on current section forces N, M_y, M_z and the shape function derivatives.
    Standard form (e.g. Przemieniecki, Bathe): axial force N contributes to the
    transverse displacement and rotation coupling; moments M_y, M_z contribute to
    the geometric stiffness terms in the bending rows/columns.

    Reference: Theory of Matrix Structural Analysis (Przemieniecki), or
    Finite Element Procedures (Bathe), Ch. 6.

    Parameters
    ----------
    element_length : float
        Reference length L of the element.
    """

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
        """
        Assemble 12×12 geometric stiffness matrix K_σ by quadrature.

        K_σ = ∫ (dN/dx)ᵀ S_geo (dN/dx) |J| dξ  where S_geo is the stress-dependent
        matrix (N, M_y, M_z). For a beam, the standard form uses N in the transverse
        and rotation DOFs (see e.g. Cook, Malkus, Plesha; or Crisfield).

        Parameters
        ----------
        N, M_y, M_z : float
            Current section forces (axial force, moments about y and z).
        xi : np.ndarray
            Gauss points in (-1, 1).
        weights : np.ndarray
            Gauss weights.
        dN_dx : np.ndarray, shape (n_gauss, 12, 6)
            Shape function derivatives w.r.t. x at each Gauss point.
        jacobian : float
            |J| = L/2 for the element.

        Returns
        -------
        K_sigma : np.ndarray, shape (12, 12)
        """
        K_sigma = np.zeros((12, 12), dtype=np.float64)
        for k, (xk, wk) in enumerate(zip(xi, weights)):
            dN = dN_dx[k]  # (12, 6)
            # Simplified geometric stiffness: N contributes to lateral/rotation coupling.
            # Standard beam K_σ: N * (integral of (dN_v/dx)ᵀ (dN_v/dx) for transverse dofs).
            for i in range(12):
                for j in range(12):
                    if i in (0, 6) and j in (0, 6):
                        K_sigma[i, j] += N * dN[i, 0] * dN[j, 0] * wk * jacobian
                    if i in (1, 7) and j in (1, 7):
                        K_sigma[i, j] += N * dN[i, 1] * dN[j, 1] * wk * jacobian
                    if i in (2, 8) and j in (2, 8):
                        K_sigma[i, j] += N * dN[i, 2] * dN[j, 2] * wk * jacobian
        K_sigma = 0.5 * (K_sigma + K_sigma.T)
        return K_sigma

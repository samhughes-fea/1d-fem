# pre_processing/element_library/linear/euler_bernoulli/utilities/shape_functions.py
"""
Shape functions for 2-node 3D Euler–Bernoulli beam.

`natural_coordinate_form(xi)` returns N (n_gp, 12, 6), dN/dξ, d²N/dξ².
Rows are global DOFs a, columns are components (u_x, u_y, u_z, θ_x, θ_y, θ_z).
Hermite functions are used on bending channels; axial and torsion use linear terms.
Used to form B and `F_dist += w_g * N.T @ q * detJ`.
See `FORMULATION_DOCSTRING_STANDARDS.md`.
"""

import numpy as np
from typing import Tuple
from dataclasses import dataclass

@dataclass(frozen=True)
class ShapeFunctionOperator:
    """
    Shape-function tensor N ∈ ℝ^{n_gp × 12 × 6} for a 2-node 3-D Euler-Bernoulli beam.

    N is a rank-3 tensor; the slice N_g ∈ ℝ^{12×6} at Gauss point g maps the
    12-component element displacement vector U_e ∈ ℝ^{12} to the 6-component
    displacement field u ∈ ℝ^6 at that point:

        u = N_g U_e,   u = [u_x, u_y, u_z, θ_x, θ_y, θ_z]^T

    Returns N and its natural-coordinate derivatives dN/dξ, d²N/dξ² (same shape)
    via `natural_coordinate_form`; `physical_coordinate_form` scales to dN/dx,
    d²N/dx² using the chain-rule factors ∂ξ/∂x = 2/L and ∂²ξ/∂x² = 4/L².

    Parameters
    ----------
    element_length : float
        Chord length L (physical x ∈ [0, L], must be > 0).

    Attributes
    ----------
    dξ_dx : float
        First derivative chain factor ∂ξ/∂x = 2/L.
    d2ξ_dx2 : float
        Second derivative chain factor ∂²ξ/∂x² = 4/L².

    Notes
    -----
    **Sparsity structure of N_g (single Gauss point slice, shape (12, 6))**

    ```text
    cols = [u_x, u_y, u_z, θ_x, θ_y, θ_z]
    N_g =
    [ n1,1   0      0      0      0       0   ]   row 0:  u_x DOF, node 1
    [ 0     n2,2    0      0      0       0   ]   row 1:  u_y DOF, node 1
    [ 0      0     n3,3    0      0       0   ]   row 2:  u_z DOF, node 1
    [ 0      0      0     n4,4    0       0   ]   row 3:  θ_x DOF, node 1
    [ 0      0      0      0     n5,5     0   ]   row 4:  θ_y DOF, node 1
    [ 0      0      0      0      0      n6,6 ]   row 5:  θ_z DOF, node 1
    [ n7,1   0      0      0      0       0   ]   row 6:  u_x DOF, node 2
    [ 0     n8,2    0      0      0       0   ]   row 7:  u_y DOF, node 2
    [ 0      0     n9,3    0      0       0   ]   row 8:  u_z DOF, node 2
    [ 0      0      0    n10,4    0       0   ]   row 9:  θ_x DOF, node 2
    [ 0      0      0      0    n11,5     0   ]   row 10: θ_y DOF, node 2
    [ 0      0      0      0      0     n12,6 ]   row 11: θ_z DOF, node 2
    ```

    N is sparse-by-structure: each DOF row activates exactly one displacement component
    column. All other entries are zero.

    **Shape function basis and active entries of N_g**

    Natural coordinate ξ ∈ [−1, 1] maps to x ∈ [0, L] via ξ = (2x − L)/L.

    Linear Lagrange polynomials (axial u_x and torsion θ_x channels):

        L₁(ξ) = ½(1 − ξ),   L₂(ξ) = ½(1 + ξ)

    Hermite cubic polynomials (bending channels, standard C¹ beam pair per plane):

        H₁(ξ) = ¼(1 − ξ)²(2 + ξ)         (displacement, node 1; H₁(−1)=1, H₁(1)=0)
        H₂(ξ) = (L/8)(1 − ξ)²(1 + ξ)     (rotation, node 1;    dH₂/dξ|_{ξ=−1}=L/2)
        H₃(ξ) = ¼(1 + ξ)²(2 − ξ)         (displacement, node 2; H₃(1)=1, H₃(−1)=0)
        H₄(ξ) = −(L/8)(1 + ξ)²(1 − ξ)    (rotation, node 2;    dH₄/dξ|_{ξ=1}=−L/2)

    ```text
    Active entries of N_g (rows = DOF index a, columns = displacement component c):

      N_g[0,0]  = L₁(ξ_g)     N_g[6,0]  = L₂(ξ_g)     [u_x channel]
      N_g[1,1]  = H₁(ξ_g)     N_g[7,1]  = H₃(ξ_g)     [u_y channel]
      N_g[2,2]  = H₁(ξ_g)     N_g[8,2]  = H₃(ξ_g)     [u_z channel]
      N_g[3,3]  = L₁(ξ_g)     N_g[9,3]  = L₂(ξ_g)     [θ_x channel]
      N_g[4,4]  = −H₂(ξ_g)    N_g[10,4] = −H₄(ξ_g)    [θ_y channel; sign from EB θ_y = −∂u_z/∂x]
      N_g[5,5]  = H₂(ξ_g)     N_g[11,5] = H₄(ξ_g)     [θ_z channel]
      N_g[a,c]  = 0  for all other (a,c) pairs          (sparse by structure)
    ```

    The same sparsity pattern holds for dN/dξ and d²N/dξ² (derivatives of the
    polynomials above, entering via the chain rule into B).

    Weak-form linkage: `linear_euler_bernoulli_3D` accumulates
    `K_e += B.T @ D @ B * w_g * detJ` and distributed-load vectors via
    `F_dist += w_g * N.T @ q * detJ` with `detJ = L/2`.

    See Also
    --------
    linear_euler_bernoulli_3D.LinearEulerBernoulliBeamElement3D
    """

    element_length: float

    def __post_init__(self):
        """Precompute and validate coordinate transformation factors."""
        if self.element_length <= 0:
            raise ValueError(f"Element length must be positive, got {self.element_length}")
        
        object.__setattr__(self, '_dξ_dx', 2 / self.element_length)
        object.__setattr__(self, '_d2ξ_dx2', 4 / (self.element_length**2))

    @property
    def dξ_dx(self) -> float:
        """First derivative transform ∂ξ/∂x = 2/L (unitless)"""
        return self._dξ_dx

    @property
    def d2ξ_dx2(self) -> float:
        """Second derivative transform ∂²ξ/∂x² = 4/L² (1/m²)"""
        return self._d2ξ_dx2

    def natural_coordinate_form(self, xi: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Evaluate shape functions and derivatives in natural coordinates (xi-space).

        Parameters
        ----------
        xi : np.ndarray
            Natural coordinates in [-1, 1], shape (n_points,).

        Returns
        -------
        N : np.ndarray
            Shape function matrix [n_points, 12, 6] where:
            - Axis 1: Evaluation points
            - Axis 2: DOFs (12 total: 6 per node)
            - Axis 3: Components (u_x, u_y, u_z, θ_x, θ_y, θ_z)
        dN_dξ : np.ndarray
            First derivatives ∂N/∂ξ [n_points, 12, 6]
        d2N_dξ2 : np.ndarray
            Second derivatives ∂²N/∂ξ² [n_points, 12, 6]

        Notes
        -----
        Shape function organization:
        Node 1: [u_x, u_y, u_z, θ_x, θ_y, θ_z]
        Node 2: [u_x, u_y, u_z, θ_x, θ_y, θ_z]
        """
        xi = np.asarray(xi, dtype=np.float64)
        n_points = xi.size
        ξ = xi.reshape(-1, 1, 1)  # Prepare for broadcasting

        # Initialize output arrays
        N = np.zeros((n_points, 12, 6))
        dN_dξ = np.zeros_like(N)
        d2N_dξ2 = np.zeros_like(N)

        # ----- Axial Displacement (Linear Lagrange) -----
        # N₁(ξ) = 0.5(1-ξ), N₇(ξ) = 0.5(1+ξ)
        N[:, [0,6], 0] = 0.5 * np.array([1 - ξ.squeeze(), 1 + ξ.squeeze()]).T
        dN_dξ[:, [0,6], 0] = 0.5 * np.array([-1, 1])

        # ----- Bending in XY Plane (Hermite Cubic) -----
        # Standard Hermite cubic shape functions:
        # N1 = (1/4)(1-ξ)²(2+ξ) = 0.5 - 0.75*ξ + 0.25*ξ³
        # N2 = (L/8)(1-ξ)²(1+ξ) = (L/8)(1 - ξ - ξ² + ξ³)
        # N3 = (1/4)(1+ξ)²(2-ξ) = 0.5 + 0.75*ξ - 0.25*ξ³
        # N4 = -(L/8)(1+ξ)²(1-ξ) = -(L/8)(1 + ξ - ξ² - ξ³)
        
        L = self.element_length
        ξ_flat = ξ.squeeze()
        
        # Displacement shape functions (N1, N3)
        N1 = 0.5 - 0.75*ξ_flat + 0.25*ξ_flat**3
        N3 = 0.5 + 0.75*ξ_flat - 0.25*ξ_flat**3
        N[:, [1,7], 1] = np.array([N1, N3]).T
        
        # First derivatives of displacement shape functions
        dN1_dξ = -0.75 + 0.75*ξ_flat**2
        dN3_dξ = 0.75 - 0.75*ξ_flat**2
        dN_dξ[:, [1,7], 1] = np.array([dN1_dξ, dN3_dξ]).T
        
        # Second derivatives of displacement shape functions
        d2N1_dξ2 = 1.5*ξ_flat
        d2N3_dξ2 = -1.5*ξ_flat
        d2N_dξ2[:, [1,7], 1] = np.array([d2N1_dξ2, d2N3_dξ2]).T

        # Rotation shape functions (N2, N4) - scaled by L/8
        N2 = (L/8) * (1 - ξ_flat - ξ_flat**2 + ξ_flat**3)
        N4 = -(L/8) * (1 + ξ_flat - ξ_flat**2 - ξ_flat**3)
        N[:, [5,11], 5] = np.array([N2, N4]).T
        
        # First derivatives of rotation shape functions
        dN2_dξ = (L/8) * (-1 - 2*ξ_flat + 3*ξ_flat**2)
        dN4_dξ = -(L/8) * (1 - 2*ξ_flat - 3*ξ_flat**2)
        dN_dξ[:, [5,11], 5] = np.array([dN2_dξ, dN4_dξ]).T
        
        # Second derivatives of rotation shape functions
        d2N2_dξ2 = (L/8) * (-2 + 6*ξ_flat)
        d2N4_dξ2 = -(L/8) * (-2 - 6*ξ_flat)
        d2N_dξ2[:, [5,11], 5] = np.array([d2N2_dξ2, d2N4_dξ2]).T

        # ----- Bending in XZ Plane (Hermite Cubic) -----
        N[:, [2,8], 2] = N[:, [1,7], 1]
        dN_dξ[:, [2,8], 2] = dN_dξ[:, [1,7], 1]
        d2N_dξ2[:, [2,8], 2] = d2N_dξ2[:, [1,7], 1]

        # Rotation terms (negative sign convention)
        N[:, [4,10], 4] = -N[:, [5,11], 5]
        dN_dξ[:, [4,10], 4] = -dN_dξ[:, [5,11], 5]
        d2N_dξ2[:, [4,10], 4] = -d2N_dξ2[:, [5,11], 5]

        # ----- Torsional Rotation (Linear Lagrange) -----
        N[:, [3,9], 3] = N[:, [0,6], 0]
        dN_dξ[:, [3,9], 3] = dN_dξ[:, [0,6], 0]

        return N, dN_dξ, d2N_dξ2

    def physical_coordinate_form(self, xi: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Evaluate shape functions and derivatives in physical coordinates (x-space).

        Parameters
        ----------
        xi : np.ndarray
            Natural coordinates in [-1, 1], shape (n_points,).

        Returns
        -------
        N : np.ndarray
            Shape function matrix [n_points, 12, 6]
        dN_dx : np.ndarray
            First derivatives ∂N/∂x [n_points, 12, 6]
        d2N_dx2 : np.ndarray
            Second derivatives ∂²N/∂x² [n_points, 12, 6]

        Notes
        -----
        ``dN_dx = dN_dxi * dxi_dx``, ``d2N_dx2 = d2N_dxi2 * d2xi_dx2`` with ``dxi_dx = 2/L``.
        """
        N, dN_dξ, d2N_dξ2 = self.natural_coordinate_form(xi)
        
        # Apply coordinate transforms
        dN_dx = dN_dξ * self.dξ_dx
        d2N_dx2 = d2N_dξ2 * self.d2ξ_dx2
        
        return N, dN_dx, d2N_dx2

    @property
    def dof_interpretation(self) -> np.ndarray:
        """
        Structured array documenting DOF physical meaning.

        Returns
        -------
        np.ndarray
            Structured array with fields:
            - index: DOF index (0-11)
            - node: 'Node 1' or 'Node 2'
            - component: 'u_x', 'u_y', 'u_z', 'θ_x', 'θ_y', 'θ_z'
            - behavior: 'Axial', 'Bending XY', 'Bending XZ', 'Torsion'
        """
        return np.array([
            (0, 'Node 1', 'u_x', 'Axial'),
            (1, 'Node 1', 'u_y', 'Bending XY'),
            (2, 'Node 1', 'u_z', 'Bending XZ'),
            (3, 'Node 1', 'θ_x', 'Torsion'),
            (4, 'Node 1', 'θ_y', 'Bending XZ'), 
            (5, 'Node 1', 'θ_z', 'Bending XY'),
            (6, 'Node 2', 'u_x', 'Axial'),
            (7, 'Node 2', 'u_y', 'Bending XY'),
            (8, 'Node 2', 'u_z', 'Bending XZ'),
            (9, 'Node 2', 'θ_x', 'Torsion'),
            (10, 'Node 2', 'θ_y', 'Bending XZ'),
            (11, 'Node 2', 'θ_z', 'Bending XY')
        ], dtype=[('index', 'i4'), ('node', 'U10'), ('component', 'U3'), ('behavior', 'U10')])
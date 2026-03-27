# pre_processing\element_library\levinson\utilities\shape_functions.py
"""
Shape functions for 2-node 3-D Levinson beam (higher-order transverse fields).

`natural_coordinate_form(xi)` returns N, dN/dξ, and d²N/dξ², each with shape (n_gp, 12, 6):
row = element DOF index, column = component (u_x, u_y, u_z, θ_x, θ_y, θ_z).

Weak-form force assembly uses `F_dist += w_g * N.T @ q * detJ` with `detJ = L/2` and ξ in [-1, 1].

Compared with EB/Timoshenko (same 12-DOF contract), Levinson keeps the same N layout but uses
quintic transverse and cubic bending-rotation fields, and B uses κ_z before κ_y in the strain vector.

See Also
--------
linear_levinson_3D.LinearLevinsonBeamElement3D
"""

import numpy as np
from typing import Tuple
from dataclasses import dataclass

@dataclass(frozen=True)
class ShapeFunctionOperator:
    """
    Shape-function tensor N ∈ ℝ^{n_gp × 12 × 6} for a 2-node 3-D Levinson beam.

    N is a rank-3 tensor; the slice N_g ∈ ℝ^{12×6} at Gauss point g maps the
    12-component element displacement vector U_e ∈ ℝ^{12} to the 6-component
    displacement field u ∈ ℝ^6 at that point:

        u = N_g U_e,   u = [u_x, u_y, u_z, θ_x, θ_y, θ_z]^T

    Returns N and its natural-coordinate derivatives dN/dξ, d²N/dξ² (same shape)
    via `natural_coordinate_form`; `physical_coordinate_form` scales to dN/dx,
    d²N/dx² using the chain-rule factors ∂ξ/∂x = 2/L and ∂²ξ/∂x² = 4/L².

    Compared with Euler-Bernoulli (same 12-DOF contract), Levinson uses quintic
    polynomials for transverse displacements (u_y, u_z) and cubic polynomials for
    bending rotations (θ_y, θ_z). Higher-order shear enters via the B operator
    through the α ∂²θ/∂x² terms, not through N directly.

    Parameters
    ----------
    element_length : float
        Chord length L (must be > 0).

    Attributes
    ----------
    dξ_dx : float
        First derivative chain rule factor ∂ξ/∂x = 2/L.
    d2ξ_dx2 : float
        Second derivative chain rule factor ∂²ξ/∂x² = 4/L².

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

    **Shape function polynomials and active entries of N_g**

    Natural coordinate ξ ∈ [−1, 1]. Chain rule: ∂N/∂x = (∂N/∂ξ)(2/L), ∂²N/∂x² = (∂²N/∂ξ²)(4/L²).

    Linear Lagrange (axial u_x and torsion θ_x channels):

        L₁(ξ) = ½(1 − ξ),   L₂(ξ) = ½(1 + ξ)

    Quintic polynomials (transverse displacement channels u_y, u_z);
    boundary conditions: value 1 at own node, value 0 at other node, zero slope at both nodes:

        Q₁(ξ) = ½ − (15/16)ξ + (5/8)ξ³ − (3/16)ξ⁵   (node 1; Q₁(−1)=1, Q₁(1)=0, Q₁'(±1)=0)
        Q₂(ξ) = ½ + (15/16)ξ − (5/8)ξ³ + (3/16)ξ⁵   (node 2; Q₂(1)=1, Q₂(−1)=0, Q₂'(±1)=0)

    Cubic polynomials (bending rotation channels θ_y, θ_z):

        C₁(ξ) = ¼(2 − 3ξ + ξ³)   (node 1; C₁(−1)=1, C₁(1)=0)
        C₂(ξ) = ¼(2 + 3ξ − ξ³)   (node 2; C₂(1)=1, C₂(−1)=0)

    ```text
    Active entries of N_g (rows = DOF index a, columns = displacement component c):

      N_g[0,0]  = L₁(ξ_g)     N_g[6,0]  = L₂(ξ_g)     [u_x channel]
      N_g[1,1]  = Q₁(ξ_g)     N_g[7,1]  = Q₂(ξ_g)     [u_y channel]
      N_g[2,2]  = Q₁(ξ_g)     N_g[8,2]  = Q₂(ξ_g)     [u_z channel]
      N_g[3,3]  = L₁(ξ_g)     N_g[9,3]  = L₂(ξ_g)     [θ_x channel]
      N_g[4,4]  = −C₁(ξ_g)    N_g[10,4] = −C₂(ξ_g)    [θ_y channel; sign from EB-type convention]
      N_g[5,5]  = C₁(ξ_g)     N_g[11,5] = C₂(ξ_g)     [θ_z channel]
      N_g[a,c]  = 0  for all other (a,c) pairs          (sparse by structure)
    ```

    The strain linkage in B_matrix uses:
        κ_z = ∂θ_z/∂x = (dC₁/dξ)(2/L)·θ_z¹ + (dC₂/dξ)(2/L)·θ_z²
        γ_xy = ∂u_y/∂x − θ_z + α ∂²θ_z/∂x²   (higher-order correction)

    Weak-form linkage: `K_e += B.T @ D @ B * w_g * detJ` and
    `F_dist += w_g * N.T @ q * detJ` with `detJ = L/2`.

    See Also
    --------
    linear_levinson_3D.LinearLevinsonBeamElement3D
    """

    element_length: float

    def __post_init__(self):
        if self.element_length <= 0:
            raise ValueError(f"Element length must be positive, got {self.element_length}")
        object.__setattr__(self, '_dξ_dx', 2 / self.element_length)
        object.__setattr__(self, '_d2ξ_dx2', 4 / (self.element_length**2))

    @property
    def dξ_dx(self) -> float:
        return self._dξ_dx

    @property
    def d2ξ_dx2(self) -> float:
        return self._d2ξ_dx2

    def natural_coordinate_form(self, ξ: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        ξ = np.asarray(ξ, dtype=np.float64).reshape(-1, 1, 1)
        n_points = ξ.size
        N = np.zeros((n_points, 12, 6))
        dN_dξ = np.zeros_like(N)
        d2N_dξ2 = np.zeros_like(N)

        # ----- Axial Displacement (Linear) -----
        N[:, [0,6], 0] = 0.5 * np.array([1 - ξ.squeeze(), 1 + ξ.squeeze()]).T
        dN_dξ[:, [0,6], 0] = 0.5 * np.array([-1, 1])

        # ----- Bending XY Plane (Quintic u_y, Cubic θ_z) -----
        # u_y: Quintic with N1_v(-1)=1, N1_v(1)=0, N2_v(-1)=0, N2_v(1)=1; zero slope at nodes
        # N1_v = 1/2 - (15/16)*ξ + (5/8)*ξ^3 - (3/16)*ξ^5, N2_v = 1/2 + (15/16)*ξ - (5/8)*ξ^3 + (3/16)*ξ^5
        ξ_flat = ξ.squeeze()
        N1_v = 0.5 - (15.0 / 16.0) * ξ_flat + (5.0 / 8.0) * ξ_flat**3 - (3.0 / 16.0) * ξ_flat**5
        N2_v = 0.5 + (15.0 / 16.0) * ξ_flat - (5.0 / 8.0) * ξ_flat**3 + (3.0 / 16.0) * ξ_flat**5
        N[:, [1, 7], 1] = np.array([N1_v, N2_v]).T

        # First derivatives
        dN1_v_dξ = -(15.0 / 16.0) + (15.0 / 8.0) * ξ_flat**2 - (15.0 / 16.0) * ξ_flat**4
        dN2_v_dξ = (15.0 / 16.0) - (15.0 / 8.0) * ξ_flat**2 + (15.0 / 16.0) * ξ_flat**4
        dN_dξ[:, [1, 7], 1] = np.array([dN1_v_dξ, dN2_v_dξ]).T

        # Second derivatives
        d2N1_v_dξ2 = (15.0 / 4.0) * ξ_flat - (15.0 / 4.0) * ξ_flat**3
        d2N2_v_dξ2 = -(15.0 / 4.0) * ξ_flat + (15.0 / 4.0) * ξ_flat**3
        d2N_dξ2[:, [1, 7], 1] = np.array([d2N1_v_dξ2, d2N2_v_dξ2]).T

        # θ_z: Cubic with N1_θ(-1)=1, N1_θ(1)=0, N2_θ(-1)=0, N2_θ(1)=1
        # N1_θ = (1/4)(2 - 3*ξ + ξ^3), N2_θ = (1/4)(2 + 3*ξ - ξ^3)
        N[:, [5, 11], 5] = 0.25 * np.array([
            2.0 - 3.0 * ξ_flat + ξ_flat**3,
            2.0 + 3.0 * ξ_flat - ξ_flat**3,
        ]).T
        dN_dξ[:, [5, 11], 5] = 0.25 * np.array([
            -3.0 + 3.0 * ξ_flat**2,
            3.0 - 3.0 * ξ_flat**2,
        ]).T
        d2N_dξ2[:, [5, 11], 5] = 0.25 * np.array([
            6.0 * ξ_flat,
            -6.0 * ξ_flat,
        ]).T

        # ----- Bending XZ Plane (Mirror XY with sign adjustments) -----
        N[:, [2,8], 2] = N[:, [1,7], 1]
        dN_dξ[:, [2,8], 2] = dN_dξ[:, [1,7], 1]
        d2N_dξ2[:, [2,8], 2] = d2N_dξ2[:, [1,7], 1]
        
        N[:, [4,10], 4] = -N[:, [5,11], 5]
        dN_dξ[:, [4,10], 4] = -dN_dξ[:, [5,11], 5]
        d2N_dξ2[:, [4,10], 4] = -d2N_dξ2[:, [5,11], 5]

        # ----- Torsional Rotation (Linear) -----
        N[:, [3,9], 3] = N[:, [0,6], 0]
        dN_dξ[:, [3,9], 3] = dN_dξ[:, [0,6], 0]

        return N, dN_dξ, d2N_dξ2

    def physical_coordinate_form(self, ξ: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        N, dN_dξ, d2N_dξ2 = self.natural_coordinate_form(ξ)
        dN_dx = dN_dξ * self.dξ_dx
        d2N_dx2 = d2N_dξ2 * self.d2ξ_dx2
        return N, dN_dx, d2N_dx2

    @property
    def dof_interpretation(self) -> np.ndarray:
        return np.array([
            (0, 'Node 1', 'u_x', 'Axial'),
            (1, 'Node 1', 'u_y', 'Bending XY (Quintic)'),
            (2, 'Node 1', 'u_z', 'Bending XZ (Quintic)'),
            (3, 'Node 1', 'θ_x', 'Torsion'),
            (4, 'Node 1', 'θ_y', 'Bending XZ (Cubic)'), 
            (5, 'Node 1', 'θ_z', 'Bending XY (Cubic)'),
            (6, 'Node 2', 'u_x', 'Axial'),
            (7, 'Node 2', 'u_y', 'Bending XY (Quintic)'),
            (8, 'Node 2', 'u_z', 'Bending XZ (Quintic)'),
            (9, 'Node 2', 'θ_x', 'Torsion'),
            (10, 'Node 2', 'θ_y', 'Bending XZ (Cubic)'),
            (11, 'Node 2', 'θ_z', 'Bending XY (Cubic)')
        ], dtype=[('index', 'i4'), ('node', 'U10'), ('component', 'U3'), ('behavior', 'U20')])
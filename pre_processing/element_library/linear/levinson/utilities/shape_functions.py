# pre_processing\element_library\levinson\utilities\shape_functions.py

import numpy as np
from typing import Tuple
from dataclasses import dataclass

@dataclass(frozen=True)
class ShapeFunctionOperator:
    """
    Operator for 3D Levinson beam with higher-order parabolic shear deformation.

    Mathematical Formulation
    -----------------------
    DOF Interpolation:
    - Axial displacement (u_x): 
      Linear Lagrange polynomials (O(ξ^1))
        N1_u = 0.5(1 - ξ), N2_u = 0.5(1 + ξ)
    
    - Transverse displacements (u_y, u_z):
      Quintic polynomials (O(ξ^5)) with C1 continuity at nodes:
        N1_v = 1/2 - (15/16)*ξ + (5/8)*ξ^3 - (3/16)*ξ^5
        N2_v = 1/2 + (15/16)*ξ - (5/8)*ξ^3 + (3/16)*ξ^5
        Satisfies: N1_v(-1)=1, N1_v(1)=0, N2_v(-1)=0, N2_v(1)=1; N1_v'=N2_v'=0 at ±1.
    
    - Rotations (θ_y, θ_z):
      Cubic polynomials (O(ξ^3)) with C1 continuity:
        N1_θ = (1/4)(2 - 3*ξ + ξ^3)
        N2_θ = (1/4)(2 + 3*ξ - ξ^3)
        Satisfies: N1_θ(-1)=1, N1_θ(1)=0, N2_θ(-1)=0, N2_θ(1)=1.
    
    - Torsional rotation (θ_x):
      Linear Lagrange polynomials (O(ξ^1)) same as axial displacement

    Strain Formulation:
    - Axial strain: ε_xx = du_x/dx
    - Shear strains: γ_xy = du_y/dx - θ_z + α(d²θ_z/dx²)
                    γ_xz = du_z/dx - θ_y + α(d²θ_y/dx²)
    - Curvatures: κ_z = dθ_z/dx (x-y plane), κ_y = dθ_y/dx (x-z plane), φ_x = dθ_x/dx (torsion)

    Element Characteristics:
    - 2 nodes, 6 DOFs per node (u_x, u_y, u_z, θ_x, θ_y, θ_z)
    - Higher-order interpolation eliminates shear correction factors
    - Parabolic shear stress distribution (O(ξ^2))
    - Consistent with Levinson's improved shear deformation theory

    Coordinate Transformation:
    - Requires second derivatives for curvature terms:
      ∂²/∂x² = (4/L²)∂²/∂ξ²
    - Jacobian scaling factors applied to all strain components
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
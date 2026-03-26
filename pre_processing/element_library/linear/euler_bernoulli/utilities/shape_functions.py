# pre_processing/element_library/linear/euler_bernoulli/utilities/shape_functions.py
"""
Shape functions for 2-node 3D Euler–Bernoulli beam.

``natural_coordinate_form(xi)`` returns ``N`` (n_gp, 12, 6), ``dN_dxi``, ``d2N_dxi2`` — row ``a`` = global DOF ``a``,
column ``c`` = ``(u_x,u_y,u_z,theta_x,theta_y,theta_z)``. Hermite on bending DOFs, linear on axial and torsion.
Used to form ``B`` and ``F_dist += w_g * N.T @ q * detJ``. See ``FORMULATION_DOCSTRING_STANDARDS.md``.
"""

import numpy as np
from typing import Tuple
from dataclasses import dataclass

@dataclass(frozen=True)
class ShapeFunctionOperator:
    """
    Evaluate 3D Euler-Bernoulli beam shape functions and natural derivatives on ``xi``.

    ``natural_coordinate_form`` returns ``N``, ``dN_dxi``, ``d2N_dxi2`` (batch ``n_points``, 12, 6);
    ``physical_coordinate_form`` scales to ``dN/dx``, ``d2N/dx2`` via ``dxi_dx``, ``d2xi_dx2``.

    Parameters
    ----------
    element_length : float
        Chord length ``L`` (physical ``x`` in ``[0, L]``, must be > 0).

    Attributes
    ----------
    dξ_dx : float
        First derivative chain factor ``2/L``.
    d2ξ_dx2 : float
        Second derivative chain factor ``4/L**2``.

    Notes
    -----
    Canonical `N` block (single Gauss point slice `N_g`, shape `(12,6)`):

    ```text
    cols = [u_x, u_y, u_z, θ_x, θ_y, θ_z]
    N_g =
    [ n1,1   0      0      0      0       0   ]
    [ 0     n2,2    0      0      0       0   ]
    [ 0      0     n3,3    0      0       0   ]
    [ 0      0      0     n4,4    0       0   ]
    [ 0      0      0      0     n5,5     0   ]
    [ 0      0      0      0      0      n6,6 ]
    [ n7,1   0      0      0      0       0   ]
    [ 0     n8,2    0      0      0       0   ]
    [ 0      0     n9,3    0      0       0   ]
    [ 0      0      0    n10,4    0       0   ]
    [ 0      0      0      0    n11,5     0   ]
    [ 0      0      0      0      0     n12,6 ]
    ```

    **N tensor contract:** ``N``/``dN_dxi``/``d2N_dxi2`` all use ``(n_gp, 12, 6)``.
    Rows are DOFs in node-major 12-DOF order; columns are components
    ``(u_x, u_y, u_z, theta_x, theta_y, theta_z)``.
    Unused row/column combinations remain zero (sparse-by-structure).

    **Formulation (natural coordinate xi in [-1, 1]):** axial ``u_x`` and torsion ``theta_x`` use linear
    Lagrange; transverse ``u_y``, ``u_z`` and bending rotations use Hermite cubics (standard beam pair per plane).
    Map ``xi = (2*x - L)/L`` on the chord; ``dN/dx = dN_dxi * dxi_dx``, ``d2N/dx2 = d2N_dxi2 * d2xi_dx2``.

    Weak-form linkage: ``linear_euler_bernoulli_3D``.

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
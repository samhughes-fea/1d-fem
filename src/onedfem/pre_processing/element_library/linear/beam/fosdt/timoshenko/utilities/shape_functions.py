# pre_processing/element_library/linear/timoshenko/utilities/shape_functions.py
"""
Shape functions for 2-node 3D Timoshenko beam.

``natural_coordinate_form(xi)`` returns ``N`` (n_gp, 12, 6), ``dN_dxi``, ``d2N_dxi2`` — row ``a`` = global DOF index ``a``,
column ``c`` = ``(u_x,u_y,u_z,theta_x,theta_y,theta_z)``. Used with ``F_dist += w_g * N.T @ q * detJ`` and to build ``B``.
See ``docs/conventions/FORMULATION_DOCSTRING_STANDARDS.md``.
"""

import numpy as np
from typing import Tuple
from dataclasses import dataclass

@dataclass(frozen=True)
class ShapeFunctionOperator:
    """
    Operator for evaluating 3D Timoshenko beam shape functions and their derivatives.
    Provides rigorous transformation between natural (ξ ∈ [-1,1]) and physical (x ∈ [0,L]) coordinates.

    Mathematical Formulation
    -----------------------
    Shape functions use independent displacement and rotation fields (no θ = du/dx):
    - Axial displacement: Linear Lagrange polynomials
    - Bending displacement and rotation: Linear Lagrange for both u_y, u_z and θ_z, θ_y
      so that shear strain γ = du/dx − θ is not forced to zero (avoids shear locking).
    - Torsional rotation: Linear Lagrange polynomials

    Coordinate Transformation:
    - Physical to natural: ξ = (2x - L)/L
    - Derivatives:
      ∂N/∂x = (∂N/∂ξ)(∂ξ/∂x) = (∂N/∂ξ)(2/L)
      ∂²N/∂x² = (∂²N/∂ξ²)(∂ξ/∂x)² = (∂²N/∂ξ²)(4/L²)

    Parameters
    ----------
    element_length : float
        Physical length of element (x ∈ [0,L], L > 0)

    Attributes
    ----------
    dξ_dx : float
        First derivative transform (∂ξ/∂x = 2/L)
    d2ξ_dx2 : float
        Second derivative transform (∂²ξ/∂x² = 4/L²)

    Notes
    -----
    Natural coordinate ``xi`` in [-1, 1]; physical ``x`` along chord length ``L``. Independent ``u`` and ``theta`` fields
    (no EB constraint ``theta = du/dx``) so Timoshenko shear is non-singular.

    See Also
    --------
    linear_timoshenko_3D.LinearTimoshenkoBeamElement3D
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

    def natural_coordinate_form(self, ξ: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Evaluate shape functions and derivatives in natural coordinates (ξ-space).

        Parameters
        ----------
        ξ : np.ndarray
            Natural coordinates ∈ [-1, 1] with shape (n_points,)

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
            Second derivatives ∂²N/∂ξ² [n_points, 12, 6] (zero for bending DOFs; B-matrix does not use them)

        Notes
        -----
        Shape function organization:
        Node 1: [u_x, u_y, u_z, θ_x, θ_y, θ_z]
        Node 2: [u_x, u_y, u_z, θ_x, θ_y, θ_z]
        Bending uses linear Lagrange so u and θ are independent; γ = du/dx − θ can be non-zero.
        """
        ξ = np.asarray(ξ, dtype=np.float64)
        n_points = ξ.size
        ξ = ξ.reshape(-1, 1, 1)  # Prepare for broadcasting
        ξ_flat = ξ.squeeze()

        # Initialize output arrays
        N = np.zeros((n_points, 12, 6))
        dN_dξ = np.zeros_like(N)
        d2N_dξ2 = np.zeros_like(N)

        # ----- Axial Displacement (Linear Lagrange) -----
        # N₁(ξ) = 0.5(1-ξ), N₇(ξ) = 0.5(1+ξ)
        N[:, [0, 6], 0] = 0.5 * np.array([1 - ξ_flat, 1 + ξ_flat]).T
        dN_dξ[:, [0, 6], 0] = 0.5 * np.array([-1, 1])

        # ----- Bending in XY Plane (Linear Lagrange – Timoshenko) -----
        # u_y and θ_z independent: N1 = (1-ξ)/2, N2 = (1+ξ)/2 for each
        N1 = 0.5 * (1 - ξ_flat)
        N2 = 0.5 * (1 + ξ_flat)
        N[:, [1, 7], 1] = np.array([N1, N2]).T
        dN_dξ[:, [1, 7], 1] = 0.5 * np.array([-1, 1])
        # d2N_dξ2 remains zero

        N[:, [5, 11], 5] = np.array([N1, N2]).T   # θ_z
        dN_dξ[:, [5, 11], 5] = 0.5 * np.array([-1, 1])

        # ----- Bending in XZ Plane (Linear Lagrange – Timoshenko) -----
        N[:, [2, 8], 2] = N[:, [1, 7], 1]   # u_z same as u_y
        dN_dξ[:, [2, 8], 2] = dN_dξ[:, [1, 7], 1]
        # θ_y sign convention consistent with B-matrix (XZ plane)
        N[:, [4, 10], 4] = -N[:, [5, 11], 5]
        dN_dξ[:, [4, 10], 4] = -dN_dξ[:, [5, 11], 5]

        # ----- Torsional Rotation (Linear Lagrange) -----
        N[:, [3, 9], 3] = N[:, [0, 6], 0]
        dN_dξ[:, [3, 9], 3] = dN_dξ[:, [0, 6], 0]

        return N, dN_dξ, d2N_dξ2

    def physical_coordinate_form(self, ξ: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Evaluate shape functions and derivatives in physical coordinates (x-space).

        Parameters
        ----------
        ξ : np.ndarray
            Natural coordinates ∈ [-1, 1] with shape (n_points,)

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
        Derivatives are transformed using:
        ∂N/∂x = (∂N/∂ξ)(∂ξ/∂x) = (∂N/∂ξ)(2/L)
        ∂²N/∂x² = (∂²N/∂ξ²)(∂ξ/∂x)² = (∂²N/∂ξ²)(4/L²)
        """
        N, dN_dξ, d2N_dξ2 = self.natural_coordinate_form(ξ)
        
        # Apply coordinate transforms
        dN_dx = dN_dξ * self.dξ_dx
        d2N_dx2 = d2N_dξ2 * self.d2ξ_dx2
        
        return N, dN_dx, d2N_dx2

    def natural_coordinate_form_coefficients(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Monomial coefficients for N, ∂N/∂ξ, and ∂²N/∂ξ² in ξ (same layout as B2 / ElementObject).

        Returns
        -------
        N_coefficients : np.ndarray
            Shape (12, 6, 4). ``N_coefficients[dof, comp, k]`` is the coefficient of ξ^k in N.
        dN_dxi_coefficients : np.ndarray
            Shape (12, 6, 4). Coefficients of ∂N/∂ξ as a polynomial in ξ.
        d2N_dxi2_coefficients : np.ndarray
            Shape (12, 6, 4). Coefficients of ∂²N/∂ξ² (zero for linear Lagrange DOFs).

        Notes
        -----
        Linear 2-node Timoshenko: axial, bending, and torsion use linear Lagrange in ξ;
        second derivatives in ξ vanish for those shape functions.
        """
        c = np.zeros((12, 6, 4), dtype=np.float64)
        dc = np.zeros((12, 6, 4), dtype=np.float64)
        d2c = np.zeros((12, 6, 4), dtype=np.float64)
        c[0, 0, 0], c[0, 0, 1] = 0.5, -0.5
        c[6, 0, 0], c[6, 0, 1] = 0.5, 0.5
        dc[0, 0, 0] = -0.5
        dc[6, 0, 0] = 0.5
        c[1, 1], c[7, 1] = c[0, 0].copy(), c[6, 0].copy()
        dc[1, 1], dc[7, 1] = dc[0, 0].copy(), dc[6, 0].copy()
        c[5, 5], c[11, 5] = c[0, 0].copy(), c[6, 0].copy()
        dc[5, 5], dc[11, 5] = dc[0, 0].copy(), dc[6, 0].copy()
        c[2, 2], c[8, 2] = c[1, 1].copy(), c[7, 1].copy()
        dc[2, 2], dc[8, 2] = dc[1, 1].copy(), dc[7, 1].copy()
        c[4, 4], c[10, 4] = -c[5, 5].copy(), -c[11, 5].copy()
        dc[4, 4], dc[10, 4] = -dc[5, 5].copy(), -dc[11, 5].copy()
        c[3, 3], c[9, 3] = c[0, 0].copy(), c[6, 0].copy()
        dc[3, 3], dc[9, 3] = dc[0, 0].copy(), dc[6, 0].copy()
        return c, dc, d2c

    def physical_coordinate_form_coefficients(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Monomial coefficients for N, ∂N/∂x, and ∂²N/∂x² expressed as polynomials in ξ.

        The polynomials are still parameterized by natural coordinate ξ ∈ [-1, 1];
        values satisfy ∂N/∂x = (∂N/∂ξ)(2/L) and ∂²N/∂x² = (∂²N/∂ξ²)(4/L²) pointwise.

        Returns
        -------
        N_coefficients : np.ndarray
            Same as :meth:`natural_coordinate_form_coefficients` (N is invariant).
        dN_dx_coefficients : np.ndarray
            Shape (12, 6, 4). Polynomial in ξ for ∂N/∂x at each (dof, component).
        d2N_dx2_coefficients : np.ndarray
            Shape (12, 6, 4). Polynomial in ξ for ∂²N/∂x².
        """
        cN, cdxi, cd2xi = self.natural_coordinate_form_coefficients()
        return cN, cdxi * self.dξ_dx, cd2xi * self.d2ξ_dx2

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

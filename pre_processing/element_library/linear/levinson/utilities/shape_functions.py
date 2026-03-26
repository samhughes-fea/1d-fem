# pre_processing\element_library\levinson\utilities\shape_functions.py
"""
Shape functions for 2-node 3-D Levinson beam (higher-order transverse fields).

**Tensors:** ``natural_coordinate_form(xi)`` returns ``N``, ``dN_dxi``, ``d2N_dxi2``, each shape ``(n_gp, 12, 6)`` —
row = element DOF index, column = component ``(u_x, u_y, u_z, theta_x, theta_y, theta_z)``.

**Weak form:** ``F_dist += w_g * N.T @ q * detJ`` with ``detJ = L/2``, ``xi in [-1, 1]``.

**Diff vs EB/Timoshenko (12 DOF contract):** same ``U_e`` length and ``N`` layout as standard beam elements; Levinson uses
quintic transverse / cubic bending rotations and feeds ``B_matrix`` with ``kappa_z`` before ``kappa_y`` in the strain vector.

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
    Evaluate Levinson beam shape functions and natural derivatives on ``xi``.

    Axial and torsion: linear Lagrange; transverse ``u_y``, ``u_z``: quintic with zero slope at nodes;
    ``theta_y``, ``theta_z``: cubic. ``physical_coordinate_form`` scales to ``dN/dx``, ``d2N/dx2`` via ``dxi_dx``, ``d2xi_dx2``.

    Parameters
    ----------
    element_length : float
        Chord length ``L`` (must be > 0).

    Attributes
    ----------
    dξ_dx : float
        ``2/L`` (first derivative chain rule factor).
    d2ξ_dx2 : float
        ``4/L**2`` (second derivative chain rule factor).

    Notes
    -----
    **Contract:** ``N`` batch shape ``(n_gp, 12, 6)`` matches standard beam load and stiffness assembly.
    **Diff:** Higher-order ``u`` / ``theta`` fields vs linear Timoshenko; shear correction factor not used in ``D`` (``G*A``);
    strain definitions and ``alpha`` terms live in ``B_matrix``.

    **Polynomial detail (natural coordinate xi):**

    - ``u_x``: ``N1_u = 0.5*(1-xi)``, ``N2_u = 0.5*(1+xi)``.
    - ``u_y``, ``u_z`` (quintic): ``N1_v = 1/2 - (15/16)*xi + (5/8)*xi**3 - (3/16)*xi**5``, etc.; zero slope at nodes.
    - ``theta_y``, ``theta_z`` (cubic): ``N1_theta = (1/4)*(2 - 3*xi + xi**3)``, ``N2_theta = (1/4)*(2 + 3*xi - xi**3)``.

    Strain linkage (in ``B_matrix``): ``gamma_xy = du_y/dx - theta_z + alpha*d2(theta_z)/dx2``, etc.;
    ``kappa_z = d(theta_z)/dx``, ``kappa_y = d(theta_y)/dx``, ``phi_x = d(theta_x)/dx``.

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
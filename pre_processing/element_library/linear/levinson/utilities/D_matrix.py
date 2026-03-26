# pre_processing/element_library/linear/levinson/utilities/D_matrix.py
"""Material stiffness ``D`` (6, 6) for Levinson beam. ``S = D @ eps``; diagonal ``EA``, ``EI_z``, ``EI_y``, ``G*A`` (shear, twice), ``GJ_t`` — no ``kappa`` factor (shear is ``G*A``).

Matches Voigt order used in ``levinson/utilities/B_matrix.py``. Parent: ``K_e += B.T @ D @ B * w_g * detJ``.
"""

import numpy as np
from typing import Dict
from dataclasses import dataclass, field

@dataclass(frozen=True)
class MaterialStiffnessOperator:
    """Constitutive operator (``D``) for 3D Levinson beam elements.

    Assembly and post-processing forms of the same 6x6 ``D``. Shear stiffness is ``G*A`` (no Timoshenko ``kappa``);
    higher-order shear in the kinematics is carried by ``B`` (``StrainDisplacementOperator``).

    Parameters
    ----------
    youngs_modulus : float
        Young's modulus (E) in Pascals (Pa)
    shear_modulus : float
        Shear modulus (G) in Pascals (Pa)
    cross_section_area : float
        Cross-sectional area (A) in m²
    moment_inertia_y : float
        Second moment of area about y-axis (I_y) in m⁴
    moment_inertia_z : float
        Second moment of area about z-axis (I_z) in m⁴
    torsion_constant : float
        Torsional constant (J_t) in m⁴
    warping_inertia_y : float, optional
        Warping constant about y-axis (I_wy) in m⁶, default=0
    warping_inertia_z : float, optional
        Warping constant about z-axis (I_wz) in m⁶, default=0

    This operator intentionally has **no** shear_correction_factor parameter;
    Levinson theory accounts for shear via higher-order terms in the
    strain-displacement relation (B-matrix), not via a constitutive κ.

    Attributes
    ----------
    has_warping_coupling : bool
        True if bending-torsion coupling exists (I_wy or I_wz ≠ 0)

    Notes
    -----
    Canonical `D` block (Levinson order, no shear correction factor):

    ```text
    S = D ε
    ε = [ε_x, κ_z, κ_y, γ_xy, γ_xz, φ_x]^T
    S = [N,   M_z, M_y, V_y,  V_z,  T  ]^T

    D =
    [ EA    0     0     0    0    0   ]
    [ 0    EI_z   0     0    0    0   ]
    [ 0     0    EI_y   0    0    0   ]
    [ 0     0     0     GA   0    0   ]
    [ 0     0     0     0    GA   0   ]
    [ 0     0     0     0    0   GJ_t ]
    ```

    **D tensor (shape (6, 6), Levinson Voigt order)**
    - row/col 0 ``eps_x``: ``D[0,0] = EA``.
    - row/col 1 ``kappa_z``: ``D[1,1] = EI_z``.
    - row/col 2 ``kappa_y``: ``D[2,2] = EI_y``.
    - row/col 3 ``gamma_xy``: ``D[3,3] = G*A`` (no ``kappa`` factor).
    - row/col 4 ``gamma_xz``: ``D[4,4] = G*A`` (no ``kappa`` factor).
    - row/col 5 ``phi_x``: ``D[5,5] = GJ_t``.
    - default off-diagonal entries are zero; optional warping coupling may fill
      ``D[1,5]``, ``D[5,1]``, ``D[2,5]``, ``D[5,2]``.

    **Resultant mapping**
    - ``S = D @ eps`` with ``eps = [eps_x, kappa_z, kappa_y, gamma_xy, gamma_xz, phi_x]``.
    - Rows of ``S`` are ``[N, M_z, M_y, V_y, V_z, T]``.

    **B/N linkage**
    - Parent weak form: ``K_e += B.T @ D @ B * w_g * detJ`` with selective bending/shear quadrature.
    - ``B`` and ``N`` tensors come from Levinson utilities with batch shape ``(n_gp, 12, 6)``.

    See Also
    --------
    linear_levinson_3D.LinearLevinsonBeamElement3D
    """

    # Material properties (immutable)
    youngs_modulus: float
    shear_modulus: float
    cross_section_area: float
    moment_inertia_y: float
    moment_inertia_z: float
    torsion_constant: float
    warping_inertia_y: float = 0.0
    warping_inertia_z: float = 0.0

    # Internal matrices
    _D_assembly: np.ndarray = field(init=False, repr=False)
    _D_postprocess: np.ndarray = field(init=False, repr=False)
    _energy_components: Dict[str, np.ndarray] = field(init=False, repr=False)

    def __post_init__(self):
        """Validate properties and build matrices immediately after construction."""
        self._validate_properties()
        self._build_constitutive_matrices()

    def assembly_form(self) -> np.ndarray:
        """
        Material matrix for ``K_e += B.T @ D @ B * w_g * detJ`` at each Gauss point.

        Returns
        -------
        np.ndarray
            Material stiffness matrix, shape (6, 6).
        """
        return self._D_assembly

    def postprocessing_form(self) -> np.ndarray:
        """
        Full section stiffness ``D`` for post-processing (same numbers as ``assembly_form``).

        Used for stress resultants ``S = D @ eps``, strain energy density, and verification.

        Returns
        -------
        np.ndarray
            Material stiffness matrix, shape (6, 6).
        """
        return self._D_postprocess

    def compute_stress_resultants(self, strain: np.ndarray) -> np.ndarray:
        """
        Stress resultants ``S = D @ strain`` (beam Voigt; Levinson strain row order).

        Parameters
        ----------
        strain : np.ndarray, shape (6,) or (6, n)
            Voigt strain ``[eps_x, kappa_z, kappa_y, gamma_xy, gamma_xz, phi_x]``
            (``kappa_z`` before ``kappa_y``, matching ``B_matrix`` / ``D`` row layout).

        Returns
        -------
        np.ndarray
            Same shape as ``strain``. Rows of ``S``: ``[N, M_z, M_y, V_y, V_z, T]``
            paired with those strain rows (``T`` torsional resultant, index 5).
        """
        return self.postprocessing_form() @ strain

    def energy_density_components(self, strain: np.ndarray) -> Dict[str, float]:
        """
        Decomposes strain energy density by deformation mode.
        
        Returns
        -------
        Dict[str, float]
            Components with keys:
            - 'total' : Total strain energy density
            - 'axial' : Axial deformation energy
            - 'bending_z' : Bending about z-axis energy
            - 'bending_y' : Bending about y-axis energy  
            - 'torsion' : Torsional energy
            - 'shear_xy' : Shear xy deformation energy
            - 'shear_xz' : Shear xz deformation energy
            - 'coupling_z' : Z-axis bending-torsion coupling energy
            - 'coupling_y' : Y-axis bending-torsion coupling energy
        """
        return {
            'total': 0.5 * strain.T @ self._D_postprocess @ strain,
            **{k: 0.5 * strain.T @ v @ strain 
               for k,v in self._energy_components.items()}
        }

    @property
    def has_warping_coupling(self) -> bool:
        """bool: True if non-zero warping constants induce bending-torsion coupling."""
        return (abs(self.warping_inertia_y) > 1e-12 or 
                abs(self.warping_inertia_z) > 1e-12)

    def _validate_properties(self) -> None:
        """Verify physical plausibility of all material parameters."""
        if not all(x > 0 for x in [
            self.youngs_modulus, self.shear_modulus,
            self.cross_section_area, self.moment_inertia_y,
            self.moment_inertia_z, self.torsion_constant
        ]):
            raise ValueError("All stiffness parameters must be strictly positive")

    def _build_constitutive_matrices(self) -> None:
        """Constructs and validates all constitutive matrices."""
        # Compute stiffness terms (consistent units)
        EA = self.youngs_modulus * self.cross_section_area
        EI_z = self.youngs_modulus * self.moment_inertia_z
        EI_y = self.youngs_modulus * self.moment_inertia_y
        GJ_t = self.shear_modulus * self.torsion_constant
        GA = self.shear_modulus * self.cross_section_area  # Levinson: no κ factor
        EIw_z = self.youngs_modulus * self.warping_inertia_z
        EIw_y = self.youngs_modulus * self.warping_inertia_y

        # Construct 6x6 matrix for Levinson (includes shear terms)
        D = np.zeros((6, 6), dtype=np.float64)
        D[0, 0] = EA          # Axial
        D[1, 1] = EI_z        # Bending-Z
        D[2, 2] = EI_y        # Bending-Y
        D[3, 3] = GA          # Shear xy (Levinson: no κ)
        D[4, 4] = GA          # Shear xz (Levinson: no κ)
        D[5, 5] = GJ_t        # Torsion
        
        # Warping coupling terms (if present)
        if self.has_warping_coupling:
            D[1, 5] = -EIw_z  # Bending-Z / Torsion coupling
            D[5, 1] = -EIw_z
            D[2, 5] = EIw_y   # Bending-Y / Torsion coupling
            D[5, 2] = EIw_y

        object.__setattr__(self, '_D_assembly', D)
        object.__setattr__(self, '_D_postprocess', D.copy())

        if self.has_warping_coupling and not np.allclose(D, D.T, atol=1e-12):
            raise ValueError("Warping terms violate D-matrix symmetry")

        object.__setattr__(self, '_energy_components', {
            'axial': np.diag([EA, 0, 0, 0, 0, 0]),
            'bending_z': np.diag([0, EI_z, 0, 0, 0, 0]),
            'bending_y': np.diag([0, 0, EI_y, 0, 0, 0]),
            'torsion': np.diag([0, 0, 0, 0, 0, GJ_t]),
            'shear_xy': np.diag([0, 0, 0, GA, 0, 0]),
            'shear_xz': np.diag([0, 0, 0, 0, GA, 0]),
            'coupling_z': np.zeros((6, 6)) if not self.has_warping_coupling else np.array([
                [0,0,0,0,0,0], [0,0,0,0,0,-EIw_z], [0,0,0,0,0,0], 
                [0,0,0,0,0,0], [0,0,0,0,0,0], [0,-EIw_z,0,0,0,0]
            ]),
            'coupling_y': np.zeros((6, 6)) if not self.has_warping_coupling else np.array([
                [0,0,0,0,0,0], [0,0,0,0,0,0], [0,0,0,0,0,EIw_y], 
                [0,0,0,0,0,0], [0,0,0,0,0,0], [0,0,EIw_y,0,0,0]
            ])
        })
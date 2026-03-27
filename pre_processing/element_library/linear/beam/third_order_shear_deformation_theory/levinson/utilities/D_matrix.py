# pre_processing/element_library/linear/levinson/utilities/D_matrix.py
"""Material stiffness D (6, 6) for Levinson beam.

S = D ε with diagonal EA, EI_z, EI_y, GA (shear, twice), and GJ_t.
No shear-correction factor is used in D.
Voigt order matches `levinson/utilities/B_matrix.py`.
Parent assembly uses `K_e += B.T @ D @ B * w_g * detJ`.
"""

import numpy as np
from typing import Dict
from dataclasses import dataclass, field

@dataclass(frozen=True)
class MaterialStiffnessOperator:
    """
    Constitutive tensor D ∈ ℝ^{6×6} for 3-D Levinson beam elements.

    D is a rank-2 symmetric tensor relating the generalised strain vector
    ε ∈ ℝ^6 to the beam section resultant vector S ∈ ℝ^6 via S = D ε.

    **Note on Voigt row order (Levinson)**

    Levinson uses a non-standard row order with κ_z before κ_y, matching B:

        ε = [ε_x,  κ_z,  κ_y,  γ_xy, γ_xz, φ_x]^T   (Voigt strains)
        S = [N,    M_z,  M_y,  V_y,  V_z,  T  ]^T   (section resultants)

    D has non-zero shear entries D[3,3] = D[4,4] = G·A (no shear correction factor κ).
    The higher-order correction to shear enters kinematically through B via the
    α ∂²θ/∂x² terms; D itself is constitutively uncorrected.

    Parameters
    ----------
    youngs_modulus : float
        Young's modulus E [Pa].
    shear_modulus : float
        Shear modulus G [Pa].
    cross_section_area : float
        Cross-sectional area A [m²].
    moment_inertia_y : float
        Second moment of area about y, I_y [m⁴].
    moment_inertia_z : float
        Second moment of area about z, I_z [m⁴].
    torsion_constant : float
        St. Venant torsional constant J_t [m⁴].
    warping_inertia_y : float, optional
        Warping constant about y, I_wy [m⁶], default = 0.
    warping_inertia_z : float, optional
        Warping constant about z, I_wz [m⁶], default = 0.

    Attributes
    ----------
    has_warping_coupling : bool
        True if bending-torsion coupling exists (I_wy or I_wz ≠ 0).

    Notes
    -----
    **Sparsity structure of D (Levinson order, no shear correction factor)**

    ```text
    S = D ε
    ε = [ε_x,  κ_z,  κ_y,  γ_xy, γ_xz, φ_x]^T
    S = [N,    M_z,  M_y,  V_y,  V_z,  T  ]^T

    D =
    [ EA    0     0     0    0    0   ]
    [ 0    EI_z   0     0    0    0   ]
    [ 0     0    EI_y   0    0    0   ]
    [ 0     0     0     GA   0    0   ]
    [ 0     0     0     0    GA   0   ]
    [ 0     0     0     0    0   GJ_t ]
    ```

    **Component definitions — D ∈ ℝ^{6×6}, diagonal, rank-2**

    ```text
    D[0,0] = E·A      (axial stiffness)
    D[1,1] = E·I_z    (bending stiffness about z; κ_z is row 1 in Levinson order)
    D[2,2] = E·I_y    (bending stiffness about y; κ_y is row 2 in Levinson order)
    D[3,3] = G·A      (shear stiffness XY; no κ factor — higher-order correction in B)
    D[4,4] = G·A      (shear stiffness XZ; no κ factor — higher-order correction in B)
    D[5,5] = G·J_t    (St. Venant torsional stiffness)
    D[i,j] = 0  for all i ≠ j   (default; warping coupling may fill D[1,5], D[2,5] etc.)
    ```

    Optional warping coupling: when warping inertia I_wy or I_wz is non-zero,
    off-diagonal entries D[1,5], D[5,1] and D[2,5], D[5,2] become non-zero.

    **Weak-form assembly linkage**

    The element stiffness is accumulated as `K_e += B.T @ D @ B * w_g * detJ`
    with ξ ∈ [−1, 1] and `detJ = L/2`. B ∈ ℝ^{6×12} comes from
    `levinson/utilities/B_matrix.py`; the shape-function tensors
    `N`, `dN_dxi`, `d2N_dxi2` of batch shape (n_gp, 12, 6) come from
    `shape_functions.py`. Selective integration is applied on shear rows.

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
            Voigt strain [ε_x, κ_z, κ_y, γ_xy, γ_xz, φ_x]
            (κ_z before κ_y, matching B_matrix / D row layout).

        Returns
        -------
        np.ndarray
            Same shape as ``strain``. Rows of S: [N, M_z, M_y, V_y, V_z, T]
            paired with those strain rows (T: torsional resultant, index 5).
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
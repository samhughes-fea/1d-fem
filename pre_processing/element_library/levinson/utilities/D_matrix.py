# pre_processing\element_library\levinson\utilities\D_matrix.py

import numpy as np
from typing import Dict
from dataclasses import dataclass, field

@dataclass(frozen=True)
class MaterialStiffnessOperator:
    """Constitutive operator for 3D Levinson beam elements.
    
    Encapsulates the material stiffness matrix (D-matrix) with dual representations:
    - Assembly form: Optimized for stiffness matrix assembly (Kᵉ = ∫BᵀDB dx)
    - Postprocessing form: Complete form for stress/strain computation and energy decomposition

    Mathematical Formulation
    -----------------------
    Levinson beam theory includes shear deformation. Higher-order kinematics
    (e.g. α(∂²θ/∂x²) in the shear strain in the B-matrix) give a better
    approximation of shear stress distribution over the section, so no
    empirical shear correction factor κ is needed; the constitutive shear
    stiffness is **GA** (κ = 1 implicitly), unlike Timoshenko which uses κGA.

    ⎡ N  ⎤   ⎡ EA     0       0       0       0       0   ⎤ ⎡ ε_x  ⎤
    ⎢ M_z⎥   ⎢ 0     EI_z     0       0       0       0   ⎥ ⎢ κ_z  ⎥
    ⎢ M_y⎥ = ⎢ 0      0     EI_y      0       0       0   ⎥ ⎢ κ_y  ⎥
    ⎢ V_y⎥   ⎢ 0      0       0      GA       0       0   ⎥ ⎢ γ_xy ⎥
    ⎢ V_z⎥   ⎢ 0      0       0       0      GA       0   ⎥ ⎢ γ_xz ⎥
    ⎣ M_x⎦   ⎣ 0      0       0       0       0    GJ_t ⎦ ⎣ φ_x  ⎦

    Symbol definitions:
        EA   = E · A       (axial stiffness)
        EI_y = E · I_y     (bending about y)
        EI_z = E · I_z     (bending about z)
        GA   = G · A       (shear stiffness; no κ, unlike Timoshenko)
        GJ_t = G · J_t     (torsional stiffness)

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
        Retrieves the material matrix optimized for stiffness matrix assembly.
        
        Used in the computation of Kᵉ = ∫BᵀDB dx where:
        - B is the strain-displacement matrix
        - D is this material stiffness matrix
        - Integration is performed over element domain

        Returns
        -------
        np.ndarray
            6×6 material stiffness matrix in assembly-optimized form
        """
        return self._D_assembly

    def postprocessing_form(self) -> np.ndarray:
        """
        Retrieves the complete material matrix for analysis and visualization.
        
        Used for:
        - Stress recovery (σ = Dε)
        - Strain energy calculations
        - Result verification and postprocessing

        Returns
        -------
        np.ndarray 
            6×6 material stiffness matrix in complete form
        """
        return self._D_postprocess

    def compute_stress_resultants(self, strain: np.ndarray) -> np.ndarray:
        """
        Compute stress resultants from strain measures using full constitutive relation.
        
        Parameters
        ----------
        strain : np.ndarray, shape (6,) or (6,n)
            Strain vector/matrix in Voigt notation [ε_x, κ_z, κ_y, γ_xy, γ_xz, φ_x]

        Returns
        -------
        np.ndarray
            Stress resultants [N, M_z, M_y, V_y, V_z, M_x] in same shape as input
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
# pre_processing/element_library/linear/timoshenko/utilities/D_matrix.py
"""Material stiffness ``D`` (6, 6) for Timoshenko beam. ``S = D @ eps``; diagonal ``EA``, ``EI_y``, ``EI_z``, ``kappa*G*A`` (shear, twice), ``GJ_t``.

Parent element: ``K_e += B.T @ D @ B * w_g * detJ`` (with selective bending/shear quadrature in ``linear_timoshenko_3D.py``).
"""

import numpy as np
from typing import Dict
from dataclasses import dataclass, field

@dataclass(frozen=True)
class MaterialStiffnessOperator:
    """
    Constitutive operator (``D`` matrix) for 3-D Timoshenko beam elements.

    Stores assembly and post-processing copies of the same 6x6 ``D``.

    Parameters
    ----------
    youngs_modulus : float
        Young's modulus *E* [Pa].
    shear_modulus : float
        Shear modulus *G* [Pa].
    cross_section_area : float
        Cross-sectional area *A* [m²].
    moment_inertia_y : float
        Second moment of area about **y** (*I_y*) [m⁴].
    moment_inertia_z : float
        Second moment of area about **z** (*I_z*) [m⁴].
    torsion_constant : float
        Torsional constant *J_t* [m⁴].
    shear_correction_factor : float, optional
        Shear correction factor κ (default: 5/6 for rectangular sections).
        Accounts for non-uniform shear stress distribution across the cross-section.

    Attributes
    ----------
    _D_assembly : ndarray (6 × 6)
        Sparse-by-design matrix used inside the element stiffness loop.
    _D_postprocess : ndarray (6 × 6)
        Copy of *D* kept intact for stress / energy work.
    _energy_components : dict[str, ndarray]
        Pre-factored diagonal blocks for axial, bending-y, bending-z,
        torsion, shear-xy and shear-xz.

    Notes
    -----
    Voigt: ``eps`` and ``S`` order per ``FORMULATION_DOCSTRING_STANDARDS.md``. ``D`` is diagonal with
    ``EA``, ``EI_y``, ``EI_z``, ``kappa*G*A`` on rows 3 and 4 (shear), ``GJ_t`` for torsion.
    Weak form: ``K_e += B.T @ D @ B * w_g * detJ`` on ``xi in [-1, 1]``; parent may use different Gauss rules for bending vs shear blocks.

    See Also
    --------
    linear_timoshenko_3D.LinearTimoshenkoBeamElement3D
    docs/conventions/FORMULATION_DOCSTRING_STANDARDS.md

    """

    # --- Material properties (immutable) ---
    youngs_modulus: float
    shear_modulus: float
    cross_section_area: float
    moment_inertia_y: float
    moment_inertia_z: float
    torsion_constant: float
    shear_correction_factor: float = 5.0 / 6.0  # Default κ = 5/6 for rectangular sections
    y_sc: float = 0.0  # Shear centre offset from centroid (y) [m]
    z_sc: float = 0.0  # Shear centre offset from centroid (z) [m]

    # --- Internal container matrices --- 

    _D_assembly: np.ndarray = field(init=False, repr=False)
    _D_postprocess: np.ndarray = field(init=False, repr=False)
    _energy_components: Dict[str, np.ndarray] = field(init=False, repr=False)

    def __post_init__(self):
        """Validate properties and build internal matrices immediately after construction."""
        self._validate_properties()
        self._build_constitutive_matrices()

    def _validate_properties(self) -> None:
        if not all(x > 0 for x in (
            self.youngs_modulus,
            self.shear_modulus,
            self.cross_section_area,
            self.moment_inertia_y,
            self.moment_inertia_z,
            self.torsion_constant,
            self.shear_correction_factor,
        )):
            raise ValueError("All stiffness parameters must be strictly positive")
        # y_sc, z_sc may be any finite float (zero or non-zero)

    def assembly_form(self) -> np.ndarray:
        """
        Material matrix for ``K_e += B.T @ D @ B * w_g * detJ`` at each Gauss point.

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
            Strain vector/matrix in Voigt notation ε = [ εₓ  κᵧ  κ_z  γ_xy  γ_xz  φₓ ]ᵀ

        Returns
        -------
        np.ndarray`
            Stress resultants [N, M_z, M_y, M_x] in same shape as input
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
        """
        return {
            'total': 0.5 * strain.T @ self._D_postprocess @ strain,
            **{k: 0.5 * strain.T @ v @ strain 
               for k,v in self._energy_components.items()}
        }

    def _build_constitutive_matrices(self) -> None:
        r"""Constructs and validates all constitutive matrices.
        
        | Voigt component | Meaning                    | Row/col in D | Filled with |
        | --------------- | -------------------------- | -------------| ----------- |
        | 0               | ε_x  (axial)               | 0            | **EA**      |
        | 1               | κ_y (bending about y)      | 1            | **EI_y**    |
        | 2               | κ_z (bending about z)      | 2            | **EI_z**    |
        | 3               | γ_xy (shear-xy)            | 3            | **κGA**     |
        | 4               | γ_xz (shear-xz)            | 4            | **κGA**     |
        | 5               | φ_x  (torsion)             | 5            | **GJ_t**    |

        """
        # Compute stiffness terms (consistent units)
        EA = self.youngs_modulus * self.cross_section_area
        EI_z = self.youngs_modulus * self.moment_inertia_z
        EI_y = self.youngs_modulus * self.moment_inertia_y
        GJ_t = self.shear_modulus * self.torsion_constant
        kappa_GA = self.shear_correction_factor * self.shear_modulus * self.cross_section_area

        # --- main constitutive matrix ---
        D = np.zeros((6, 6), dtype=np.float64)
        D[0, 0] = EA
        D[1, 1] = EI_y          # bending about y
        D[2, 2] = EI_z          # bending about z
        D[3, 3] = kappa_GA      # shear xy (Timoshenko)
        D[4, 4] = kappa_GA      # shear xz (Timoshenko)
        D[5, 5] = GJ_t          # torsion
        # Shear-centre coupling: when y_sc or z_sc != 0, bending–torsion coupling
        if self.y_sc != 0.0 or self.z_sc != 0.0:
            D[1, 5] = D[5, 1] = -EI_y * self.z_sc
            D[2, 5] = D[5, 2] = EI_z * self.y_sc

        object.__setattr__(self, '_D_assembly', D)
        object.__setattr__(self, '_D_postprocess', D.copy())

        object.__setattr__(self, '_energy_components', {
            'axial'     : np.diag([EA, 0, 0, 0, 0, 0]),
            'bending_y' : np.diag([0, EI_y, 0, 0, 0, 0]),
            'bending_z' : np.diag([0, 0, EI_z, 0, 0, 0]),
            'torsion'   : np.diag([0, 0, 0, 0, 0, GJ_t]),
            'shear_xy'  : np.diag([0, 0, 0, kappa_GA, 0, 0]),
            'shear_xz'  : np.diag([0, 0, 0, 0, kappa_GA, 0]),
        })
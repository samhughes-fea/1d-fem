# pre_processing/element_library/linear/timoshenko/utilities/D_matrix.py
"""Material stiffness **D** (6, 6) — **linear elastic**, infinitesimal strain (**Cauchy** resultants).

**Constitutive (Voigt):** \\(\\mathbf{S} = \\mathbf{D}\\,\\boldsymbol{\\varepsilon}\\), diagonal \\(EA\\), \\(EI_y\\), \\(EI_z\\),
\\(\\kappa GA\\) (shear, twice), \\(GJ_t\\).

**Stiffness:** \\(\\mathbf{K} = \\int \\mathbf{B}^\\top \\mathbf{D}\\,\\mathbf{B}\\,\\mathrm{d}x\\) — **constant** w.r.t. \\(\\mathbf{U}_e\\).
Assembly: `K_e += B.T @ D @ B * w_g * detJ`; selective orders in `linear_timoshenko_3D.py`.
"""

import numpy as np
from typing import Dict
from dataclasses import dataclass, field

@dataclass(frozen=True)
class MaterialStiffnessOperator:
    """
    Constitutive tensor **D** ∈ ℝ^{6×6} for 3-D Timoshenko beam elements (**linear isotropic**, infinitesimal strain).

    **Governing relation:** \\(\\mathbf{S} = \\mathbf{D}\\,\\boldsymbol{\\varepsilon}\\) with \\(\\boldsymbol{\\varepsilon}\\) from `B_matrix.py`.

    D is a rank-2 symmetric tensor relating the generalised strain vector
    ε ∈ ℝ^6 to the beam section resultant vector S ∈ ℝ^6 via S = D ε:

        ε = [ε_x,  κ_y,  κ_z,  γ_xy, γ_xz, φ_x]^T   (Voigt strains)
        S = [N,    M_y,  M_z,  V_y,  V_z,  T  ]^T   (section resultants)

    Unlike the Euler-Bernoulli theory, the shear diagonal entries D[3,3] and
    D[4,4] are non-zero; they are set to κ·G·A where κ is the shear correction
    factor that accounts for the non-uniform shear stress distribution across the
    cross-section (κ = 5/6 for rectangular sections by default). Stores assembly
    and post-processing copies of the same 6×6 D.

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
        Saint-Venant torsion constant J_t [m⁴], a geometric constant; ``D`` holds the
        torsional stiffness ``G J_t`` [N m²]. Not the polar moment ``I_p``.
    shear_correction_factor : float, optional
        Shear correction factor κ (default: 5/6 for rectangular sections).

    Attributes
    ----------
    _D_assembly : ndarray, shape (6, 6)
        Sparse-by-design matrix used inside the element stiffness loop.
    _D_postprocess : ndarray, shape (6, 6)
        Copy of D kept intact for stress/energy post-processing.
    _energy_components : dict[str, ndarray]
        Pre-factored diagonal blocks for each deformation mode.

    Notes
    -----
    **Sparsity structure of D (Timoshenko, Voigt order)**

    ```text
    S = D ε
    ε = [ε_x,  κ_y,  κ_z,  γ_xy, γ_xz, φ_x]^T
    S = [N,    M_y,  M_z,  V_y,  V_z,  T  ]^T

    D =
    [ EA    0     0     0    0    0   ]
    [ 0    EI_y   0     0    0    0   ]
    [ 0     0    EI_z   0    0    0   ]
    [ 0     0     0    κGA   0    0   ]
    [ 0     0     0     0   κGA   0   ]
    [ 0     0     0     0    0   GJ_t ]
    ```

    **Component definitions — D ∈ ℝ^{6×6}, diagonal, rank-2**

    ```text
    D[0,0] = E·A      (axial stiffness)
    D[1,1] = E·I_y    (bending stiffness about y)
    D[2,2] = E·I_z    (bending stiffness about z)
    D[3,3] = κ·G·A    (shear stiffness in XY plane; κ is shear correction factor)
    D[4,4] = κ·G·A    (shear stiffness in XZ plane; κ is shear correction factor)
    D[5,5] = G·J_t    (St. Venant torsional stiffness)
    D[i,j] = 0  for all other i ≠ j
    ```

    The numbers come from the section provider (``SectionAnalysis.from_isotropic``), so
    this operator and the element ``section_analysis`` hold the same ``D``. A shear centre
    ``(y_sc, z_sc)`` [m] off the centroid couples shear and twist about it::

        D[3,5] = -z_sc κGA     D[4,5] = y_sc κGA     D[5,5] = G J_t + κGA (y_sc² + z_sc²)

    so a shear force through the shear centre twists nothing and the torsional stiffness
    about the shear centre is ``G J_t``. Axial and bending stay uncoupled. (Before the
    provider took over, a non-zero offset added ``D[1,5] = -EI_y z_sc`` and
    ``D[2,5] = EI_z y_sc``, which are not in units of a stiffness.)

    **Weak-form assembly linkage**

    The element stiffness is accumulated as `K_e += B.T @ D @ B * w_g * detJ`
    with ξ ∈ [−1, 1] and `detJ = L/2`. B ∈ ℝ^{6×12} comes from
    `timoshenko/utilities/B_matrix.py`; the shape-function tensors
    `N`, `dN_dxi` of batch shape (n_gp, 12, 6) come from `shape_functions.py`.
    Selective reduced integration on shear rows avoids shear locking.

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
        Stress resultants ``S = D @ strain`` (beam Voigt, not 3D Cauchy stress).

        Parameters
        ----------
        strain : np.ndarray, shape (6,) or (6, n)
            Voigt strain ``[eps_x, kappa_y, kappa_z, gamma_xy, gamma_xz, phi_x]``
            (same row order as ``FORMULATION_DOCSTRING_STANDARDS.md``).

        Returns
        -------
        np.ndarray
            Same shape as ``strain``. Rows of ``S``: ``[N, M_y, M_z, V_y, V_z, T]``
            paired with strain rows.
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
        # The section D comes from the section provider (``SectionAnalysis.from_isotropic``);
        # 1D-FEM computes no section properties of its own.
        from onedfem.pre_processing.element_library.isotropic_section import isotropic_stiffness

        D = isotropic_stiffness(
            E=self.youngs_modulus, G=self.shear_modulus, A=self.cross_section_area,
            I_y=self.moment_inertia_y, I_z=self.moment_inertia_z, J_t=self.torsion_constant,
            kappa=self.shear_correction_factor, y_sc=self.y_sc, z_sc=self.z_sc,
        )
        EA, EI_y, EI_z = D[0, 0], D[1, 1], D[2, 2]
        kappa_GA, GJ_t = D[3, 3], D[5, 5]

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
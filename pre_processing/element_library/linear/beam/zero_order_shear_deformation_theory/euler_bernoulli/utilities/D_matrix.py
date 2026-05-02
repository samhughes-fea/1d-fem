# pre_processing/element_library/linear/euler_bernoulli/utilities/D_matrix.py
"""Material stiffness **D** (6, 6) — **linear elastic**, infinitesimal strain (**Cauchy** resultants).

**Constitutive (Voigt):** \\(\\mathbf{S} = \\mathbf{D}\\,\\boldsymbol{\\varepsilon}\\) with \\(\\boldsymbol{\\varepsilon}\\) from `B_matrix.py`
(\\(\\mathbf{S} = [N, M_y, M_z, V_y, V_z, T]^T\\)). Diagonal **D**: \\(EA\\), \\(EI_y\\), \\(EI_z\\), zero shear rows, \\(GJ_t\\).

**Stiffness:** contributes \\(\\mathbf{K} = \\int \\mathbf{B}^\\top \\mathbf{D}\\,\\mathbf{B}\\,\\mathrm{d}x\\) — **constant** w.r.t. \\(\\mathbf{U}_e\\).
Assembly: `K_e += B.T @ D @ B * w_g * detJ` (`linear_euler_bernoulli_3D.py`).
"""

import numpy as np
from typing import Dict
from dataclasses import dataclass, field

@dataclass(frozen=True)
class MaterialStiffnessOperator:
    """
    Constitutive tensor **D** ∈ ℝ^{6×6} for 3-D Euler-Bernoulli beam elements (**linear isotropic**, infinitesimal strain).

    **Governing relation:** \\(\\mathbf{S} = \\mathbf{D}\\,\\boldsymbol{\\varepsilon}\\) with \\(\\boldsymbol{\\varepsilon}\\) the infinitesimal engineering strain vector (same Voigt as TL nonlinear uses for \\(\\mathbf{E}_\\mathrm{lin}\\), but here \\(\\boldsymbol{\\varepsilon}\\) is strictly linear in \\(\\mathbf{U}_e\\)).

    D is a rank-2 symmetric tensor relating the generalised strain vector
    ε ∈ ℝ^6 to the beam section resultant vector S ∈ ℝ^6 via S = D ε:

        ε = [ε_x,  κ_y,  κ_z,  γ_xy, γ_xz, φ_x]^T   (Voigt strains)
        S = [N,    M_y,  M_z,  V_y,  V_z,  T  ]^T   (section resultants)

    D is diagonal; shear rows D[3,:] and D[4,:] are identically zero because
    the EB hypothesis imposes γ_xy = γ_xz = 0 kinematically. Shear forces
    V_y and V_z are recovered from equilibrium V = dM/dx. Stores assembly and
    post-processing copies (identical for EB; split supports other theories).

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
    **Sparsity structure of D (Euler-Bernoulli, Voigt order)**

    ```text
    S = D ε
    ε = [ε_x,  κ_y,  κ_z,  γ_xy, γ_xz, φ_x]^T
    S = [N,    M_y,  M_z,  V_y,  V_z,  T  ]^T

    D =
    [ EA    0     0     0    0    0   ]
    [ 0    EI_y   0     0    0    0   ]
    [ 0     0    EI_z   0    0    0   ]
    [ 0     0     0     0    0    0   ]
    [ 0     0     0     0    0    0   ]
    [ 0     0     0     0    0   GJ_t ]
    ```

    **Component definitions — D ∈ ℝ^{6×6}, diagonal, rank-2**

    ```text
    D[0,0] = E·A      (axial stiffness)
    D[1,1] = E·I_y    (bending stiffness about y)
    D[2,2] = E·I_z    (bending stiffness about z)
    D[3,3] = 0        (no constitutive shear; γ_xy = 0 by EB kinematic constraint)
    D[4,4] = 0        (no constitutive shear; γ_xz = 0 by EB kinematic constraint)
    D[5,5] = G·J_t    (St. Venant torsional stiffness)
    D[i,j] = 0  for all i ≠ j   (uncoupled; no shear-centre offset terms)
    ```

    **Weak-form assembly linkage**

    The element stiffness is accumulated as `K_e += B.T @ D @ B * w_g * detJ`
    with ξ ∈ [−1, 1] and `detJ = L/2`. B ∈ ℝ^{6×12} comes from
    `euler_bernoulli/utilities/B_matrix.py`; the shape-function tensors
    `N`, `dN_dxi`, `d2N_dxi2` of batch shape (n_gp, 12, 6) come from
    `shape_functions.py`.

    See Also
    --------
    linear_euler_bernoulli_3D.LinearEulerBernoulliBeamElement3D
    docs/conventions/FORMULATION_DOCSTRING_STANDARDS.md

    """

    # --- Material properties (immutable) ---
    youngs_modulus: float
    shear_modulus: float
    cross_section_area: float
    moment_inertia_y: float
    moment_inertia_z: float
    torsion_constant: float

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
        )):
            raise ValueError("All stiffness parameters must be strictly positive")

    def assembly_form(self) -> np.ndarray:
        """
        Material matrix for stiffness assembly: ``K_e += B.T @ D @ B * w_g * detJ`` at each Gauss point.

        Returns
        -------
        np.ndarray
            Material stiffness matrix, shape (6, 6).
        """
        return self._D_assembly

    def postprocessing_form(self) -> np.ndarray:
        """
        Full section stiffness ``D`` for post-processing (same numbers as ``assembly_form`` for EB).

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
            Voigt strain vector [ε_x, κ_y, κ_z, γ_xy, γ_xz, φ_x]
            (EB: γ_xy = γ_xz = 0).

        Returns
        -------
        np.ndarray
            Same shape as ``strain``. Rows of S: [N, M_y, M_z, V_y, V_z, T]
            paired with strain rows (EB: V_y, V_z from D are zero).
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
        """
        return {
            'total': 0.5 * strain.T @ self._D_postprocess @ strain,
            **{k: 0.5 * strain.T @ v @ strain 
               for k,v in self._energy_components.items()}
        }

    def _build_constitutive_matrices(self) -> None:
        """Constructs and validates all constitutive matrices.
        
        | Voigt component | Meaning                    | Row/col in D | Filled with |
        | --------------- | -------------------------- | -------------| ----------- |
        | 0               | ε_x  (axial)               | 0            | **EA**      |
        | 1               | κ_y (bending about y)      | 1            | **EI_y**    |
        | 2               | κ_z (bending about z)      | 2            | **EI_z**    |
        | 3               | γ_xy (shear-xy)            | 3            | 0           |
        | 4               | γ_xz (shear-xz)            | 4            | 0           |
        | 5               | φ_x  (torsion)             | 5            | **GJ_t**    |

        """
        # Compute stiffness terms (consistent units)
        EA = self.youngs_modulus * self.cross_section_area
        EI_z = self.youngs_modulus * self.moment_inertia_z
        EI_y = self.youngs_modulus * self.moment_inertia_y
        GJ_t = self.shear_modulus * self.torsion_constant

        # --- main constitutive matrix ---
        D = np.zeros((6, 6), dtype=np.float64)
        D[0, 0] = EA
        D[1, 1] = EI_y          # bending about y
        D[2, 2] = EI_z          # bending about z
        D[5, 5] = GJ_t          # torsion

        object.__setattr__(self, '_D_assembly', D)
        object.__setattr__(self, '_D_postprocess', D.copy())

        object.__setattr__(self, '_energy_components', {
            'axial'     : np.diag([EA, 0, 0, 0, 0, 0]),
            'bending_y' : np.diag([0, EI_y, 0, 0, 0, 0]),
            'bending_z' : np.diag([0, 0, EI_z, 0, 0, 0]),
            'torsion'   : np.diag([0, 0, 0, 0, 0, GJ_t]),
            'shear_xy'  : np.zeros((6, 6)),
            'shear_xz'  : np.zeros((6, 6)),
        })
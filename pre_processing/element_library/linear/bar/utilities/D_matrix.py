# pre_processing/element_library/linear/bar/utilities/D_matrix.py
"""
Material stiffness ``D`` (2, 2) for bar: ``diag(EA, GJ_t)``; ``S_section = D @ eps`` (axial force, torque).

Used in ``K_e += B.T @ D @ B * w_g * detJ`` in ``linear_bar_3D.py``. Optional 6x6 stress view zeros out unused rows for post-processing.
"""

import numpy as np
from typing import Dict
from dataclasses import dataclass, field


@dataclass(frozen=True)
class MaterialStiffnessOperator:
    """
    Constitutive operator for 3-D Bar elements (axial + torsion only).

    Strain: ε = [ ε_axial  φ_torsion ]ᵀ
    Stress resultants: N = [ N_axial  M_torsion ]ᵀ = D @ ε

    D = diag(EA, GJ_t)
    - EA = Young's modulus × area (axial stiffness)
    - GJ_t = Shear modulus × torsion constant (torsional stiffness)

    Parameters
    ----------
    youngs_modulus : float
        Young's modulus E [Pa].
    shear_modulus : float
        Shear modulus G [Pa].
    cross_section_area : float
        Cross-sectional area A [m²].
    torsion_constant : float
        Torsional constant J_t [m⁴].

    Notes
    -----
    Weak-form linkage: same Gauss accumulation as other 1-D elements; ``detJ = L/2`` on ``xi in [-1, 1]``.

    See Also
    --------
    linear_bar_3D.LinearBarElement3D
    """

    youngs_modulus: float
    shear_modulus: float
    cross_section_area: float
    torsion_constant: float

    _D_assembly: np.ndarray = field(init=False, repr=False)
    _D_postprocess: np.ndarray = field(init=False, repr=False)
    _energy_components: Dict[str, np.ndarray] = field(init=False, repr=False)

    def __post_init__(self):
        self._validate_properties()
        self._build_constitutive_matrices()

    def _validate_properties(self) -> None:
        if not all(x > 0 for x in (
            self.youngs_modulus,
            self.shear_modulus,
            self.cross_section_area,
            self.torsion_constant,
        )):
            raise ValueError("All stiffness parameters must be strictly positive")

    def _build_constitutive_matrices(self) -> None:
        EA = self.youngs_modulus * self.cross_section_area
        GJ_t = self.shear_modulus * self.torsion_constant
        D = np.diag([EA, GJ_t]).astype(np.float64)
        object.__setattr__(self, '_D_assembly', D)
        object.__setattr__(self, '_D_postprocess', D.copy())
        object.__setattr__(self, '_energy_components', {
            'axial': np.diag([EA, 0.0]),
            'torsion': np.diag([0.0, GJ_t]),
        })

    def assembly_form(self) -> np.ndarray:
        """Return 2×2 D-matrix for stiffness assembly."""
        return self._D_assembly

    def postprocessing_form(self) -> np.ndarray:
        """Return 2×2 D-matrix for stress recovery and post-processing."""
        return self._D_postprocess

    def compute_stress_resultants(self, strain: np.ndarray) -> np.ndarray:
        """Compute stress resultants N = D @ ε. strain shape (2,) or (2, n)."""
        return self.postprocessing_form() @ np.asarray(strain, dtype=np.float64)

    def energy_density_components(self, strain: np.ndarray) -> Dict[str, float]:
        """Strain energy density breakdown: total, axial, torsion. strain shape (2,)."""
        e = np.asarray(strain, dtype=np.float64).reshape(2)
        out = {'total': float(0.5 * e.T @ self._D_postprocess @ e)}
        for k, v in self._energy_components.items():
            out[k] = float(0.5 * e.T @ v @ e)
        return out

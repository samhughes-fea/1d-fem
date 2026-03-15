# pre_processing/element_library/linear/truss/utilities/D_matrix.py
"""Material stiffness D for 2-node 3-D Truss. N = D @ ε; D shape (3, 3), diag(EA, κGA, GJ_t). Full 6-component view (N, M_y, M_z, V_y, V_z, T) uses zeros for bending/shear not modelled."""

import numpy as np
from typing import Dict
from dataclasses import dataclass, field


@dataclass(frozen=True)
class MaterialStiffnessOperator:
    """
    Constitutive operator for 3-D Truss elements (axial + transverse shear + torsion).

    Strain: ε = [ ε_axial  γ_transverse  φ_torsion ]ᵀ
    Stress resultants: N = D @ ε

    D = diag(EA, κGA, GJ_t)
    - EA = axial stiffness
    - κGA = shear correction × G × A (transverse)
    - GJ_t = torsional stiffness

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
    shear_correction_factor : float, default 5/6
        κ for transverse shear (κGA).
    """

    youngs_modulus: float
    shear_modulus: float
    cross_section_area: float
    torsion_constant: float
    shear_correction_factor: float = 5.0 / 6.0

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
        if not (0 < self.shear_correction_factor <= 1.0):
            raise ValueError("shear_correction_factor must be in (0, 1]")

    def _build_constitutive_matrices(self) -> None:
        EA = self.youngs_modulus * self.cross_section_area
        kappa_GA = self.shear_correction_factor * self.shear_modulus * self.cross_section_area
        GJ_t = self.shear_modulus * self.torsion_constant
        D = np.diag([EA, kappa_GA, GJ_t]).astype(np.float64)
        object.__setattr__(self, '_D_assembly', D)
        object.__setattr__(self, '_D_postprocess', D.copy())
        object.__setattr__(self, '_energy_components', {
            'axial': np.diag([EA, 0.0, 0.0]),
            'transverse': np.diag([0.0, kappa_GA, 0.0]),
            'torsion': np.diag([0.0, 0.0, GJ_t]),
        })

    def assembly_form(self) -> np.ndarray:
        """Return 3×3 D-matrix for stiffness assembly."""
        return self._D_assembly

    def postprocessing_form(self) -> np.ndarray:
        """Return 3×3 D-matrix for stress recovery and post-processing."""
        return self._D_postprocess

    def compute_stress_resultants(self, strain: np.ndarray) -> np.ndarray:
        """Compute stress resultants N = D @ ε. strain shape (3,) or (3, n)."""
        return self.postprocessing_form() @ np.asarray(strain, dtype=np.float64)

    def energy_density_components(self, strain: np.ndarray) -> Dict[str, float]:
        """Strain energy density breakdown: total, axial, transverse, torsion. strain shape (3,)."""
        e = np.asarray(strain, dtype=np.float64).reshape(3)
        out = {'total': float(0.5 * e.T @ self._D_postprocess @ e)}
        for k, v in self._energy_components.items():
            out[k] = float(0.5 * e.T @ v @ e)
        return out

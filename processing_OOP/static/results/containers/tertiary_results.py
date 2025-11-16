# processing_OOP\static\results\containers\tertiary_results.py

from dataclasses import dataclass
from typing import List, Optional
import numpy as np

# ─────────────────────────────────────────────────────────────
# Tertiary-level results (highly derived/post-processed)
# These are engineering quantities computed from primary/secondary
# results, typically for visualization and design verification.
# ─────────────────────────────────────────────────────────────

@dataclass
class TertiaryResults:
    """
    Container for tertiary (highly derived) engineering results.
    
    These are computed from secondary results and represent
    design-critical quantities like section forces, principal
    stresses, failure criteria, etc.
    
    Storage follows resolution-first organization (consistent with
    primary/secondary results):
    - Gaussian resolution: `tertiary_results/gaussian/`
    - Elemental resolution: `tertiary_results/elemental/`
    
    Attributes
    ----------
    section_forces : Optional[List[List[np.ndarray]]]
        Section force resultants [N, Vy, Vz, T, My, Mz] at Gauss points
        Shape: List[element] -> List[gauss_point] -> np.ndarray(shape=(6,))
        Storage: `tertiary_results/gaussian/section_forces/`
        Native Resolution: Gaussian
    
    principal_stresses : Optional[List[List[np.ndarray]]]
        Principal stress components [σ1, σ2, σ3] at Gauss points
        Shape: List[element] -> List[gauss_point] -> np.ndarray(shape=(3,))
        Storage: `tertiary_results/gaussian/principal_stress/`
        Native Resolution: Gaussian
    
    von_mises_stress : Optional[List[List[float]]]
        Von Mises equivalent stress at Gauss points
        Shape: List[element] -> List[gauss_point] -> float
        Storage: Summary CSV only (`tertiary_results/tertiary_summary.csv`)
        Native Resolution: Gaussian
    
    max_shear_stress : Optional[List[List[float]]]
        Maximum shear stress at Gauss points
        Shape: List[element] -> List[gauss_point] -> float
        Storage: Summary CSV only (`tertiary_results/tertiary_summary.csv`)
        Native Resolution: Gaussian
    
    failure_index : Optional[List[List[float]]]
        Material failure index (stress/yield) at Gauss points
        Shape: List[element] -> List[gauss_point] -> float
        Storage: Summary CSV only (`tertiary_results/tertiary_summary.csv`)
        Native Resolution: Gaussian
    
    total_strain_energy : Optional[List[float]] = None
        Total strain energy per element (integrated from energy density)
        Shape: List[element] -> float
        Units: Joules (J)
        Storage: `tertiary_results/elemental/total_strain_energy.csv`
        Native Resolution: Elemental (integrated from Gaussian)
    
    integrated_section_forces : Optional[List[np.ndarray]] = None
        Integrated section force resultants per element [N, Vy, Vz, T, My, Mz]
        Shape: List[element] -> np.ndarray(shape=(6,))
        Units: [N, N, N, N⋅m, N⋅m, N⋅m]
        Storage: `tertiary_results/elemental/integrated_section_forces.csv`
        Native Resolution: Elemental (integrated from Gaussian)
        Note: For beam elements, these represent average or representative
        section forces over the element length
    """
    
    section_forces: Optional[List[List[np.ndarray]]] = None
    principal_stresses: Optional[List[List[np.ndarray]]] = None
    von_mises_stress: Optional[List[List[float]]] = None
    max_shear_stress: Optional[List[List[float]]] = None
    failure_index: Optional[List[List[float]]] = None
    
    # Integrated elemental results
    total_strain_energy: Optional[List[float]] = None
    integrated_section_forces: Optional[List[np.ndarray]] = None


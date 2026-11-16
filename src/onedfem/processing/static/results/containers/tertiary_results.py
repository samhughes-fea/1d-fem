# processing\static\results\containers\tertiary_results.py

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

    nodal_section_forces : Optional[np.ndarray] = None
        Section force resultants [N, Vy, Vz, T, My, Mz] projected to nodes
        Shape: (n_nodes, 6), column order [N, Vy, Vz, T, My, Mz]
        Storage: `tertiary_results/nodal/nodal_section_forces.csv`
        Native Resolution: Nodal (projected from Gaussian via shape functions)

    element_error_eta : Optional[List[float]] = None
        Zienkiewicz–Zhu energy-norm indicator per element from Gauss resultants
        versus the nodal-projected recovered field. No adaptivity.
        Storage: `tertiary_results/elemental/element_error_estimate.csv`

    element_error_eta_rel : Optional[List[float]] = None
        Relative indicator ``eta_e / sqrt(2 U_e)`` per element.

    element_error_eta_global : Optional[float] = None
        ``sqrt(sum eta_e**2)`` over the mesh.
    """
    
    section_forces: Optional[List[List[np.ndarray]]] = None
    principal_stresses: Optional[List[List[np.ndarray]]] = None
    von_mises_stress: Optional[List[List[float]]] = None
    max_shear_stress: Optional[List[List[float]]] = None
    
    # Integrated elemental results
    total_strain_energy: Optional[List[float]] = None
    integrated_section_forces: Optional[List[np.ndarray]] = None

    # Nodal resolution (projected from Gaussian)
    nodal_section_forces: Optional[np.ndarray] = None

    # Element error estimate (from Gaussian + nodal-projected resultants; no adaptivity)
    element_error_eta: Optional[List[float]] = None
    element_error_eta_rel: Optional[List[float]] = None
    element_error_eta_global: Optional[float] = None


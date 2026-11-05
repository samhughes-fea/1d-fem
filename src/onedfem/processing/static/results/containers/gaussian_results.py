# processing\static\results\containers\gaussian_results.py

from dataclasses import dataclass
from typing import List, Optional
import numpy as np

@dataclass
class GaussianResults:
    """
    Container for Gaussian (integration point) resolution results.
    
    **Native Resolution**: Gaussian - these quantities are FIRST computed at 
    integration points (Gauss points) within elements.
    
    These are pure field quantities computed directly at integration points:
    - Strain: kinematic field derived from displacements via `ε = B @ U_e`
    - Stress: constitutive response from strain via `σ = D @ ε`
    - Energy density: scalar energy field via `w = 0.5 * ε^T @ σ`
    
    These quantities can be:
    - Projected to nodal resolution using shape functions
    - Integrated to elemental resolution via quadrature
    
    Note: Section forces [N, Vy, Vz, T, My, Mz] are integrated stress
    resultants and belong in TertiaryResults, not here.
    """
    
    gauss_coords: Optional[List[np.ndarray]] = None
    # Shape: List[element] -> np.ndarray(n_gauss_points, dim)
    # Physical coordinates of integration points

    gauss_points: Optional[List[tuple]] = None
    # One (element_id, gauss_point, xi, weight) per Gauss point, in the order of the strain rows

    strain: Optional[List[List[np.ndarray]]] = None
    # Shape: List[element] -> List[gauss_point] -> np.ndarray(shape=(n_strain_components,))
    # Beam generalised strains [ε_x, κ_y, κ_z, γ_xy, γ_xz, φ_x] (bar 2, truss 3 components)

    stress: Optional[List[List[np.ndarray]]] = None
    # Shape: List[element] -> List[gauss_point] -> np.ndarray(shape=(n_stress_components,))
    # Beam stress resultants [N, M_y, M_z, V_y, V_z, T], conjugate to the strains

    internal_energy_density: Optional[List[List[np.ndarray]]] = None
    # Shape: List[element] -> List[gauss_point] -> np.ndarray(shape=(1,)) or float
    # Strain energy density u = 0.5 * σ^T * ε
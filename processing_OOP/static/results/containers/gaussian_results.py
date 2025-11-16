# processing_OOP\static\results\containers\gaussian_results.py

from dataclasses import dataclass
from typing import List, Optional
import numpy as np

@dataclass
class GaussianResults:
    """
    Container for Gaussian (integration point) resolution results.
    
    These are pure field quantities computed directly at integration points:
    - Strain: kinematic field derived from displacements
    - Stress: constitutive response from strain
    - Energy density: scalar energy field
    
    Note: Section forces [N, Vy, Vz, T, My, Mz] are integrated stress
    resultants and belong in TertiaryResults, not here.
    """
    
    gauss_coords: Optional[List[np.ndarray]] = None
    # Shape: List[element] -> np.ndarray(n_gauss_points, dim)
    # Physical coordinates of integration points

    strain: Optional[List[List[np.ndarray]]] = None
    # Shape: List[element] -> List[gauss_point] -> np.ndarray(shape=(n_strain_components,))
    # Strain tensor components [ε_xx, ε_yy, ε_zz, γ_xy, γ_yz, γ_xz]

    stress: Optional[List[List[np.ndarray]]] = None
    # Shape: List[element] -> List[gauss_point] -> np.ndarray(shape=(n_stress_components,))
    # Stress tensor components [σ_xx, σ_yy, σ_zz, τ_xy, τ_yz, τ_xz]

    internal_energy_density: Optional[List[List[np.ndarray]]] = None
    # Shape: List[element] -> List[gauss_point] -> np.ndarray(shape=(1,)) or float
    # Strain energy density u = 0.5 * σ^T * ε
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
    
    Attributes
    ----------
    section_forces : Optional[List[List[np.ndarray]]]
        Section force resultants [N, Vy, Vz, T, My, Mz] at Gauss points
        Shape: List[element] -> List[gauss_point] -> np.ndarray(shape=(6,))
    
    principal_stresses : Optional[List[List[np.ndarray]]]
        Principal stress components [σ1, σ2, σ3] at Gauss points
        Shape: List[element] -> List[gauss_point] -> np.ndarray(shape=(3,))
    
    von_mises_stress : Optional[List[List[float]]]
        Von Mises equivalent stress at Gauss points
        Shape: List[element] -> List[gauss_point] -> float
    
    max_shear_stress : Optional[List[List[float]]]
        Maximum shear stress at Gauss points
        Shape: List[element] -> List[gauss_point] -> float
    
    failure_index : Optional[List[List[float]]]
        Material failure index (stress/yield) at Gauss points
        Shape: List[element] -> List[gauss_point] -> float
    """
    
    section_forces: Optional[List[List[np.ndarray]]] = None
    principal_stresses: Optional[List[List[np.ndarray]]] = None
    von_mises_stress: Optional[List[List[float]]] = None
    max_shear_stress: Optional[List[List[float]]] = None
    failure_index: Optional[List[List[float]]] = None


# pre_processing/element_library/gauss_point_data.py

"""
Gauss Point Formulation Cache Data Structures

This module defines data containers for caching element formulation data
at the Gauss integration point level. This enables post-processing of
stresses, strains, and section forces without recomputing shape functions.

Key Features:
- Stores intermediate matrices (B, D, J) at each Gauss point
- Tracks element formulation type for multi-element meshes
- Supports all classical beam theories (Euler-Bernoulli, Timoshenko, Levinson)
- Enables efficient stress recovery and engineering analysis

Classes:
- GaussPointData: Single Gauss point stiffness formulation data
- ElementObject: Complete element stiffness with all GP data
- ForceGaussPointData: Single Gauss point force formulation data
- ForceObject: Complete element force with all GP data
"""

from dataclasses import dataclass, field
from typing import List, Optional
import numpy as np


@dataclass
class GaussPointData:
    """
    Stores formulation data at a single Gauss integration point.
    
    This includes all matrices and vectors needed to compute element
    stiffness and force contributions at this specific Gauss point.
    
    Attributes
    ----------
    xi : float
        Natural coordinate of the Gauss point (-1 to 1)
    weight : float
        Integration weight for this Gauss point
    B_matrix : np.ndarray
        Strain-displacement matrix at this Gauss point (n_strain x n_dof)
    D_matrix : np.ndarray
        Constitutive (material stiffness) matrix (n_strain x n_strain)
    jacobian : float
        Jacobian determinant for coordinate transformation
    shape_functions : Optional[np.ndarray]
        Shape function values at this point
    shape_derivatives : Optional[np.ndarray]
        Shape function derivatives at this point
    """
    xi: float
    weight: float
    B_matrix: np.ndarray
    D_matrix: np.ndarray
    jacobian: float
    shape_functions: Optional[np.ndarray] = None
    shape_derivatives: Optional[np.ndarray] = None


@dataclass
class ElementObject:
    """
    Caches element stiffness formulation data at the Gauss point level.
    
    Instead of just storing the final K_e matrix, this object preserves
    the intermediate formulation data (B, D, J) at each Gauss point,
    enabling post-processing of stresses, strains, and section forces
    without recomputing shape functions.
    
    Attributes
    ----------
    element_id : int
        Unique identifier for this element
    element_type : str
        Beam formulation type (e.g., "Euler-Bernoulli-3D", "Timoshenko-3D", "Levinson-3D")
    K_e : np.ndarray
        Assembled element stiffness matrix (n_dof x n_dof)
    gauss_data : List[GaussPointData]
        Formulation data at each Gauss integration point
    integration_scheme : str
        Description of integration scheme used (e.g., "Gauss-Legendre")
    """
    element_id: int
    element_type: str
    K_e: np.ndarray
    gauss_data: List[GaussPointData] = field(default_factory=list)
    integration_scheme: str = "Gauss-Legendre"
    
    def get_gauss_point(self, index: int) -> GaussPointData:
        """Retrieve data for a specific Gauss point by index."""
        return self.gauss_data[index]
    
    @property
    def n_gauss_points(self) -> int:
        """Number of integration points."""
        return len(self.gauss_data)


@dataclass
class ForceGaussPointData:
    """
    Stores force formulation data at a single Gauss integration point.
    
    Attributes
    ----------
    xi : float
        Natural coordinate of the Gauss point
    weight : float
        Integration weight
    shape_functions : np.ndarray
        Shape function values N(xi)
    jacobian : float
        Jacobian determinant
    distributed_load : Optional[np.ndarray]
        Distributed load vector at this point
    """
    xi: float
    weight: float
    shape_functions: np.ndarray
    jacobian: float
    distributed_load: Optional[np.ndarray] = None


@dataclass
class ForceObject:
    """
    Caches element force vector formulation data at the Gauss point level.
    
    Preserves the load distribution and shape function values at each
    Gauss point, enabling detailed force/load analysis.
    
    Attributes
    ----------
    element_id : int
        Unique identifier for this element
    element_type : str
        Beam formulation type (e.g., "Euler-Bernoulli-3D", "Timoshenko-3D", "Levinson-3D")
    F_e : np.ndarray
        Assembled element force vector (n_dof,)
    gauss_data : List[ForceGaussPointData]
        Force formulation data at each Gauss point
    point_loads : Optional[np.ndarray]
        Concentrated loads applied to element nodes
    """
    element_id: int
    element_type: str
    F_e: np.ndarray
    gauss_data: List[ForceGaussPointData] = field(default_factory=list)
    point_loads: Optional[np.ndarray] = None
    
    def get_gauss_point(self, index: int) -> ForceGaussPointData:
        """Retrieve force data for a specific Gauss point by index."""
        return self.gauss_data[index]
    
    @property
    def n_gauss_points(self) -> int:
        """Number of integration points."""
        return len(self.gauss_data)


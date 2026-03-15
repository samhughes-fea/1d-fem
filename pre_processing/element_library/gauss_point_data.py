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
- StiffnessGaussPointData: Single Gauss point stiffness formulation data
- ElementObject: Complete element stiffness with all GP data
- ForceGaussPointData: Single Gauss point force formulation data
- ForceObject: Complete element force with all GP data
"""

from dataclasses import dataclass, field
from typing import Callable, List, Optional, Tuple
import numpy as np


@dataclass
class StiffnessGaussPointData:
    """
    Stores stiffness formulation data at a single Gauss integration point.
    
    This includes all matrices and vectors needed to compute element
    stiffness at this specific Gauss point. When used in the formulation
    cache for the results pipeline, ``shape_functions`` and
    ``shape_derivatives`` must be populated at every Gauss point so that
    result computers and projection use the same formulation that built K_e.
    
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
        Shape function values N(xi) at this point. Required when used in
        the results formulation cache.
    shape_derivatives : Optional[np.ndarray]
        Shape function derivatives dN/dxi at this point. Required when
        used in the results formulation cache.
    strain_at_convergence : Optional[np.ndarray]
        When set (e.g. after nonlinear solve), strain E_lin + E_nl at this
        Gauss point. Secondary results use this instead of B_matrix @ U_e
        when present. Default None.
    """
    xi: float
    weight: float
    B_matrix: np.ndarray
    D_matrix: np.ndarray
    jacobian: float
    shape_functions: Optional[np.ndarray] = None
    shape_derivatives: Optional[np.ndarray] = None
    strain_at_convergence: Optional[np.ndarray] = None


@dataclass
class ElementObject:
    """
    Caches element stiffness formulation data at the Gauss point level.
    
    Instead of just storing the final K_e matrix, this object preserves
    the intermediate formulation data (B, D, J) at each Gauss point,
    enabling post-processing of stresses, strains, and section forces
    without recomputing shape functions. When used in the formulation
    cache for the results pipeline, every entry in ``gauss_data`` must
    have ``shape_functions`` and ``shape_derivatives`` populated.
    
    Attributes
    ----------
    element_id : int
        Unique identifier for this element
    element_type : str
        Beam formulation type (e.g., "Euler-Bernoulli-3D", "Timoshenko-3D", "Levinson-3D")
    K_e : np.ndarray
        Assembled element stiffness matrix (n_dof x n_dof)
    gauss_data : List[StiffnessGaussPointData]
        Formulation data at each Gauss integration point
    integration_scheme : str
        Description of integration scheme used (e.g., "Gauss-Legendre")
    evaluate_shape_functions : Optional[Callable]
        Optional callable for evaluating N(ξ) at arbitrary ξ. When set, signature is
        (xi: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray] returning
        (N, dN_dξ, d2N_dξ2) in natural coordinates. Used for post-processing between
        Gauss points. Not serialized; in-memory only.
    shape_function_N_coefficients : Optional[np.ndarray]
        (B2) Polynomial coefficients for N(ξ) in monomial basis. Shape (12, 6, 4):
        N_coefficients[dof, comp, k] = coefficient of ξ^k. Enables evaluation after
        save/load without the element class. Serializable.
    shape_function_dN_dxi_coefficients : Optional[np.ndarray]
        (B2) Coefficients for dN/dξ, shape (12, 6, 4). Same convention.
    shape_function_d2N_dxi2_coefficients : Optional[np.ndarray]
        (B2) Coefficients for d²N/dξ², shape (12, 6, 4). Same convention.
    """
    element_id: int
    element_type: str
    K_e: np.ndarray
    gauss_data: List[StiffnessGaussPointData] = field(default_factory=list)
    integration_scheme: str = "Gauss-Legendre"
    evaluate_shape_functions: Optional[
        Callable[[np.ndarray], Tuple[np.ndarray, np.ndarray, np.ndarray]]
    ] = None
    shape_function_N_coefficients: Optional[np.ndarray] = None
    shape_function_dN_dxi_coefficients: Optional[np.ndarray] = None
    shape_function_d2N_dxi2_coefficients: Optional[np.ndarray] = None
    
    def get_gauss_point(self, index: int) -> StiffnessGaussPointData:
        """Retrieve data for a specific Gauss point by index."""
        return self.gauss_data[index]
    
    @property
    def n_gauss_points(self) -> int:
        """Number of integration points."""
        return len(self.gauss_data)


# Backward-compatibility alias (prefer StiffnessGaussPointData for new code).
GaussPointData = StiffnessGaussPointData


@dataclass
class ForceGaussPointData:
    """
    Stores force formulation data at a single Gauss integration point.
    
    When used in the formulation cache for the results pipeline,
    ``shape_functions`` must be set at every Gauss point (it is a
    required field for force assembly and for result computers).
    
    Attributes
    ----------
    xi : float
        Natural coordinate of the Gauss point
    weight : float
        Integration weight
    shape_functions : np.ndarray
        Shape function values N(xi). Required for the results pipeline.
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
    Gauss point, enabling detailed force/load analysis. When used in the
    formulation cache for the results pipeline, every entry in
    ``gauss_data`` must have ``shape_functions`` set.
    
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


@dataclass
class MassObject:
    """
    Caches element mass matrix for modal/dynamic assembly.

    Attributes
    ----------
    element_id : int
        Unique identifier for this element
    element_type : str
        Element type name (e.g., "Bar-3D")
    M_e : np.ndarray
        Assembled element mass matrix (n_dof x n_dof), same shape and DOF ordering as K_e.
    """
    element_id: int
    element_type: str
    M_e: np.ndarray


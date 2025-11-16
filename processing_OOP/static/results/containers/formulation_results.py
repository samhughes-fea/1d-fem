"""
FormulationResultSet Container

Caches element formulation data with Gauss-level resolution for efficient
post-processing stress/strain recovery without recomputing shape functions.

This container stores:
- ElementObject: Cached K_e and Gauss point data (B, D, J, N) for each element
- ForceObject: Cached F_e and Gauss point data (N, J, distributed loads) for each element

Benefits:
- Performance: No shape function recomputation during results computation
- Traceability: Element type tracking for mixed-element meshes
- Accuracy: Use exact same B, D matrices as used in stiffness assembly
"""

from dataclasses import dataclass, field
from typing import List
from pre_processing.element_library.gauss_point_data import ElementObject, ForceObject


@dataclass
class FormulationResultSet:
    """
    Caches element formulation data with Gauss-level resolution.
    
    Attributes:
        element_objects: List of ElementObject containing K_e and Gauss data
        force_objects: List of ForceObject containing F_e and force Gauss data
    """
    element_objects: List[ElementObject] = field(default_factory=list)
    force_objects: List[ForceObject] = field(default_factory=list)


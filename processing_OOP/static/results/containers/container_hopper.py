# processing_OOP\static\results\containers\container_hopper.py

from dataclasses import dataclass
from typing import Optional, List
from .global_results import GlobalResults
from .elemental_results import ElementalResults
from .nodal_results import NodalResults
from .gaussian_results import GaussianResults
from .tertiary_results import TertiaryResults
from .map_results import MapEntry

# ─────────────────────────────────────────────────────────────
# Primary, secondary, and tertiary results stratified by resolution
# ─────────────────────────────────────────────────────────────

@dataclass
class PrimaryResultSet:
    """
    Holds all primary (first-order) simulation results:
    - Global system-level outputs
    - Element-wise displacements and reactions
    - Optional nodal and Gaussian-level results (if projected)
    """
    global_results: GlobalResults
    elemental_results: ElementalResults
    nodal_results: Optional[NodalResults] = None
    gaussian_results: Optional[GaussianResults] = None

@dataclass
class SecondaryResultSet:
    """
    Holds all secondary (derived) results, typically post-processed from primary data.
    Fields are optional to allow for partial or staged computation.
    """
    global_results: Optional[GlobalResults] = None
    elemental_results: Optional[ElementalResults] = None
    nodal_results: Optional[NodalResults] = None
    gaussian_results: Optional[GaussianResults] = None

@dataclass
class TertiaryResultSet:
    """
    Holds all tertiary (highly derived) engineering results.
    
    These are design-critical quantities computed from secondary results:
    - Section force resultants
    - Principal stresses
    - Von Mises stress
    - Failure indices
    """
    tertiary_results: TertiaryResults

# ─────────────────────────────────────────────────────────────
# All intermediate DOF space transformation maps
# ─────────────────────────────────────────────────────────────

@dataclass
class IndexMapSet:
    assembly_map: Optional[List[MapEntry]] = None
    modification_map: Optional[List[MapEntry]] = None
    condensation_map: Optional[List[MapEntry]] = None
    reconstruction_map: Optional[List[MapEntry]] = None

# processing\static\results\containers\container_hopper.py

from dataclasses import dataclass
from typing import Optional, List
from .global_results import GlobalResults
from .elemental_results import ElementalResults
from .nodal_results import NodalResults
from .gaussian_results import GaussianResults
from .tertiary_results import TertiaryResults
from .map_results import MapEntry

"""
Resolution Levels
----------------
The results system uses four resolution levels:

Global: System-wide quantities (primary results)
    - Quantities computed at the assembled system level
    - Examples: U_global, R_global, K_global, F_global
    - Visualization: Typically shown as summary values or system-wide plots

Gaussian: Native integration point quantities (where strain, stress, etc. are first computed)
    - Quantities computed at integration points within elements
    - Examples: ε (strain), σ (stress), w (energy density)
    - Visualization: Discrete markers at Gauss point physical locations (within elements, not at nodes)
    - Note: Gauss points are typically located inside elements (e.g., at natural coordinates like ξ = ±1/√3)

Nodal: Projected quantities (all nodal results are projections from Gaussian)
    - Quantities interpolated/extrapolated to nodes using shape functions
    - Examples: nodal strain, nodal stress, nodal energy density
    - Visualization: Markers at node locations, with continuous field interpolated using shape functions
    - The interpolated line/field follows the shape function path between nodes

Elemental: Per-element quantities (can be native, disassembled, or integrated)
    - Native: Formulated at element level (K_e, F_e)
    - Disassembled: Extracted from global (U_e, R_e)
    - Integrated: Computed via quadrature from Gaussian (total_strain_energy)
    - Visualization: Typically shown as bar charts, element-wise color maps, or summary values per element

Visualization Approach
----------------------
When plotting quantities across resolutions:
1. Gaussian: Plot discrete markers at Gauss point physical locations
2. Nodal: Plot markers at node locations, then interpolate using shape functions to create smooth continuous field
3. The interpolated visualization uses the same shape functions that were used for projection, creating a consistent representation
"""

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

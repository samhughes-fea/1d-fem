# processing_OOP\static\results\containers\nodal_results.py

from dataclasses import dataclass
from typing import Optional
import numpy as np

# ─────────────────────────────────────────────────────────────
# Nodal-level results (interpolated field quantities at nodes)
# These results are typically obtained by extrapolating from
# integration points or directly solving for nodal quantities.
# ─────────────────────────────────────────────────────────────

@dataclass
class NodalResults:
    """
    Container for nodal resolution results.
    
    **Native Resolution**: None - all nodal results are PROJECTIONS from Gaussian.
    
    Contains pure field quantities interpolated/extrapolated to nodes using
    cached shape functions:
    - Strain: projected from Gaussian resolution
    - Stress: projected from Gaussian resolution
    - Energy density: projected from Gaussian resolution
    
    These quantities are derived from GaussianResults via shape function
    interpolation/extrapolation, not computed directly at nodes.
    
    Note: Section forces [N, Vy, Vz, T, My, Mz] are integrated stress
    resultants and belong in TertiaryResults, not as nodal field values.
    """
    
    strain: Optional[np.ndarray] = None
    # Shape: (n_nodes, 6)
    # Components: [ε_xx, ε_yy, ε_zz, γ_xy, γ_yz, γ_xz]

    stress: Optional[np.ndarray] = None
    # Shape: (n_nodes, 6)
    # Components: [σ_xx, σ_yy, σ_zz, τ_xy, τ_yz, τ_xz]

    strain_energy_density: Optional[np.ndarray] = None
    # Shape: (n_nodes,)
    # Scalar strain energy density per node (typically in J/m³ or similar)

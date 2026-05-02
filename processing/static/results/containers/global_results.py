# processing\static\results\containers\global_results.py

from dataclasses import dataclass
from typing import Optional
import numpy as np
import scipy.sparse as sp

# ─────────────────────────────────────────────────────────────
# Global-level results (system-wide matrices and vectors)
# 
# Native Resolution: Global - quantities are FIRST computed at the
# assembled system level (smallest resolution for system-wide quantities).
# ─────────────────────────────────────────────────────────────

@dataclass
class GlobalResults:
    # --------------------------------AssembleGlobalSystem outputs
    F_global: Optional[np.ndarray] = None
    K_global: Optional[sp.csr_matrix] = None
    # --------------------------------ModifyGlobalSystem outputs
    F_mod: Optional[np.ndarray] = None
    K_mod: Optional[sp.csr_matrix] = None
    # --------------------------------CondenseModifiedSystem outputs
    F_cond: Optional[np.ndarray] = None
    K_cond: Optional[sp.csr_matrix] = None
    # --------------------------------SolveCondensedSystem outputs
    U_cond: Optional[np.ndarray] = None
    # --------------------------------ReconstructGlobalSystem outputs
    U_global: Optional[np.ndarray] = None
    # --------------------------------PrimaryResultsOrchestrator outputs
    R_global: Optional[np.ndarray] = None
    R_residual: Optional[np.ndarray] = None
    # --------------------------------Nonlinear static only (optional summary columns)
    newton_converged: Optional[bool] = None
    newton_iterations_total: Optional[int] = None
    load_increments_completed: Optional[int] = None
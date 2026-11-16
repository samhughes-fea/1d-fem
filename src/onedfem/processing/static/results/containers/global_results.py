# processing\static\results\containers\global_results.py

from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Union
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


@dataclass
class ModelTotals:
    """Whole-model scalars of a linear static solve.

    ``strain_energy`` is the sum of the per-element integrated strain energies.
    ``external_work`` is ``0.5 U^T F`` of the assembled load. ``equilibrium_residual``
    is the unscaled free-DOF residual ``||(F - K U)[free]||``.
    """

    strain_energy: float
    external_work: float
    equilibrium_residual: float

    def to_csv(self, path: Union[str, Path]) -> Path:
        """Write ``quantity,value`` rows. Returns the path written."""
        out = Path(path)
        out.parent.mkdir(parents=True, exist_ok=True)
        with out.open("w", encoding="utf-8") as handle:
            handle.write("quantity,value\n")
            handle.write(f"strain_energy,{self.strain_energy:.16e}\n")
            handle.write(f"external_work,{self.external_work:.16e}\n")
            handle.write(f"equilibrium_residual,{self.equilibrium_residual:.16e}\n")
        return out
# processing/harmonic/operations/modify_harmonic_structural.py
"""BCs on K, M and zero constrained entries on F for harmonic."""

from __future__ import annotations

from pathlib import Path
from typing import Optional, Sequence, Tuple, Union

import numpy as np
from scipy.sparse import csr_matrix

from processing.eigen.boundary_conditions import apply_boundary_conditions
from processing.spectral.spectral_diagnostics import log_spectral_diagnostics

from processing.harmonic.operations._logging import init_stage_logger


class ModifyHarmonicStructuralMatrices:
    def __init__(
        self,
        prescribed_displacements: Optional[dict] = None,
        job_results_dir: Optional[Union[str, Path]] = None,
        fixed_dofs: Optional[Sequence[int]] = None,
    ):
        self.prescribed_displacements = prescribed_displacements
        self.fixed_dofs = fixed_dofs
        self.job_results_dir = Path(job_results_dir) if job_results_dir else None
        self._log = init_stage_logger(
            "ModifyHarmonicStructuralMatrices", self.job_results_dir
        )

    def run(
        self,
        K_global: csr_matrix,
        M_global: csr_matrix,
        F_global: np.ndarray,
    ) -> Tuple[csr_matrix, csr_matrix, np.ndarray, np.ndarray]:
        self._log.info("Applying harmonic boundary conditions")
        K_mod, M_mod, bc_dofs = apply_boundary_conditions(
            K_global,
            M_global,
            fixed_dofs=self.fixed_dofs,
            prescribed_displacements=self.prescribed_displacements,
        )
        log_spectral_diagnostics(
            K_mod,
            M_mod,
            str(self.job_results_dir) if self.job_results_dir else None,
        )
        F_mod = np.asarray(F_global, dtype=np.float64).copy()
        if bc_dofs.size:
            bci = np.asarray(bc_dofs, dtype=np.int64)
            F_mod[bci] = 0.0
        self._log.info("Harmonic BCs: %s constrained DOFs", bc_dofs.size)
        return K_mod, M_mod, bc_dofs, F_mod

# processing/spectral/operations/modify_spectral_global.py
"""Apply BCs to global K and M for undamped vibration (penalty on K, M diagonal)."""

from __future__ import annotations

from pathlib import Path
from typing import Optional, Sequence, Tuple, Union

import numpy as np
from scipy.sparse import csr_matrix

from processing.eigen.boundary_conditions import apply_boundary_conditions
from processing.spectral.spectral_diagnostics import (
    log_spectral_constrained_dofs,
    log_spectral_diagnostics,
)

from processing.common.stage_logging import init_stage_logger


class ModifySpectralGlobalSystem:
    """Wrap :func:`processing.eigen.boundary_conditions.apply_boundary_conditions`."""

    def __init__(
        self,
        job_results_dir: Optional[Union[str, Path]] = None,
        fixed_dofs: Optional[Sequence[int]] = None,
        prescribed_displacements: Optional[dict] = None,
    ):
        self.job_results_dir = Path(job_results_dir) if job_results_dir else None
        self.fixed_dofs = fixed_dofs
        self.prescribed_displacements = prescribed_displacements
        self._log = init_stage_logger("ModifySpectralGlobalSystem", self.job_results_dir)

    def run(
        self,
        K_global: csr_matrix,
        M_global: csr_matrix,
    ) -> Tuple[csr_matrix, csr_matrix, np.ndarray]:
        self._log.info("Applying spectral boundary conditions")
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
        log_spectral_constrained_dofs(
            bc_dofs,
            total_dof=int(K_mod.shape[0]),
            job_results_dir=str(self.job_results_dir) if self.job_results_dir else None,
        )
        self._log.info("Spectral BCs applied: %s constrained DOFs", np.asarray(bc_dofs).size)
        return K_mod, M_mod, bc_dofs

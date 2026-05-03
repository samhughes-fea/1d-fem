# processing/dynamic/operations/modify_dynamic_global.py
"""BCs on K, M, C for transient dynamics (penalty / diagonal treatment)."""

from __future__ import annotations

from pathlib import Path
from typing import Optional, Sequence, Tuple, Union

import numpy as np
from scipy.sparse import csr_matrix

from processing.dynamic.boundary_conditions import apply_boundary_conditions

from processing.dynamic.operations._logging import init_stage_logger


class ModifyDynamicGlobalSystem:
    def __init__(
        self,
        fixed_dofs: Optional[Sequence[int]] = None,
        prescribed_displacements: Optional[dict] = None,
        job_results_dir: Optional[Union[str, Path]] = None,
    ):
        self.fixed_dofs = fixed_dofs
        self.prescribed_displacements = prescribed_displacements
        self.job_results_dir = Path(job_results_dir) if job_results_dir else None
        self._log = init_stage_logger("ModifyDynamicGlobalSystem", self.job_results_dir)

    def run(
        self,
        K_global: csr_matrix,
        M_global: csr_matrix,
        C_global: Optional[csr_matrix],
    ) -> Tuple[csr_matrix, csr_matrix, Optional[csr_matrix], np.ndarray]:
        self._log.info("Applying dynamic boundary conditions")
        K_mod, M_mod, C_mod, bc_dofs = apply_boundary_conditions(
            K_global,
            M_global,
            C_global,
            fixed_dofs=self.fixed_dofs,
            prescribed_displacements=self.prescribed_displacements,
        )
        self._log.info("Dynamic BCs: %s constrained DOFs", np.asarray(bc_dofs).size)
        return K_mod, M_mod, C_mod, bc_dofs

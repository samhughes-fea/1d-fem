# processing/dynamic/operations/assemble_dynamic_global.py
"""Global K, M, C assembly for transient dynamics."""

from __future__ import annotations

from pathlib import Path
from typing import List, Optional, Union

from scipy.sparse import csr_matrix

from processing.dynamic.assembly import assemble_global_system

from processing.dynamic.operations._logging import init_stage_logger


class AssembleDynamicGlobalSystem:
    def __init__(
        self,
        elements: List,
        element_stiffness_matrices,
        element_mass_matrices,
        total_dof: int,
        job_results_dir: Optional[Union[str, Path]] = None,
        element_damping_matrices=None,
    ):
        self.elements = list(elements)
        self.element_stiffness_matrices = element_stiffness_matrices
        self.element_mass_matrices = element_mass_matrices
        self.element_damping_matrices = element_damping_matrices
        self.total_dof = int(total_dof)
        self.job_results_dir = Path(job_results_dir) if job_results_dir else None
        self._log = init_stage_logger("AssembleDynamicGlobalSystem", self.job_results_dir)

    def run(self) -> tuple[csr_matrix, csr_matrix, Optional[csr_matrix], list]:
        self._log.info("Dynamic global assembly: n_el=%s dof=%s", len(self.elements), self.total_dof)
        K_global, M_global, C_global, local_map = assemble_global_system(
            elements=self.elements,
            element_stiffness_matrices=self.element_stiffness_matrices,
            element_mass_matrices=self.element_mass_matrices,
            element_damping_matrices=self.element_damping_matrices,
            total_dof=self.total_dof,
            job_results_dir=str(self.job_results_dir) if self.job_results_dir else None,
        )
        self._log.info(
            "K nnz=%s M nnz=%s C=%s",
            K_global.nnz,
            M_global.nnz,
            "None" if C_global is None else C_global.nnz,
        )
        return K_global, M_global, C_global, local_map

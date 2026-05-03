# processing/spectral/operations/assemble_spectral_global.py
"""Global K / M assembly for §2 spectral problems."""

from __future__ import annotations

from pathlib import Path
from typing import List, Optional, Union

from scipy.sparse import csr_matrix

from processing.eigen.assembly import assemble_global_matrices
from processing.spectral.spectral_diagnostics import log_spectral_diagnostics

from processing.common.stage_logging import init_stage_logger


class AssembleSpectralGlobalSystem:
    """
    Scatter element matrices to global ``K`` and ``M`` (delegates to
    :func:`processing.eigen.assembly.assemble_global_matrices`), then log diagnostics.
    """

    def __init__(
        self,
        elements: List,
        element_stiffness_matrices: List,
        element_mass_matrices: List,
        total_dof: int,
        job_results_dir: Optional[Union[str, Path]] = None,
    ):
        self.elements = list(elements)
        self.element_stiffness_matrices = element_stiffness_matrices
        self.element_mass_matrices = element_mass_matrices
        self.total_dof = int(total_dof)
        self.job_results_dir = Path(job_results_dir) if job_results_dir else None
        self._log = init_stage_logger("AssembleSpectralGlobalSystem", self.job_results_dir)

    def run(
        self,
    ) -> tuple[csr_matrix, csr_matrix, list]:
        self._log.info(
            "Assembling spectral global K/M: n_el=%s total_dof=%s",
            len(self.elements),
            self.total_dof,
        )
        K_global, M_global, local_map = assemble_global_matrices(
            elements=self.elements,
            element_stiffness_matrices=self.element_stiffness_matrices,
            element_mass_matrices=self.element_mass_matrices,
            total_dof=self.total_dof,
            job_results_dir=str(self.job_results_dir) if self.job_results_dir else None,
        )
        if K_global is None or M_global is None:
            raise ValueError("Spectral assembly produced None K or M")
        log_spectral_diagnostics(
            K_global,
            M_global,
            str(self.job_results_dir) if self.job_results_dir else None,
        )
        self._log.info("Spectral assembly complete: K nnz=%s M nnz=%s", K_global.nnz, M_global.nnz)
        return K_global, M_global, local_map

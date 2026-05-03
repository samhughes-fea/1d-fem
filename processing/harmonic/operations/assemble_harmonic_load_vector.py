# processing/harmonic/operations/assemble_harmonic_load_vector.py
"""Global force vector for harmonic loads (static-style assembly)."""

from __future__ import annotations

from pathlib import Path
from typing import List, Optional, Union

import numpy as np

from processing.static.operations.assembly import AssembleGlobalSystem

from processing.harmonic.operations._logging import init_stage_logger


class AssembleHarmonicLoadVector:
    def __init__(
        self,
        elements: List,
        element_stiffness_matrices: List,
        element_force_vectors: List,
        total_dof: int,
        job_results_dir: Optional[Union[str, Path]] = None,
    ):
        self.elements = list(elements)
        self.element_stiffness_matrices = list(element_stiffness_matrices)
        self.element_force_vectors = list(element_force_vectors)
        self.total_dof = int(total_dof)
        self.job_results_dir = Path(job_results_dir) if job_results_dir else None
        self._log = init_stage_logger("AssembleHarmonicLoadVector", self.job_results_dir)

    def run(self) -> np.ndarray:
        self._log.info("Assembling harmonic global load vector")
        assembler = AssembleGlobalSystem(
            elements=self.elements,
            element_stiffness_matrices=self.element_stiffness_matrices,
            element_force_vectors=self.element_force_vectors,
            total_dof=self.total_dof,
            job_results_dir=str(self.job_results_dir) if self.job_results_dir else None,
        )
        _, F_global, _, _ = assembler.assemble()
        F = np.asarray(F_global).ravel()
        self._log.info("F_global assembled: norm=%.3e", float(np.linalg.norm(F)))
        return F

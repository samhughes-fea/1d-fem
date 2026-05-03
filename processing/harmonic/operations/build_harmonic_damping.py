# processing/harmonic/operations/build_harmonic_damping.py
"""Rayleigh / modal-proportional damping matrix for harmonic sweep."""

from __future__ import annotations

from pathlib import Path
from typing import Optional, Union

from scipy.sparse import csr_matrix

from processing.harmonic.frequency_response import harmonic_damping_matrix

from processing.harmonic.operations._logging import init_stage_logger


class BuildHarmonicDampingMatrix:
    def __init__(self, job_results_dir: Optional[Union[str, Path]] = None):
        self.job_results_dir = Path(job_results_dir) if job_results_dir else None
        self._log = init_stage_logger("BuildHarmonicDampingMatrix", self.job_results_dir)

    def run(
        self,
        M_mod: csr_matrix,
        K_mod: csr_matrix,
        zeta: float,
        omega_ref: float,
        rayleigh_alpha: float,
        rayleigh_beta: float,
    ) -> csr_matrix:
        self._log.info(
            "Building harmonic damping: zeta=%s omega_ref=%s ra=%s rb=%s",
            zeta,
            omega_ref,
            rayleigh_alpha,
            rayleigh_beta,
        )
        C_mod = harmonic_damping_matrix(
            M_mod, K_mod, zeta, omega_ref, rayleigh_alpha, rayleigh_beta
        )
        self._log.info("C_mod assembled: nnz=%s", C_mod.nnz)
        return C_mod

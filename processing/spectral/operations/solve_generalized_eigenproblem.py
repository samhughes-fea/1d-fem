# processing/spectral/operations/solve_generalized_eigenproblem.py
"""Smallest generalized eigenpairs for K x = λ M x (§2 vibration)."""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Optional, Tuple, Union

import numpy as np
from scipy.sparse import csr_matrix

from processing.eigen.smallest_generalized_eigenpairs import smallest_generalized_eigenpairs

from processing.common.stage_logging import init_stage_logger

logger = logging.getLogger(__name__)


class SolveGeneralizedEigenproblem:
    """
    Solve :math:`K x = \\lambda M x` for the smallest-magnitude eigenvalues
    (delegates to :func:`processing.eigen.smallest_generalized_eigenpairs`).
    """

    def __init__(
        self,
        num_modes: int,
        context: str = "eigen vibration",
        dense_threshold: int = 512,
        job_results_dir: Optional[Union[str, Path]] = None,
    ):
        self.num_modes = int(num_modes)
        self.context = str(context)
        self.dense_threshold = int(dense_threshold)
        self.job_results_dir = Path(job_results_dir) if job_results_dir else None
        self._log = init_stage_logger("SolveGeneralizedEigenproblem", self.job_results_dir)

    def run(
        self,
        K_mod: csr_matrix,
        M_mod: csr_matrix,
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Returns
        -------
        eigenvalues, eigenvectors, frequencies_hz
            Natural frequencies in Hz are ``sqrt(|λ|)/(2π)`` in the historical convention.
        """
        self._log.info("Solving generalized eigenproblem: num_modes=%s", self.num_modes)
        eigenvalues, eigenvectors = smallest_generalized_eigenpairs(
            K_mod,
            M_mod,
            self.num_modes,
            dense_threshold=self.dense_threshold,
            context=self.context,
        )
        frequencies_hz = np.sqrt(np.abs(eigenvalues)) / (2 * np.pi)
        self._log.info("Eigen solve finished: %s modes", self.num_modes)
        return eigenvalues, eigenvectors, frequencies_hz

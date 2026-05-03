# processing/spectral/operations/prepare_spectral_local.py
"""Normalize element K_e / M_e for spectral global assembly (COO float64)."""

from __future__ import annotations

from pathlib import Path
from typing import List, Optional, Sequence, Union

import numpy as np
from scipy.sparse import coo_matrix, issparse

from processing.common.stage_logging import init_stage_logger


class PrepareSpectralLocalMatrices:
    """
    Validate and convert element stiffness / mass matrices to COO float64 lists,
    mirroring :class:`processing.static.operations.preparation.PrepareLocalSystem`
    for the stiffness path (mass has no force analogue).
    """

    def __init__(
        self,
        element_stiffness_matrices: Optional[Sequence] = None,
        element_mass_matrices: Optional[Sequence] = None,
        job_results_dir: Optional[Union[str, Path]] = None,
    ):
        self.job_results_dir = Path(job_results_dir) if job_results_dir else None
        self._log = init_stage_logger(
            "PrepareSpectralLocalMatrices", self.job_results_dir
        )
        self._Ke_raw = element_stiffness_matrices
        self._Me_raw = element_mass_matrices

    @staticmethod
    def _to_coo_list(matrices: Optional[Sequence]) -> List[coo_matrix]:
        if matrices is None:
            return []
        seq = list(matrices) if not isinstance(matrices, np.ndarray) else list(matrices.ravel())
        out: List[coo_matrix] = []
        for mat in seq:
            if mat is None:
                raise ValueError("Spectral prepare: encountered None element matrix")
            if issparse(mat):
                out.append(mat.astype(np.float64).tocoo())
            else:
                dense = np.asarray(mat, dtype=np.float64)
                out.append(coo_matrix(dense))
        return out

    def run(self) -> tuple[List[coo_matrix], List[coo_matrix]]:
        """Return ``(Ke_coo_list, Me_coo_list)``."""
        ke = self._to_coo_list(self._Ke_raw)
        me = self._to_coo_list(self._Me_raw)
        self._log.info("Prepared %s stiffness and %s mass element matrices", len(ke), len(me))
        return ke, me

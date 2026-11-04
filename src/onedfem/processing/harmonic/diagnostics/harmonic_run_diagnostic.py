"""Lightweight assembled-system logging for harmonic jobs (Section 4)."""

from __future__ import annotations

import logging
import os
from pathlib import Path
from typing import Optional, Union

import numpy as np
from scipy.sparse import csr_matrix

_logger = logging.getLogger(__name__)


def log_harmonic_structural_summary(
    K_mod: csr_matrix,
    M_mod: csr_matrix,
    C_mod: csr_matrix,
    F_mod: np.ndarray,
    *,
    n_bc_dofs: int,
    job_results_dir: Optional[Union[str, Path]] = None,
) -> None:
    """Log norms of modified K, M, C and load vector after BCs."""

    def _frob(a: csr_matrix) -> float:
        if a is None or a.nnz == 0:
            return 0.0
        return float(np.sqrt(np.sum(np.asarray(a.data, dtype=np.float64) ** 2)))

    fv = np.asarray(F_mod, dtype=np.complex128).ravel()
    fnorm = float(np.sqrt(np.real(np.vdot(fv, fv))))
    _logger.info(
        "Harmonic structural (modified): frob(K)=%.3e frob(M)=%.3e frob(C)=%.3e ||F||=%.3e n_bc=%s",
        _frob(K_mod),
        _frob(M_mod),
        _frob(C_mod),
        fnorm,
        n_bc_dofs,
    )
    if job_results_dir and os.path.isdir(str(job_results_dir)):
        root = Path(job_results_dir).resolve().parent
        log_path = root / "logs" / "harmonic_run_diagnostic.log"
        try:
            log_path.parent.mkdir(parents=True, exist_ok=True)
            with open(log_path, "a", encoding="utf-8") as fh:
                fh.write(
                    f"frob_K={_frob(K_mod):.6e} frob_M={_frob(M_mod):.6e} frob_C={_frob(C_mod):.6e} "
                    f"fnorm={fnorm:.6e} n_bc={n_bc_dofs}\n"
                )
        except OSError:
            pass

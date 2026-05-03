"""Lightweight assembled-system logging for transient jobs (Section 3)."""

from __future__ import annotations

import logging
import os
from pathlib import Path
from typing import Optional, Union

import numpy as np
from scipy.sparse import csr_matrix

_logger = logging.getLogger(__name__)


def log_transient_modified_system(
    K_mod: csr_matrix,
    M_mod: csr_matrix,
    C_mod: Optional[csr_matrix],
    *,
    n_bc_dofs: int,
    job_results_dir: Optional[Union[str, Path]] = None,
) -> None:
    """Log Frobenius norms and BC counts after dynamic BC application."""

    def _frob(a: csr_matrix) -> float:
        if a is None or a.nnz == 0:
            return 0.0
        return float(np.sqrt(np.sum(np.asarray(a.data, dtype=np.float64) ** 2)))

    _logger.info(
        "Transient modified system: nnz(K)=%s nnz(M)=%s frob(K)=%.3e frob(M)=%.3e "
        "C_present=%s frob(C)=%.3e n_bc=%s",
        getattr(K_mod, "nnz", 0),
        getattr(M_mod, "nnz", 0),
        _frob(K_mod),
        _frob(M_mod),
        C_mod is not None and getattr(C_mod, "nnz", 0) > 0,
        _frob(C_mod) if C_mod is not None else 0.0,
        n_bc_dofs,
    )
    if job_results_dir and os.path.isdir(str(job_results_dir)):
        # ``job_results_dir`` is typically ``<job_root>/primary_results`` for runners.
        root = Path(job_results_dir).resolve().parent
        log_path = root / "logs" / "transient_run_diagnostic.log"
        try:
            log_path.parent.mkdir(parents=True, exist_ok=True)
            with open(log_path, "a", encoding="utf-8") as fh:
                fh.write(
                    f"K nnz={getattr(K_mod, 'nnz', 0)} M nnz={getattr(M_mod, 'nnz', 0)} "
                    f"frob_K={_frob(K_mod):.6e} frob_M={_frob(M_mod):.6e} "
                    f"frob_C={_frob(C_mod) if C_mod is not None else 0.0:.6e} n_bc={n_bc_dofs}\n"
                )
        except OSError:
            pass

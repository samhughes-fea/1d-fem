# processing/spectral/spectral_diagnostics.py
"""Diagnostics for assembled spectral matrices (norms, dimensions, nnz)."""

import logging
import os

import numpy as np

_logger = logging.getLogger(__name__)


def log_spectral_constrained_dofs(bc_dofs, total_dof: int | None = None, job_results_dir=None):
    """Log penalty / fixed DOF count after spectral matrix modification."""
    n_bc = int(np.asarray(bc_dofs).size) if bc_dofs is not None else 0
    _logger.info(
        "Spectral BC: %s constrained global DOF(s)%s",
        n_bc,
        f", total_dof={total_dof}" if total_dof is not None else "",
    )
    if job_results_dir and os.path.isdir(str(job_results_dir)):
        log_path = os.path.join(os.path.dirname(str(job_results_dir)), "logs", "spectral_bc_diagnostic.log")
        os.makedirs(os.path.dirname(log_path), exist_ok=True)
        try:
            with open(log_path, "a", encoding="utf-8") as f:
                f.write(f"constrained_dofs={n_bc}\n")
                if total_dof is not None:
                    f.write(f"total_dof={int(total_dof)}\n")
        except OSError:
            pass


def log_spectral_diagnostics(K, M, job_results_dir=None):
    """Log basic diagnostics for global K and M (shape, nnz, Frobenius norm of data)."""
    if K is not None:
        _logger.info(
            "Spectral K: shape=%s, nnz=%s, norm=%.3e",
            K.shape,
            K.nnz,
            (K.data**2).sum() ** 0.5 if K.nnz else 0.0,
        )
    if M is not None:
        _logger.info(
            "Spectral M: shape=%s, nnz=%s, norm=%.3e",
            M.shape,
            M.nnz,
            (M.data**2).sum() ** 0.5 if M.nnz else 0.0,
        )
    if job_results_dir and os.path.isdir(job_results_dir):
        log_path = os.path.join(os.path.dirname(job_results_dir), "logs", "spectral_diagnostic.log")
        os.makedirs(os.path.dirname(log_path), exist_ok=True)
        try:
            with open(log_path, "a", encoding="utf-8") as f:
                if K is not None:
                    f.write(f"K shape={K.shape} nnz={K.nnz}\n")
                if M is not None:
                    f.write(f"M shape={M.shape} nnz={M.nnz}\n")
        except OSError:
            pass

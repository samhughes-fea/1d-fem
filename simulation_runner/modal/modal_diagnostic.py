# simulation_runner/modal/modal_diagnostic.py
"""Diagnostics for modal matrices (norms, dimensions, nnz)."""

import logging
import os

_logger = logging.getLogger(__name__)


def log_modal_diagnostics(K, M, job_results_dir=None):
    """Log basic diagnostics for global K and M (shape, nnz, norm)."""
    if K is not None:
        _logger.info(
            "Modal K: shape=%s, nnz=%s, norm=%.3e",
            K.shape,
            K.nnz,
            (K.data ** 2).sum() ** 0.5 if K.nnz else 0.0,
        )
    if M is not None:
        _logger.info(
            "Modal M: shape=%s, nnz=%s, norm=%.3e",
            M.shape,
            M.nnz,
            (M.data ** 2).sum() ** 0.5 if M.nnz else 0.0,
        )
    if job_results_dir and os.path.isdir(job_results_dir):
        log_path = os.path.join(os.path.dirname(job_results_dir), "logs", "modal_diagnostic.log")
        os.makedirs(os.path.dirname(log_path), exist_ok=True)
        try:
            with open(log_path, "a", encoding="utf-8") as f:
                if K is not None:
                    f.write(f"K shape={K.shape} nnz={K.nnz}\n")
                if M is not None:
                    f.write(f"M shape={M.shape} nnz={M.nnz}\n")
        except OSError:
            pass

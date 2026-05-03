# processing/dynamic/operations/integrate_transient_system.py
"""Newmark time integration for M u'' + C u' + K u = F(t)."""

from __future__ import annotations

from pathlib import Path
from typing import Callable, Optional, Union

import numpy as np
from scipy.sparse import csr_matrix

from processing.dynamic.time_integration import newmark_integrate

from processing.dynamic.operations._logging import init_stage_logger


class IntegrateTransientSystem:
    """
    Wrap :func:`processing.dynamic.time_integration.newmark_integrate`.

    ``force_func`` maps time ``t`` (scalar) to global force vector length ``n_dof``.
    """

    def __init__(
        self,
        t_grid: np.ndarray,
        force_func: Callable[[float], np.ndarray],
        job_results_dir: Optional[Union[str, Path]] = None,
    ):
        self.t_grid = np.asarray(t_grid, dtype=np.float64).ravel()
        self.force_func = force_func
        self.job_results_dir = Path(job_results_dir) if job_results_dir else None
        self._log = init_stage_logger("IntegrateTransientSystem", self.job_results_dir)

    def run(
        self,
        K_mod: csr_matrix,
        M_mod: csr_matrix,
        C_mod: Optional[csr_matrix],
        u0: np.ndarray,
        v0: np.ndarray,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        self._log.info(
            "Newmark integration: %s time steps, n_dof=%s",
            self.t_grid.size,
            u0.size,
        )
        U, V, A = newmark_integrate(
            K_mod, M_mod, C_mod, u0, v0, self.t_grid, self.force_func
        )
        peak_u = float(np.max(np.abs(U)))
        peak_v = float(np.max(np.abs(V)))
        peak_a = float(np.max(np.abs(A)))
        drift = float(np.max(np.abs(U[-1] - U[0]))) if U.shape[0] > 1 else 0.0
        self._log.info(
            "Transient stability snapshot: max|u|=%.3e max|v|=%.3e max|a|=%.3e end-start|u|=%.3e",
            peak_u,
            peak_v,
            peak_a,
            drift,
        )
        self._log.info("Time integration complete")
        return U, V, A

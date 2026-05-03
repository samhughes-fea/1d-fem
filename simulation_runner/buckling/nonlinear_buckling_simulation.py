# simulation_runner/buckling/nonlinear_buckling_simulation.py
"""§5 Nonlinear buckling — MVP orchestration shell (solver deferred)."""

from __future__ import annotations

import logging
import os

logger = logging.getLogger(__name__)


class NonlinearBucklingSimulationRunner:
    """
    Placeholder runner for future **nonlinear** buckling / path-following analysis.

    This is **not** linearized :math:`(K + \\lambda K_g)\\phi = 0` (see
    :class:`~simulation_runner.buckling.buckling_simulation.LinearBucklingSimulationRunner`).
    When ``[Buckling] nonlinear_buckling = true``, :func:`workflow_orchestrator.run_job.process_job`
    dispatches here. The MVP writes a diagnostics marker and exits successfully so CI can lock wiring.
    """

    def __init__(self, settings: dict, job_name: str):
        self.settings = settings
        self.job_name = job_name
        self.simulation_settings = settings.get("simulation_settings") or {}

    def run(self) -> None:
        job_results_dir = self.settings.get("job_results_dir")
        if not job_results_dir:
            raise ValueError("NonlinearBucklingSimulationRunner requires job_results_dir in settings")
        diag = os.path.join(job_results_dir, "diagnostics")
        os.makedirs(diag, exist_ok=True)
        marker = os.path.join(diag, "nonlinear_buckling_mvp_stub.txt")
        with open(marker, "w", encoding="utf-8") as f:
            f.write(
                "NonlinearBucklingSimulationRunner MVP: path-following / arc-length not implemented.\n"
                "See docs/conventions/NONLINEAR_BUCKLING_MVP.md\n"
            )
        logger.info(
            "Nonlinear buckling MVP stub wrote %s (solver not implemented).",
            marker,
        )

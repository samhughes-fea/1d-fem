# simulation_runner/buckling/buckling_simulation.py
"""§5 Linear buckling with prestress."""

from __future__ import annotations

import logging
import warnings

from simulation_runner.spectral.vibration_buckling_backend import VibrationBucklingBackend

logger = logging.getLogger(__name__)


class LinearBucklingSimulationRunner(VibrationBucklingBackend):
    """Linear buckling eigenproblem about a prestressed state (§5 linearized theory)."""

    def run(self):
        try:
            self.setup_simulation()
            self._run_buckling_analysis()
        except Exception as exc:
            logger.exception("Buckling simulation failed")
            raise RuntimeError("Buckling simulation aborted") from exc


class BucklingSimulationRunner(LinearBucklingSimulationRunner):
    """Deprecated alias for :class:`LinearBucklingSimulationRunner`."""

    def __init__(self, *args, **kwargs):
        warnings.warn(
            "BucklingSimulationRunner is deprecated; use LinearBucklingSimulationRunner.",
            DeprecationWarning,
            stacklevel=2,
        )
        super().__init__(*args, **kwargs)

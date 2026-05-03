# simulation_runner/buckling/buckling_simulation.py
"""§5 Linear buckling with prestress."""

import logging

from simulation_runner.spectral.vibration_buckling_backend import VibrationBucklingBackend

logger = logging.getLogger(__name__)


class BucklingSimulationRunner(VibrationBucklingBackend):
    """Linear buckling eigenproblem about a prestressed state."""

    def run(self):
        try:
            self.setup_simulation()
            self._run_buckling_analysis()
        except Exception as exc:
            logger.exception("Buckling simulation failed")
            raise RuntimeError("Buckling simulation aborted") from exc

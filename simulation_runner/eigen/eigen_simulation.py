# simulation_runner/eigen/eigen_simulation.py
"""§2 Eigen / undamped free vibration."""

import logging

from simulation_runner.spectral.vibration_buckling_backend import VibrationBucklingBackend

logger = logging.getLogger(__name__)


class EigenSimulationRunner(VibrationBucklingBackend):
    """Natural frequencies and mode shapes (same physics pipeline as legacy modal vibration)."""

    def run(self):
        try:
            self.setup_simulation()
            self._run_vibration_analysis()
        except Exception as exc:
            logger.exception("Eigen simulation failed")
            raise RuntimeError("Eigen simulation aborted") from exc

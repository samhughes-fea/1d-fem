# simulation_runner/eigen/eigen_simulation.py
"""§2 Eigen / undamped free vibration."""

import logging

from onedfem.simulation_runner.spectral.vibration_buckling_backend import VibrationBucklingBackend
from onedfem.stage_log import stage

logger = logging.getLogger(__name__)


class EigenSimulationRunner(VibrationBucklingBackend):
    """Natural frequencies and mode shapes (same physics pipeline as legacy modal vibration)."""

    def run(self):
        """Solve the blade beam and return its vibration modes.

        The return is the eigenpairs this run stored. It does not assemble a
        rotor and it does not solve a second eigenproblem.
        """
        try:
            with stage("onedfem", "beam_vibration"):
                self.setup_simulation()
                self._run_vibration_analysis()
        except Exception as exc:
            logger.exception("Eigen simulation failed")
            raise RuntimeError("Eigen simulation aborted") from exc
        from onedfem.post_processing.limit_states import blade_vibration_modes

        return blade_vibration_modes(self)

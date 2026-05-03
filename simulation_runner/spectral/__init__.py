# simulation_runner/spectral/__init__.py
"""
Shared spectral analysis backend for eigen vibration and linear buckling (analysis sections 2 and 5).

``EigenSimulationRunner`` and ``BucklingSimulationRunner`` subclass
:class:`VibrationBucklingBackend` from ``vibration_buckling_backend``.
"""

from simulation_runner.spectral.spectral_diagnostics import (
    log_spectral_constrained_dofs,
    log_spectral_diagnostics,
)
from simulation_runner.spectral.vibration_buckling_backend import VibrationBucklingBackend

__all__ = [
    "VibrationBucklingBackend",
    "log_spectral_constrained_dofs",
    "log_spectral_diagnostics",
]

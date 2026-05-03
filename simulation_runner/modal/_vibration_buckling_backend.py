"""Deprecated shim — import from ``simulation_runner.spectral`` instead."""

from __future__ import annotations

import warnings

warnings.warn(
    "simulation_runner.modal._vibration_buckling_backend is deprecated; use "
    "simulation_runner.spectral.vibration_buckling_backend (or simulation_runner.spectral.VibrationBucklingBackend).",
    DeprecationWarning,
    stacklevel=2,
)

from simulation_runner.spectral.vibration_buckling_backend import VibrationBucklingBackend

__all__ = ["VibrationBucklingBackend"]

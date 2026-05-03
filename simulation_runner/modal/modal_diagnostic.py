"""Deprecated shim — import from ``simulation_runner.spectral.spectral_diagnostics``."""

from __future__ import annotations

import warnings

warnings.warn(
    "simulation_runner.modal.modal_diagnostic is deprecated; use "
    "simulation_runner.spectral.spectral_diagnostics.log_spectral_diagnostics.",
    DeprecationWarning,
    stacklevel=2,
)

from simulation_runner.spectral.spectral_diagnostics import log_spectral_diagnostics as log_modal_diagnostics

__all__ = ["log_modal_diagnostics"]

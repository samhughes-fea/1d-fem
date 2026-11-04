# simulation_runner/spectral/spectral_diagnostics.py
"""Re-export spectral diagnostics from :mod:`onedfem.processing.spectral` (single source of truth)."""

from onedfem.processing.spectral.spectral_diagnostics import (
    log_spectral_constrained_dofs,
    log_spectral_diagnostics,
)

__all__ = ["log_spectral_constrained_dofs", "log_spectral_diagnostics"]

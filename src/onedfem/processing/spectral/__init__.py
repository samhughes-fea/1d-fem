# processing/spectral/__init__.py
"""
Spectral analysis helpers (Sections 2 and 5): diagnostics and staged operations.

Stage classes for eigen vibration and linear buckling live in
:mod:`onedfem.processing.spectral.operations`. Low-level scatter/BC kernels remain in
:mod:`onedfem.processing.eigen` and :mod:`onedfem.processing.buckling`.
"""

from onedfem.processing.spectral.spectral_diagnostics import (
    log_spectral_constrained_dofs,
    log_spectral_diagnostics,
)

__all__ = ["log_spectral_constrained_dofs", "log_spectral_diagnostics"]

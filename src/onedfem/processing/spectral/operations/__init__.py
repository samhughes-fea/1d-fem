# processing/spectral/operations/__init__.py
"""Staged spectral pipeline (eigen vibration + linear buckling matrix steps)."""

from onedfem.processing.spectral.operations.assemble_spectral_global import AssembleSpectralGlobalSystem
from onedfem.processing.spectral.operations.buckling_stages import (
    AssembleBucklingGeometricStiffness,
    ModifyBucklingGlobalMatrices,
    SolveLinearBucklingEigenpairs,
)
from onedfem.processing.spectral.operations.modify_spectral_global import ModifySpectralGlobalSystem
from onedfem.processing.spectral.operations.prepare_spectral_local import PrepareSpectralLocalMatrices
from onedfem.processing.spectral.operations.solve_generalized_eigenproblem import SolveGeneralizedEigenproblem

__all__ = [
    "PrepareSpectralLocalMatrices",
    "AssembleSpectralGlobalSystem",
    "ModifySpectralGlobalSystem",
    "SolveGeneralizedEigenproblem",
    "AssembleBucklingGeometricStiffness",
    "ModifyBucklingGlobalMatrices",
    "SolveLinearBucklingEigenpairs",
]

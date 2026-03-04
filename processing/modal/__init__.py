# processing/modal/__init__.py
"""Modal analysis: assembly of K/M and boundary conditions. No imports from processing.static."""

from processing.modal.assembly import assemble_global_matrices
from processing.modal.boundary_conditions import apply_boundary_conditions

__all__ = ["assemble_global_matrices", "apply_boundary_conditions"]

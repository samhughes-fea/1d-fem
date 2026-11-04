# processing/eigen/__init__.py
"""
Eigenproblems and spectral assembly helpers (taxonomy §2).

Global **K** / **M** assembly and penalty BCs for undamped vibration live in
``onedfem.processing.eigen.assembly`` and ``onedfem.processing.eigen.boundary_conditions``.
§5 linear buckling reuses the same scatter/BC utilities via ``onedfem.processing.buckling``.
Smallest generalized eigenpairs (**K** x = λ **M** x) for sparse §2 pencils:
``onedfem.processing.eigen.smallest_generalized_eigenpairs``.
"""

from onedfem.processing.eigen.assembly import (
    _compute_local_global_dof_map,
    _scatter_element_matrix,
    assemble_global_matrices,
)
from onedfem.processing.eigen.boundary_conditions import PENALTY, apply_boundary_conditions

__all__ = [
    "PENALTY",
    "_compute_local_global_dof_map",
    "_scatter_element_matrix",
    "assemble_global_matrices",
    "apply_boundary_conditions",
]

# processing/dynamic/__init__.py
"""Dynamic analysis: assembly, boundary conditions, time integration. No imports from processing.static."""

from processing.dynamic.assembly import assemble_global_system
from processing.dynamic.boundary_conditions import apply_boundary_conditions
from processing.dynamic.time_integration import newmark_integrate

__all__ = ["assemble_global_system", "apply_boundary_conditions", "newmark_integrate"]

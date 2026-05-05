"""Transient analysis: canonical Section 3 package name."""

from processing.transient.assembly import assemble_global_system
from processing.transient.boundary_conditions import apply_boundary_conditions
from processing.transient.time_integration import newmark_integrate

__all__ = [
    "assemble_global_system",
    "apply_boundary_conditions",
    "newmark_integrate",
]

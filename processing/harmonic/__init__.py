# processing/harmonic/__init__.py
"""
Frequency-domain dynamics (taxonomy §4): harmonic / steady-state response.

Kernels live in ``frequency_response``; see ``docs/conventions/HARMONIC_FREQUENCY_DOMAIN.md``
and the pointer ``docs/element_library/harmonic_frequency_domain_design.md``.
"""

from processing.harmonic.frequency_response import (
    frequency_grid_hz,
    harmonic_damping_matrix,
    mass_proportional_damping,
    rayleigh_damping,
    solve_one_frequency,
    sweep_displacements,
)
from processing.harmonic.modal_superposition import (
    harmonic_displacement_modal_superposition,
    undamped_natural_modes,
)

__all__ = [
    "frequency_grid_hz",
    "harmonic_damping_matrix",
    "harmonic_displacement_modal_superposition",
    "mass_proportional_damping",
    "rayleigh_damping",
    "solve_one_frequency",
    "sweep_displacements",
    "undamped_natural_modes",
]

# processing\static\results\compute_tertiary\__init__.py

"""
Tertiary results computation module.

This module computes highly derived engineering quantities from
secondary results, including section forces, principal stresses,
and failure criteria.
"""

from .section_force import ComputeSectionForce
from .principal_stress import ComputePrincipalStress
from .compute_tertiary_results import TertiaryResultsOrchestrator

__all__ = [
    "ComputeSectionForce",
    "ComputePrincipalStress",
    "TertiaryResultsOrchestrator",
]


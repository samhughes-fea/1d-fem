"""Staged transient dynamics pipeline."""

from processing.transient.operations.assemble_transient_global import AssembleTransientGlobalSystem
from processing.transient.operations.integrate_transient_system import IntegrateTransientSystem
from processing.transient.operations.modify_transient_global import ModifyTransientGlobalSystem

AssembleDynamicGlobalSystem = AssembleTransientGlobalSystem
ModifyDynamicGlobalSystem = ModifyTransientGlobalSystem

__all__ = [
    "AssembleTransientGlobalSystem",
    "ModifyTransientGlobalSystem",
    "IntegrateTransientSystem",
    "AssembleDynamicGlobalSystem",
    "ModifyDynamicGlobalSystem",
]

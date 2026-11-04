"""Staged transient dynamics pipeline."""

from onedfem.processing.transient.operations.assemble_transient_global import AssembleTransientGlobalSystem
from onedfem.processing.transient.operations.integrate_transient_system import IntegrateTransientSystem
from onedfem.processing.transient.operations.modify_transient_global import ModifyTransientGlobalSystem

AssembleDynamicGlobalSystem = AssembleTransientGlobalSystem
ModifyDynamicGlobalSystem = ModifyTransientGlobalSystem

__all__ = [
    "AssembleTransientGlobalSystem",
    "ModifyTransientGlobalSystem",
    "IntegrateTransientSystem",
    "AssembleDynamicGlobalSystem",
    "ModifyDynamicGlobalSystem",
]

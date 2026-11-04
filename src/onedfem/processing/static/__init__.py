"""Static FEM processing: assembly, boundary conditions, condensation, solve, results."""

from onedfem.processing.static.operations.assembly import AssembleGlobalSystem
from onedfem.processing.static.operations.condensation import CondenseModifiedSystem
from onedfem.processing.static.operations.disassembly import DisassembleGlobalSystem
from onedfem.processing.static.operations.modification import ModifyGlobalSystem
from onedfem.processing.static.operations.preparation import PrepareLocalSystem
from onedfem.processing.static.operations.reconstruction import ReconstructGlobalSystem
from onedfem.processing.static.operations.solver import SolveCondensedSystem

__all__ = [
    "PrepareLocalSystem",
    "AssembleGlobalSystem",
    "ModifyGlobalSystem",
    "CondenseModifiedSystem",
    "SolveCondensedSystem",
    "ReconstructGlobalSystem",
    "DisassembleGlobalSystem",
]

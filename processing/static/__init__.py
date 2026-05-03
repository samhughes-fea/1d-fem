"""Static FEM processing: assembly, boundary conditions, condensation, solve, results."""

from processing.static.operations.assembly import AssembleGlobalSystem
from processing.static.operations.condensation import CondenseModifiedSystem
from processing.static.operations.disassembly import DisassembleGlobalSystem
from processing.static.operations.modification import ModifyGlobalSystem
from processing.static.operations.preparation import PrepareLocalSystem
from processing.static.operations.reconstruction import ReconstructGlobalSystem
from processing.static.operations.solver import SolveCondensedSystem

__all__ = [
    "PrepareLocalSystem",
    "AssembleGlobalSystem",
    "ModifyGlobalSystem",
    "CondenseModifiedSystem",
    "SolveCondensedSystem",
    "ReconstructGlobalSystem",
    "DisassembleGlobalSystem",
]

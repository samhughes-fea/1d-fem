# tosdt = third-order shear deformation theory

Shared operators live in `utilities/` (Hermite **N**, **B**, **D**). Reddy is variational;
Levinson is Petrov–Galerkin. Both are spanwise 1D beams (`third_order_beam.py` orchestrator).
The Timoshenko shape-function registry entry is only for parent-class 12-DOF hooks;
stiffness uses the Hermite operators in `utilities/`.

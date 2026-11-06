# Nonlinear tosdt (third-order shear deformation theory)

Total Lagrangian and Updated Lagrangian Reddy and Levinson beams (16 DOF). Linear
operators come from `linear/beam/tosdt/utilities/`. Green–Lagrange on the 12
standard DOFs (TL) comes from Timoshenko TL. Reddy is variational; Levinson is
Petrov–Galerkin.

UL lives under `{reddy,levinson}/updated_lagrangian/`: current-chord linear
`K_0` plus beam-column `K_sigma`, with the runner calling
`update_reference_configuration` after each converged increment.

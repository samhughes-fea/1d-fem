# Linear Levinson and Reddy beam elements (3D)

Both elements are `ThirdOrderShearBeamElement3D` (`../third_order_beam.py`): Timoshenko
section motion plus the Saint-Venant shear warping slaved to the shear strain, with
cubic Hermite `v`, `w` and their slopes as two extra nodal DOFs (`Q1 = v'`, `Q2 = w'`).
Reddy is the variational weak form; Levinson keeps only the classical equilibrium
equations (Petrov-Galerkin).

Derivation, section terms, boundary conditions and checks:
[docs/element_library/higher_order_shear.md](../../../../../../docs/element_library/higher_order_shear.md).

The `utilities/` operators in this folder belong to the previous formulation and are
no longer used by the elements.

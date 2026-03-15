# pre_processing/element_library/nonlinear/gebt_shear/__init__.py
"""
Shear-deformable GEBT (Geometrically Exact Beam Theory) element package.

Phase 3a: 2-node 3D beam with shear deformation; tangent stiffness K_T(U_e)
and internal force F_int(U_e). At U_e=0, K_T equals linear Timoshenko K_e.
"""

from pre_processing.element_library.nonlinear.gebt_shear.gebt_shear_3D import (
    GEBTShearBeamElement3D,
)

__all__ = ["GEBTShearBeamElement3D"]

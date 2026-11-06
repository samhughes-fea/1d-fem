# pre_processing/element_library/nonlinear/timoshenko/total_lagrangian/utilities/geometric_stiffness.py
"""
\\(\\mathbf{K}_\\sigma\\) for nonlinear Timoshenko — beam-column Gauss form (**re-export of EB operator**).

**When to use ``GeometricStiffnessOperator.assemble_K_sigma``:**
axial ``N`` and moments ``M_y``, ``M_z`` acting on ``dN/dx`` (beam-column / linearised prestress).
Used by TL Euler–Bernoulli, GESDB ``native``, and linear buckling prestress twins.

**When the Hessian path is required (TL Timoshenko main tangent):**
``K_sigma`` on ``NonlinearTimoshenkoBeamElement3D`` is ``sum_k S_k d²E_k/dU²``
(``_geometric_stiffness_at``). Green–Lagrange shear/rotation remainders in ``E_nl``
(Voigt rows 3–4 when ``include_shear=True``, plus axial–curvature coupling on κ)
are **not** contained in the ``N, M_y, M_z`` beam-column operator. Selective
integration of the shear block of ``D`` also means ``S`` in the Hessian is
split by ``_material_rules``; the operator form has no shear-block split.

The operator instance on the Timoshenko element is kept for the beam-column
prestress path (``linear_geometric_stiffness_matrix``), not the Newton ``K_T``.
"""
from onedfem.pre_processing.element_library.nonlinear.euler_bernoulli.total_lagrangian.utilities.geometric_stiffness import (
    GeometricStiffnessOperator,
)

__all__ = ["GeometricStiffnessOperator"]

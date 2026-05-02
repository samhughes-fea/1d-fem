# pre_processing/element_library/nonlinear/timoshenko/utilities/stress_resultant.py
"""
Stress resultants — **Total Lagrangian** Timoshenko beam (**St. Venant–Kirchhoff-style** beam reduction).

**Governing constitutive (Voigt):** \\(\\mathbf{S} = \\mathbf{D}\\,\\mathbf{E}\\) with \\(\\mathbf{E}\\) Green–Lagrange-type from ``green_lagrange_strain``;
\\(\\mathbf{S}\\) **2PK-type** section resultants (6,). ``section_forces_from_strain`` returns ``(N, M_y, M_z)`` for ``GeometricStiffnessOperator`` (\\(\\mathbf{K}_\\sigma\\)); full ``S`` includes shear and torsion rows.
"""

from dataclasses import dataclass
from typing import Tuple

import numpy as np


@dataclass(frozen=True)
class StressResultantOperator:
    """
    Section forces (**2PK-type resultants**) from Green–Lagrange strain \\(\\mathbf{E}\\) in the **reference** configuration.

    **Governing relation:** \\(\\mathbf{S} = \\mathbf{D}\\,\\mathbf{E}\\).

    S = D @ E  at each Gauss point. Full 6-component view: [N, M_y, M_z, V_y, V_z, T]. For Timoshenko,
    D has non-zero shear rows so V_y, V_z and T are produced by D @ E; this method returns N, M_y, M_z;
    the complete data structure holds (N, M_y, M_z, V_y, V_z, T) for all formulations.

    Parameters
    ----------
    None (stateless; D and E are passed to methods).

    Notes
    -----
    Same API as EB ``StressResultantOperator``; Timoshenko ``D`` gives non-zero ``V_y``, ``V_z`` in full ``S`` but this method only exposes ``N``, ``M_y``, ``M_z`` for geometric stiffness.

    See Also
    --------
    nonlinear_timoshenko_3D.NonlinearTimoshenkoBeamElement3D
    nonlinear.euler_bernoulli.utilities.geometric_stiffness.GeometricStiffnessOperator
    """

    def section_forces_from_strain(
        self,
        E: np.ndarray,
        D: np.ndarray,
    ) -> Tuple[float, float, float]:
        """
        Compute section force resultants from strain and material stiffness at one point.

        For a single Gauss point or average: N = (D @ E)[0], M_y = (D @ E)[1], M_z = (D @ E)[2].
        Full integration over the element would sum over Gauss points with weights.

        Parameters
        ----------
        E : np.ndarray, shape (6,)
            Green–Lagrange-type strain vector **E** (Voigt; finite-strain Timoshenko entries).
        D : np.ndarray, shape (6, 6)
            Material stiffness (same as linear D matrix).

        Returns
        -------
        N : float
            Axial force.
        M_y : float
            Bending moment about y.
        M_z : float
            Bending moment about z.
        """
        S = D @ E
        N = float(S[0])
        M_y = float(S[1])
        M_z = float(S[2])
        return N, M_y, M_z

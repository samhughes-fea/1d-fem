# pre_processing/element_library/nonlinear/euler_bernoulli/utilities/stress_resultant.py
"""
Section force (stress resultant) operator for 2-node 3D Euler–Bernoulli beam (Total Lagrangian).
S = D @ E; N, M_y, M_z at Gauss points.
"""

from dataclasses import dataclass
from typing import Tuple

import numpy as np


@dataclass(frozen=True)
class StressResultantOperator:
    """
    Section forces (stress resultants) from 2nd Piola–Kirchhoff stress in reference configuration.

    S = D @ E  at each Gauss point. Full 6-component view: [N, M_y, M_z, V_y, V_z, T].
    For Euler–Bernoulli, D has rows 4–5 zero so V_y = V_z = 0 from constitutive; shear from
    equilibrium (V = dM/dx) if needed. This operator returns N, M_y, M_z; T and V_y, V_z
    are in the complete data structure (zero for EB from D).

    Parameters
    ----------
    None (stateless; D and E are passed to methods).
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
            Strain vector [ε_x, κ_y, κ_z, γ_xy, γ_xz, φ_x].
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

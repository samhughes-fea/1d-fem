# pre_processing/element_library/nonlinear/timoshenko/utilities/stress_resultant.py
"""
Stress resultants for TL Timoshenko beam: ``S = D @ E`` with ``E``, ``S`` (6,) Voigt.

``section_forces_from_strain`` returns ``(N, M_y, M_z)`` = ``(S[0], S[1], S[2])`` for ``K_sigma`` assembly; full ``S`` retains Timoshenko shear and torsion rows.
"""

from dataclasses import dataclass
from typing import Tuple

import numpy as np


@dataclass(frozen=True)
class StressResultantOperator:
    """
    Section forces (stress resultants) from 2nd Piola–Kirchhoff stress in reference configuration.

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

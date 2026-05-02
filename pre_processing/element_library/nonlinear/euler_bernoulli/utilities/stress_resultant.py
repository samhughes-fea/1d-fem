# pre_processing/element_library/nonlinear/euler_bernoulli/utilities/stress_resultant.py
"""
Stress resultants — **Total Lagrangian** 2-node 3D Euler–Bernoulli beam (**St. Venant–Kirchhoff-style** beam reduction).

**Governing constitutive (Voigt):** \\(\\mathbf{S} = \\mathbf{D}\\,\\mathbf{E}\\) where \\(\\mathbf{S}\\) are **2nd Piola–Kirchhoff–type** section resultants
conjugate to the Green–Lagrange strain vector \\(\\mathbf{E}\\) from ``green_lagrange_strain`` (same \\(\\mathbf{D}\\) structure as linear EB, evaluated on \\(\\mathbf{E}\\)).

``section_forces_from_strain`` returns ``N``, ``M_y``, ``M_z`` for ``GeometricStiffnessOperator`` (\\(\\mathbf{K}_\\sigma\\)); full ``S`` (6,) follows Voigt order.
"""

from dataclasses import dataclass
from typing import Tuple

import numpy as np


@dataclass(frozen=True)
class StressResultantOperator:
    """
    Section forces (**2PK-type resultants**) from Green–Lagrange strain \\(\\mathbf{E}\\) in the **reference** configuration.

    **Governing relation:** \\(\\mathbf{S} = \\mathbf{D}\\,\\mathbf{E}\\) at each Gauss point (StVK beam reduction; not Cauchy \\(\\boldsymbol{\\sigma}\\) on infinitesimal \\(\\boldsymbol{\\varepsilon}\\)).

    S = D @ E  at each Gauss point. Full 6-component view: [N, M_y, M_z, V_y, V_z, T].
    For Euler–Bernoulli, D has rows 4–5 zero so V_y = V_z = 0 from constitutive; shear from
    equilibrium (V = dM/dx) if needed. This operator returns N, M_y, M_z; T and V_y, V_z
    are in the complete data structure (zero for EB from D).

    Parameters
    ----------
    None (stateless; D and E are passed to methods).

    Notes
    -----
    Stateless: does not integrate. The element supplies ``E`` and ``D`` per Gauss point; ``GeometricStiffnessOperator`` uses ``N``, ``M_y``, ``M_z``.
    EB ``D`` zeros shear rows — constitutive ``V_y``, ``V_z`` are zero; equilibrium shear is separate.

    See Also
    --------
    nonlinear_euler_bernoulli_3D.NonlinearEulerBernoulliBeamElement3D
    geometric_stiffness.GeometricStiffnessOperator
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
            Green–Lagrange-type strain vector **E** (Voigt; same row labels as infinitesimal ε but finite-strain values).
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

# processing_OOP\static\results\compute_tertiary\principal_stress.py

import numpy as np
from typing import List, Tuple
import logging

logger = logging.getLogger(__name__)


class ComputePrincipalStress:
    """
    Computes principal stresses and derived stress quantities from stress tensors.

    Principal stresses are the eigenvalues of the stress tensor and represent
    the maximum/minimum normal stresses acting on the material. They are critical
    for failure analysis and design verification.

    Also computes:
    - Von Mises equivalent stress (for ductile failure)
    - Maximum shear stress (Tresca criterion)

    Parameters
    ----------
    stress_gauss : List[List[np.ndarray]]
        Stress tensors at Gauss points for all elements
        Shape: List[element] -> List[gauss_point] -> np.ndarray(6,)
        Components: [σ_xx, σ_yy, σ_zz, τ_xy, τ_yz, τ_xz]
    """

    def __init__(self, stress_gauss: List[List[np.ndarray]]):
        self.stress_gauss = stress_gauss

    def compute_all(self) -> Tuple[
        List[List[np.ndarray]],  # principal stresses
        List[List[float]],        # von Mises stress
        List[List[float]]         # max shear stress
    ]:
        """
        Compute all derived stress quantities.

        Returns
        -------
        principal_stresses : List[List[np.ndarray]]
            Principal stress components [σ1, σ2, σ3] per Gauss point
            Shape: List[element] -> List[gauss_point] -> np.ndarray(3,)
        
        von_mises : List[List[float]]
            Von Mises equivalent stress per Gauss point
            Shape: List[element] -> List[gauss_point] -> float
        
        max_shear : List[List[float]]
            Maximum shear stress per Gauss point
            Shape: List[element] -> List[gauss_point] -> float
        """
        principal_stresses = []
        von_mises_stresses = []
        max_shear_stresses = []

        for elem_idx, elem_stresses in enumerate(self.stress_gauss):
            elem_principal = []
            elem_von_mises = []
            elem_max_shear = []

            for gp_idx, stress in enumerate(elem_stresses):
                stress = np.asarray(stress)
                # Principal / Von Mises / max shear require 6-component stress (full 3D tensor).
                # Bar (2) and truss (3) components: skip tensor computation; use placeholder values.
                if stress.size != 6:
                    principals = np.zeros(3, dtype=stress.dtype)
                    von_mises = 0.0
                    max_shear = 0.0
                else:
                    principals = self._compute_principal_stresses(stress)
                    von_mises = self._compute_von_mises(stress)
                    max_shear = self._compute_max_shear(principals)
                elem_principal.append(principals)
                elem_von_mises.append(von_mises)
                elem_max_shear.append(max_shear)

            principal_stresses.append(elem_principal)
            von_mises_stresses.append(elem_von_mises)
            max_shear_stresses.append(elem_max_shear)

        logger.info(f"✅ Computed principal stresses for {len(principal_stresses)} elements")
        return principal_stresses, von_mises_stresses, max_shear_stresses

    def _compute_principal_stresses(self, stress: np.ndarray) -> np.ndarray:
        """
        Compute principal stresses (eigenvalues of stress tensor).

        Parameters
        ----------
        stress : np.ndarray
            Stress vector [σ_xx, σ_yy, σ_zz, τ_xy, τ_yz, τ_xz]

        Returns
        -------
        np.ndarray
            Principal stresses [σ1, σ2, σ3] sorted in descending order
        """
        # Construct 3x3 stress tensor from Voigt notation
        σ_xx, σ_yy, σ_zz = stress[0], stress[1], stress[2]
        τ_xy, τ_yz, τ_xz = stress[3], stress[4], stress[5]

        stress_tensor = np.array([
            [σ_xx, τ_xy, τ_xz],
            [τ_xy, σ_yy, τ_yz],
            [τ_xz, τ_yz, σ_zz]
        ])

        # Compute eigenvalues (principal stresses)
        eigenvalues = np.linalg.eigvalsh(stress_tensor)
        
        # Sort in descending order: σ1 ≥ σ2 ≥ σ3
        principal_stresses = np.sort(eigenvalues)[::-1]
        
        return principal_stresses

    def _compute_von_mises(self, stress: np.ndarray) -> float:
        """
        Compute Von Mises equivalent stress.

        The Von Mises stress is used to predict yielding of ductile materials
        under complex loading:

            σ_vm = sqrt( 0.5 * [(σ1-σ2)² + (σ2-σ3)² + (σ3-σ1)²] )

        Or equivalently from stress components:

            σ_vm = sqrt( σ_xx² + σ_yy² + σ_zz² - σ_xx*σ_yy - σ_yy*σ_zz 
                        - σ_zz*σ_xx + 3*(τ_xy² + τ_yz² + τ_xz²) )

        Parameters
        ----------
        stress : np.ndarray
            Stress vector [σ_xx, σ_yy, σ_zz, τ_xy, τ_yz, τ_xz]

        Returns
        -------
        float
            Von Mises equivalent stress
        """
        σ_xx, σ_yy, σ_zz = stress[0], stress[1], stress[2]
        τ_xy, τ_yz, τ_xz = stress[3], stress[4], stress[5]

        von_mises = np.sqrt(
            σ_xx**2 + σ_yy**2 + σ_zz**2
            - σ_xx*σ_yy - σ_yy*σ_zz - σ_zz*σ_xx
            + 3.0 * (τ_xy**2 + τ_yz**2 + τ_xz**2)
        )

        return float(von_mises)

    def _compute_max_shear(self, principal_stresses: np.ndarray) -> float:
        """
        Compute maximum shear stress (Tresca criterion).

        For a 3D stress state:
            τ_max = (σ1 - σ3) / 2

        where σ1 ≥ σ2 ≥ σ3 are the principal stresses.

        Parameters
        ----------
        principal_stresses : np.ndarray
            Principal stresses [σ1, σ2, σ3] in descending order

        Returns
        -------
        float
            Maximum shear stress
        """
        σ1, σ3 = principal_stresses[0], principal_stresses[2]
        τ_max = abs(σ1 - σ3) / 2.0
        return float(τ_max)


class ComputeFailureIndex:
    """
    Computes material failure indices based on yield/ultimate strength criteria.
    
    The failure index represents the ratio of actual stress to allowable stress:
        FI = σ_equivalent / σ_allowable
    
    FI > 1.0 indicates failure/yielding.
    """

    def __init__(
        self,
        von_mises_stress: List[List[float]],
        yield_strength: float,
        safety_factor: float = 1.0
    ):
        """
        Parameters
        ----------
        von_mises_stress : List[List[float]]
            Von Mises stress at each Gauss point
        yield_strength : float
            Material yield strength
        safety_factor : float
            Design safety factor (default 1.0)
        """
        self.von_mises_stress = von_mises_stress
        self.allowable_stress = yield_strength / safety_factor

    def compute(self) -> List[List[float]]:
        """
        Compute failure index at all Gauss points.

        Returns
        -------
        List[List[float]]
            Failure index per Gauss point
            Shape: List[element] -> List[gauss_point] -> float
        """
        failure_indices = []

        for elem_von_mises in self.von_mises_stress:
            elem_failure = []
            for von_mises in elem_von_mises:
                fi = von_mises / self.allowable_stress
                elem_failure.append(float(fi))
            failure_indices.append(elem_failure)

        max_fi = max(max(elem) for elem in failure_indices if elem)
        logger.info(f"✅ Computed failure indices (max FI = {max_fi:.3f})")
        
        return failure_indices


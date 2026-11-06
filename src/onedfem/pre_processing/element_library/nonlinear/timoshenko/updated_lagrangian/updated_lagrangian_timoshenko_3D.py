# Updated-Lagrangian 3D Timoshenko beam (12 DOF).

"""Small-strain Timoshenko operators on the **current** chord (selective ``K_0``).
The runner calls ``update_reference_configuration`` after each converged increment.
"""

from __future__ import annotations

import numpy as np

from onedfem.pre_processing.element_library.linear.beam.fosdt.timoshenko.linear_timoshenko_3D import (
    LinearTimoshenkoBeamElement3D,
)
from onedfem.pre_processing.element_library.linear.beam.fosdt.timoshenko.utilities.B_matrix import (
    StrainDisplacementOperator,
)
from onedfem.pre_processing.element_library.nonlinear.euler_bernoulli.updated_lagrangian.utilities.current_chord import (
    CurrentChordFrameOperator,
    UpdatedLagrangianChordMixin,
)
from onedfem.pre_processing.element_library.nonlinear.timoshenko.updated_lagrangian.utilities.config_update import (
    ConfigurationUpdateOperator,
)
from onedfem.pre_processing.element_library.nonlinear.timoshenko.updated_lagrangian.utilities.current_strain import (
    CurrentConfigurationStrainOperator,
)
from onedfem.pre_processing.element_library.nonlinear.timoshenko.updated_lagrangian.utilities.geometric_stiffness import (
    GeometricStiffnessOperator,
)
from onedfem.pre_processing.element_library.nonlinear.timoshenko.updated_lagrangian.utilities.stress_resultant import (
    StressResultantOperator,
)
from onedfem.pre_processing.element_library.shape_function_registry import get_shape_function_operator


class UpdatedLagrangianTimoshenkoBeamElement3D(UpdatedLagrangianChordMixin, LinearTimoshenkoBeamElement3D):
    """UL Timoshenko: selective linear ``K_0`` on ``L_curr`` plus beam-column ``K_sigma``."""

    element_type_name = "Timoshenko-3D-UpdatedLagrangian"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.configuration_update_operator = ConfigurationUpdateOperator()
        self.current_chord_frame_operator = CurrentChordFrameOperator()
        self.current_strain_operator = CurrentConfigurationStrainOperator()
        self.ul_stress_resultant_operator = StressResultantOperator()
        self.ul_geometric_stiffness_operator = GeometricStiffnessOperator(element_length=self.L)
        self._U_ref = np.zeros(int(self._n_dof), dtype=np.float64)
        self._rebuild_ul_T()

    def _increment(self, U_e: np.ndarray) -> np.ndarray:
        u = np.asarray(U_e, dtype=np.float64).ravel()
        n = self._U_ref.size
        if u.size < n:
            u = np.pad(u, (0, n - u.size))
        return u[:n] - self._U_ref

    def _rebuild_current_operators(self) -> None:
        self.L = self.configuration_update_operator.chord_length(self.node_coords)
        self.shape_function_operator = get_shape_function_operator(
            "LinearTimoshenkoBeamElement3D", self.L
        )
        self.strain_displacement_operator = StrainDisplacementOperator(element_length=self.L)
        self.ul_geometric_stiffness_operator = GeometricStiffnessOperator(element_length=self.L)
        if hasattr(self, "_K_0"):
            self._K_0 = None
        self._rebuild_ul_T()

    def update_reference_configuration(self, U_e: np.ndarray) -> None:
        du = self._increment(U_e)
        self.node_coords = self.configuration_update_operator.apply_translations(self.node_coords, du)
        u = np.asarray(U_e, dtype=np.float64).ravel()
        n = self._U_ref.size
        self._U_ref = (u[:n] if u.size >= n else np.pad(u, (0, n - u.size))).copy()
        self._rebuild_current_operators()

    def _get_K_0_current(self) -> np.ndarray:
        return np.asarray(self.element_stiffness_matrix().K_e, dtype=np.float64)


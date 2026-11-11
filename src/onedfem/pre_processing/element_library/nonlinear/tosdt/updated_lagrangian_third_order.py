# Shared Updated-Lagrangian pipeline for Reddy / Levinson (16 DOF).

"""Linear tosdt Hermite ``B``/``D`` on the **current** chord; runner updates the reference."""

from __future__ import annotations

import numpy as np

from onedfem.pre_processing.element_library.linear.beam.tosdt.third_order_beam import (
    N_DOF,
    ThirdOrderShearBeamElement3D,
)
from onedfem.pre_processing.element_library.linear.beam.tosdt.utilities.B_matrix import (
    StrainDisplacementOperator as ThirdOrderStrainDisplacementOperator,
)
from onedfem.pre_processing.element_library.linear.beam.tosdt.utilities.shape_functions import (
    ShapeFunctionOperator as ThirdOrderShapeFunctionOperator,
)
from onedfem.pre_processing.element_library.nonlinear.euler_bernoulli.updated_lagrangian.utilities.config_update import (
    ConfigurationUpdateOperator,
)
from onedfem.pre_processing.element_library.nonlinear.euler_bernoulli.updated_lagrangian.utilities.current_chord import (
    CurrentChordFrameOperator,
    UpdatedLagrangianChordMixin,
)
from onedfem.pre_processing.element_library.nonlinear.euler_bernoulli.updated_lagrangian.utilities.current_strain import (
    CurrentConfigurationStrainOperator,
)
from onedfem.pre_processing.element_library.nonlinear.euler_bernoulli.updated_lagrangian.utilities.geometric_stiffness import (
    GeometricStiffnessOperator,
)
from onedfem.pre_processing.element_library.nonlinear.euler_bernoulli.updated_lagrangian.utilities.stress_resultant import (
    StressResultantOperator,
)


class UpdatedLagrangianThirdOrderShearBeamElement3D(UpdatedLagrangianChordMixin, ThirdOrderShearBeamElement3D):
    """UL tosdt: current-chord linear ``K_0`` plus beam-column ``K_sigma`` on 16 DOFs."""

    #: The tangent and internal force are not rotated by ``framed_element``: refuse rotated or
    #: curved meshes rather than inherit the linear parent's frame support.
    supports_element_frame = False

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.configuration_update_operator = ConfigurationUpdateOperator()
        self.current_chord_frame_operator = CurrentChordFrameOperator()
        self.current_strain_operator = CurrentConfigurationStrainOperator()
        self.ul_stress_resultant_operator = StressResultantOperator()
        self.ul_geometric_stiffness_operator = GeometricStiffnessOperator(element_length=self.L)
        self._U_ref = np.zeros(N_DOF, dtype=np.float64)
        self._rebuild_ul_T()

    def _increment(self, U_e: np.ndarray) -> np.ndarray:
        u = np.asarray(U_e, dtype=np.float64).ravel()
        if u.size < N_DOF:
            u = np.pad(u, (0, N_DOF - u.size))
        return u[:N_DOF] - self._U_ref

    def _rebuild_current_operators(self) -> None:
        self.L = self.configuration_update_operator.chord_length(self.node_coords)
        self.third_order_strain_displacement_operator = ThirdOrderStrainDisplacementOperator(
            element_length=self.L, rate=self._R
        )
        self.third_order_shape_function_operator = ThirdOrderShapeFunctionOperator(
            element_length=self.L
        )
        self.ul_geometric_stiffness_operator = GeometricStiffnessOperator(element_length=self.L)
        self._rebuild_ul_T()

    def update_reference_configuration(self, U_e: np.ndarray) -> None:
        du = self._increment(U_e)
        self.node_coords = self.configuration_update_operator.apply_translations(self.node_coords, du)
        u = np.asarray(U_e, dtype=np.float64).ravel()
        self._U_ref = (u[:N_DOF] if u.size >= N_DOF else np.pad(u, (0, N_DOF - u.size))).copy()
        self._rebuild_current_operators()

    def _get_K_0_current(self) -> np.ndarray:
        return np.asarray(self.element_stiffness_matrix().K_e, dtype=np.float64)


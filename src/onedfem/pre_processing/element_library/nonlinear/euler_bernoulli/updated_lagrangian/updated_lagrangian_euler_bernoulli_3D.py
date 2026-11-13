# Updated-Lagrangian 3D Euler-Bernoulli beam (12 DOF).

"""Small-strain operators on the **current** chord; after each converged increment the
runner calls ``update_reference_configuration`` so the last geometry is the new reference.

**Tensors:** ``eps = B U_e`` (linear EB ``B`` at ``L_curr``), ``S = D eps``,
``K_mat = int B.T D B dx``, ``K_sigma`` from ``GeometricStiffnessOperator`` /
linear geometric mixin, ``K_T = K_mat + K_sigma``, ``detJ = L_curr / 2``.
"""

from __future__ import annotations

import numpy as np

from onedfem.pre_processing.element_library.linear.beam.zosdt.euler_bernoulli.linear_euler_bernoulli_3D import (
    LinearEulerBernoulliBeamElement3D,
)
from onedfem.pre_processing.element_library.linear.beam.zosdt.euler_bernoulli.utilities.B_matrix import (
    StrainDisplacementOperator,
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
from onedfem.pre_processing.element_library.shape_function_registry import get_shape_function_operator


class UpdatedLagrangianEulerBernoulliBeamElement3D(UpdatedLagrangianChordMixin, LinearEulerBernoulliBeamElement3D):
    """UL Euler–Bernoulli: linear ``B``, ``D`` on the current chord plus beam-column ``K_sigma``."""

    #: The tangent and internal force are not rotated by ``framed_element``: refuse rotated or
    #: curved meshes rather than inherit the linear parent's frame support.
    supports_element_frame = False
    bending_interpolation = "hermite"

    element_type_name = "EulerBernoulli-3D-UpdatedLagrangian"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.configuration_update_operator = ConfigurationUpdateOperator()
        self.current_chord_frame_operator = CurrentChordFrameOperator()
        self.current_strain_operator = CurrentConfigurationStrainOperator()
        self.ul_stress_resultant_operator = StressResultantOperator()
        self.ul_geometric_stiffness_operator = GeometricStiffnessOperator(element_length=self.L)
        self._U_ref = np.zeros(12, dtype=np.float64)
        self._rebuild_ul_T()

    def _increment(self, U_e: np.ndarray) -> np.ndarray:
        return np.asarray(U_e, dtype=np.float64).ravel()[:12] - self._U_ref

    def _rebuild_current_operators(self) -> None:
        self.L = self.configuration_update_operator.chord_length(self.node_coords)
        self.shape_function_operator = get_shape_function_operator(
            "LinearEulerBernoulliBeamElement3D", self.L
        )
        self.strain_displacement_operator = StrainDisplacementOperator(element_length=self.L)
        self.ul_geometric_stiffness_operator = GeometricStiffnessOperator(element_length=self.L)
        self._rebuild_ul_T()

    def update_reference_configuration(self, U_e: np.ndarray) -> None:
        """Runner hook: bake the increment into the chord; later Newton steps use ``U - U_ref``."""
        du = self._increment(U_e)
        self.node_coords = self.configuration_update_operator.apply_translations(self.node_coords, du)
        self._U_ref = np.asarray(U_e, dtype=np.float64).ravel()[:12].copy()
        self._rebuild_current_operators()

    def _get_K_0_current(self) -> np.ndarray:
        return np.asarray(self.element_stiffness_matrix().K_e, dtype=np.float64)


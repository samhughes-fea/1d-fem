# pre_processing/element_library/nonlinear/geometrically_exact_shear_deformable_beam/geometrically_exact_shear_deformable_beam_3D.py
"""
Classical geometrically exact **shear-deformable** 3D beam (finite rotations, director kinematics) — STUB.

This is **not** the chord-based Total Lagrangian Timoshenko element; see
``NonlinearTimoshenkoBeamElement3D`` and ``assemble_timoshenko_K0`` for that stack.

Implementation milestones: literature weak form (Simo & Vu-Quoc–type), consistent tangent,
Gauss quadrature tied to that formulation. All assembly methods currently raise
``NotImplementedError``.
"""

from __future__ import annotations

import logging

import numpy as np

from pre_processing.element_library.element_1D_base import Element1DBase

logger = logging.getLogger(__name__)


class GeometricallyExactShearDeformableBeam3D(Element1DBase):
    """
    Placeholder registration type for a future Simo / Vu-Quoc–class shear-deformable GEBT element.

    Notes
    -----
    Factory type string: ``GeometricallyExactShearDeformableBeam3D``.
    Shape functions temporarily reuse the 2-node Timoshenko operator from the registry only
    so logging smoke-tests pass; kernels are not implemented.
    """

    element_type_name = "GeometricallyExactShearDeformable-3D"

    def __init__(
        self,
        *,
        element_id: int,
        element_dictionary: dict,
        grid_dictionary: dict,
        section_dictionary: dict,
        material_dictionary: dict,
        point_load_array: np.ndarray,
        distributed_load_array: np.ndarray,
        job_results_dir: str,
        quadrature_order: int | None = None,
    ):
        super().__init__(
            element_id=element_id,
            element_dictionary=element_dictionary,
            grid_dictionary=grid_dictionary,
            material_dictionary=material_dictionary,
            section_dictionary=section_dictionary,
            point_load_array=point_load_array,
            distributed_load_array=distributed_load_array,
            job_results_dir=job_results_dir,
            dof_per_node=6,
        )
        self.node_coords = self.grid_array
        self.L = float(np.linalg.norm(self.node_coords[1] - self.node_coords[0]))
        if quadrature_order is not None:
            self.quadrature_order = max(int(quadrature_order), 2)
        else:
            axial_order = int(self.element_array[3])
            bending_y_order = int(self.element_array[4])
            bending_z_order = int(self.element_array[5])
            shear_y_order = int(self.element_array[6])
            shear_z_order = int(self.element_array[7])
            torsion_order = int(self.element_array[8])
            load_order = int(self.element_array[9])
            self.quadrature_order = max(
                axial_order,
                bending_y_order,
                bending_z_order,
                shear_y_order,
                shear_z_order,
                torsion_order,
                load_order,
                2,
            )
        self._validate_element_properties()
        self._assert_logging_ready()
        logger.debug(
            "GeometricallyExactShearDeformableBeam3D id=%s stub initialised (L=%.4e)",
            self.element_id,
            self.L,
        )

    def _validate_element_properties(self) -> None:
        if self.L <= 0:
            raise ValueError(f"Invalid element length {self.L:.2e} for element {self.element_id}")
        if self.material_array.size != 4 or self.section_array.size not in (5, 7, 9, 10):
            raise ValueError("Material/section arrays not properly initialised")

    def element_stiffness_matrix(self):
        raise NotImplementedError(
            f"{self.__class__.__name__}: classical GEBT stiffness not implemented — see docs/element_library/"
            "geometrically_exact_shear_deformable_beam_formulation.md"
        )

    def element_force_vector(self):
        raise NotImplementedError(
            f"{self.__class__.__name__}: classical GEBT force vector not implemented."
        )

    def element_mass_matrix(self):
        raise NotImplementedError(
            f"{self.__class__.__name__}: classical GEBT mass matrix not implemented."
        )

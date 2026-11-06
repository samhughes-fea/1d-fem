# pre_processing/element_library/nonlinear/large_rotations/geometrically_exact_shear_deformable_beam/geometrically_exact_shear_deformable_beam_3D.py
"""
Geometrically exact **shear-deformable** beam registration (12-DOF Timoshenko-class kinematics).

**Kernel selection** (``simulation_settings['nonlinear']['gesdb_kernel']``, constructor):

- ``tl_locked`` — chord-frame Green–Lagrange Voigt (same stress path as
  :class:`NonlinearTimoshenkoBeamElement3D`) via :mod:`gesdb_kinematics`.
- ``native`` — engineering axial stretch plus linear Timoshenko bending/shear/torsion rows;
  :meth:`internal_force_vector` / :meth:`tangent_stiffness_matrix` use the consistent Jacobian
  ``B_eng`` (see ``docs/element_library/gesdb_weak_form.md``).

Optional ``gesdb_tl_fallback`` forces the parent TL strain hook regardless of ``gesdb_kernel``.

See ``docs/element_library/gesdb_weak_form.md``,
``docs/element_library/geometrically_exact_shear_deformable_beam_formulation.md``, and
``docs/element_library/large_rotation_vs_total_lagrangian.md``.
"""

from __future__ import annotations

import numpy as np

from onedfem.pre_processing.element_library.nonlinear.timoshenko.total_lagrangian.nonlinear_timoshenko_3D import (
    NonlinearTimoshenkoBeamElement3D,
)

_NL_STD = 12


class GeometricallyExactShearDeformableBeam3D(NonlinearTimoshenkoBeamElement3D):
    """
    Factory-facing **GeometricallyExactShearDeformableBeam3D** element type.

    Notes
    -----
    Strains at Gauss points follow :mod:`gesdb_kinematics` unless ``gesdb_tl_fallback`` is set.
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
        quadrature_order: int | None = 3,
        gesdb_tl_fallback: bool = False,
        gesdb_kernel: str = "tl_locked",
    ):
        super().__init__(
            element_id=element_id,
            element_dictionary=element_dictionary,
            grid_dictionary=grid_dictionary,
            section_dictionary=section_dictionary,
            material_dictionary=material_dictionary,
            point_load_array=point_load_array,
            distributed_load_array=distributed_load_array,
            job_results_dir=job_results_dir,
            quadrature_order=quadrature_order,
        )
        self._gesdb_tl_fallback = bool(gesdb_tl_fallback)
        gk = str(gesdb_kernel).strip().lower().replace("-", "_")
        if gk in ("tl_locked", "tl", "locked", "chord_tl"):
            self._gesdb_kernel = "tl_locked"
        elif gk in ("native", "engineering", "gesdb_native"):
            self._gesdb_kernel = "native"
        else:
            raise ValueError(
                "gesdb_kernel must be 'tl_locked' or 'native', "
                f"got {gesdb_kernel!r}"
            )

    def _gesdb_native_active(self) -> bool:
        return (
            not self._gesdb_tl_fallback
            and self._n_dof == _NL_STD
            and self._gesdb_kernel == "native"
        )

    def _tl_voigt_strain_at_gauss(self, U_e: np.ndarray, xi_g: float) -> np.ndarray:
        """GESDB pathway; warping meshes defer to the TL chord-frame strain unchanged."""
        if self._gesdb_tl_fallback or self._n_dof != _NL_STD:
            return super()._tl_voigt_strain_at_gauss(U_e, xi_g)
        if self._gesdb_kernel == "native":
            from onedfem.pre_processing.element_library.nonlinear.large_rotations.geometrically_exact_shear_deformable_beam.utilities.gesdb_kinematics import (
                native_engineering_voigt_strain_timoshenko_12,
            )

            return native_engineering_voigt_strain_timoshenko_12(self, U_e, xi_g)
        from onedfem.pre_processing.element_library.nonlinear.large_rotations.geometrically_exact_shear_deformable_beam.utilities.gesdb_kinematics import (
            gesdb_director_voigt_strain_timoshenko_12,
        )

        return gesdb_director_voigt_strain_timoshenko_12(self, U_e, xi_g)

    def internal_force_vector(self, U_e: np.ndarray) -> np.ndarray:
        """Consistent ``B_engᵀ S`` Gauss sum for the native kernel; else TL parent."""
        if not self._gesdb_native_active():
            return super().internal_force_vector(U_e)
        if self._n_dof != _NL_STD:
            return super().internal_force_vector(U_e)

        from onedfem.pre_processing.element_library.nonlinear.large_rotations.geometrically_exact_shear_deformable_beam.utilities.native_weak_form import (
            native_internal_force,
        )

        return native_internal_force(self, U_e)

    def tangent_stiffness_matrix(self, U_e: np.ndarray) -> np.ndarray:
        """``K_mat + K_sigma`` for native (no TL ``K_delta`` split); else parent TL tangent."""
        if not self._gesdb_native_active():
            return super().tangent_stiffness_matrix(U_e)
        if self._n_dof != _NL_STD:
            return super().tangent_stiffness_matrix(U_e)

        from onedfem.pre_processing.element_library.nonlinear.large_rotations.geometrically_exact_shear_deformable_beam.utilities.native_weak_form import (
            native_tangent_stiffness,
        )

        return native_tangent_stiffness(self, U_e)

    def _gesdb_strain_voigt_at_gauss(self, U_e: np.ndarray, xi_g: float) -> np.ndarray:
        """
        Voigt strain vector at natural coordinate ``xi_g`` (tests and diagnostics).

        Uses the same hook as the Gauss loops for ``F_int`` / ``K_T``.
        """
        U_e = np.asarray(U_e, dtype=np.float64).ravel()
        return np.asarray(self._tl_voigt_strain_at_gauss(U_e, float(xi_g)), dtype=np.float64).ravel()

"""GESDB milestone 1: strain scaffold vs TL parent and small-strain linear limit."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))


def _minimal_linear_ts_mesh():
    """Shared dictionaries for a single straight Timoshenko element."""
    grid_dictionary = {
        "ids": np.array([0, 1], dtype=np.int32),
        "coordinates": np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]], dtype=np.float64),
    }
    element_dictionary = {
        "ids": np.array([0], dtype=np.int32),
        "connectivity": np.array([[0, 1]], dtype=np.int32),
        "types": np.array(["LinearTimoshenkoBeamElement3D"]),
        "integration_orders": {
            "axial": np.array([3], dtype=np.int64),
            "bending_y": np.array([3], dtype=np.int64),
            "bending_z": np.array([3], dtype=np.int64),
            "shear_y": np.array([2], dtype=np.int64),
            "shear_z": np.array([2], dtype=np.int64),
            "torsion": np.array([3], dtype=np.int64),
            "load": np.array([2], dtype=np.int64),
        },
    }
    material_dictionary = {
        "E": np.array([210e9]),
        "G": np.array([80e9]),
        "nu": np.array([0.3]),
        "rho": np.array([7850.0]),
    }
    section_dictionary = {
        "A": np.array([0.01]),
        "I_x": np.array([1e-9]),
        "I_y": np.array([1e-6]),
        "I_z": np.array([1e-6]),
        "J_t": np.array([2e-9]),
        "kappa": np.array([5.0 / 6.0]),
        "alpha": np.array([0.0]),
    }
    return grid_dictionary, element_dictionary, material_dictionary, section_dictionary


def test_gesdb_strain_matches_gp_tl_kinematics():
    import tempfile
    import os
    import shutil
    from pre_processing.element_library.nonlinear.large_rotations.geometrically_exact_shear_deformable_beam.geometrically_exact_shear_deformable_beam_3D import (
        GeometricallyExactShearDeformableBeam3D,
    )
    from pre_processing.element_library.nonlinear.timoshenko.nonlinear_timoshenko_3D import (
        NonlinearTimoshenkoBeamElement3D,
    )

    gd, ed, md, sd = _minimal_linear_ts_mesh()
    td = tempfile.mkdtemp()
    jrd = os.path.join(td, "j")
    os.makedirs(os.path.join(jrd, "element_stiffness_matrices"), exist_ok=True)
    os.makedirs(os.path.join(jrd, "element_force_vectors"), exist_ok=True)
    try:
        ges = GeometricallyExactShearDeformableBeam3D(
            element_id=0,
            element_dictionary={**ed, "types": np.array(["GeometricallyExactShearDeformableBeam3D"])},
            grid_dictionary=gd,
            material_dictionary=md,
            section_dictionary=sd,
            point_load_array=np.empty((0, 9)),
            distributed_load_array=np.empty((0, 9)),
            job_results_dir=jrd,
        )
        nl = NonlinearTimoshenkoBeamElement3D(
            element_id=0,
            element_dictionary=ed,
            grid_dictionary=gd,
            material_dictionary=md,
            section_dictionary=sd,
            point_load_array=np.empty((0, 9)),
            distributed_load_array=np.empty((0, 9)),
            job_results_dir=jrd,
        )
        np.random.seed(0)
        u = np.random.randn(12) * 1e-5
        xi_g = 0.33
        _, _, _, _, E_nl, _ = nl._gp_tl_kinematics(u, xi_g)
        E_g = ges._gesdb_strain_voigt_at_gauss(u, xi_g)
        np.testing.assert_allclose(E_g, np.asarray(E_nl).ravel(), rtol=1e-12)
    finally:
        shutil.rmtree(td, ignore_errors=True)


def test_gesdb_small_displacement_strain_near_linear_B_u():
    """Green–Lagrange strain ≈ linearised B u for ||u|| → 0."""
    import tempfile
    import os
    import shutil
    from pre_processing.element_library.linear.beam.first_order_shear_deformation_theory.timoshenko.linear_timoshenko_3D import (
        LinearTimoshenkoBeamElement3D,
    )
    from pre_processing.element_library.nonlinear.large_rotations.geometrically_exact_shear_deformable_beam.geometrically_exact_shear_deformable_beam_3D import (
        GeometricallyExactShearDeformableBeam3D,
    )

    gd, ed, md, sd = _minimal_linear_ts_mesh()
    td = tempfile.mkdtemp()
    jrd = os.path.join(td, "j")
    os.makedirs(os.path.join(jrd, "element_stiffness_matrices"), exist_ok=True)
    os.makedirs(os.path.join(jrd, "element_force_vectors"), exist_ok=True)
    try:
        lin = LinearTimoshenkoBeamElement3D(
            element_id=0,
            element_dictionary=ed,
            grid_dictionary=gd,
            material_dictionary=md,
            section_dictionary=sd,
            point_load_array=np.empty((0, 9)),
            distributed_load_array=np.empty((0, 9)),
            job_results_dir=jrd,
        )
        ges = GeometricallyExactShearDeformableBeam3D(
            element_id=0,
            element_dictionary={**ed, "types": np.array(["GeometricallyExactShearDeformableBeam3D"])},
            grid_dictionary=gd,
            material_dictionary=md,
            section_dictionary=sd,
            point_load_array=np.empty((0, 9)),
            distributed_load_array=np.empty((0, 9)),
            job_results_dir=jrd,
        )
        np.random.seed(2)
        u = np.random.randn(12) * 1e-9
        xi_g = 0.0
        N, dN_dξ, d2N_dξ2 = lin.shape_function_operator.natural_coordinate_form(np.array([xi_g]))
        B = lin.strain_displacement_operator.physical_coordinate_form(dN_dξ, d2N_dξ2, N)[0]
        eps_lin = B @ u
        E_tl = ges._gesdb_strain_voigt_at_gauss(u, xi_g)
        np.testing.assert_allclose(E_tl, eps_lin, rtol=5e-5, atol=1e-15)
    finally:
        shutil.rmtree(td, ignore_errors=True)


def test_gesdb_internal_force_matches_parent_delegate():
    """Subclass shares TL delegate — smoke equality."""
    import tempfile
    import os
    import shutil
    from pre_processing.element_library.nonlinear.large_rotations.geometrically_exact_shear_deformable_beam.geometrically_exact_shear_deformable_beam_3D import (
        GeometricallyExactShearDeformableBeam3D,
    )
    from pre_processing.element_library.nonlinear.timoshenko.nonlinear_timoshenko_3D import (
        NonlinearTimoshenkoBeamElement3D,
    )

    gd, ed, md, sd = _minimal_linear_ts_mesh()
    td = tempfile.mkdtemp()
    jrd = os.path.join(td, "j")
    os.makedirs(os.path.join(jrd, "element_stiffness_matrices"), exist_ok=True)
    os.makedirs(os.path.join(jrd, "element_force_vectors"), exist_ok=True)
    try:
        ges = GeometricallyExactShearDeformableBeam3D(
            element_id=0,
            element_dictionary={**ed, "types": np.array(["GeometricallyExactShearDeformableBeam3D"])},
            grid_dictionary=gd,
            material_dictionary=md,
            section_dictionary=sd,
            point_load_array=np.empty((0, 9)),
            distributed_load_array=np.empty((0, 9)),
            job_results_dir=jrd,
        )
        np.random.seed(7)
        u = np.random.randn(12) * 1e-4
        F1 = ges.internal_force_vector(u)
        F2 = NonlinearTimoshenkoBeamElement3D(
            element_id=0,
            element_dictionary=ed,
            grid_dictionary=gd,
            material_dictionary=md,
            section_dictionary=sd,
            point_load_array=np.empty((0, 9)),
            distributed_load_array=np.empty((0, 9)),
            job_results_dir=jrd,
        ).internal_force_vector(u)
        np.testing.assert_allclose(np.asarray(F1).ravel(), np.asarray(F2).ravel(), rtol=1e-12)
    finally:
        shutil.rmtree(td, ignore_errors=True)


def test_gesdb_tangent_matches_parent_delegate():
    import tempfile
    import os
    import shutil
    from pre_processing.element_library.nonlinear.large_rotations.geometrically_exact_shear_deformable_beam.geometrically_exact_shear_deformable_beam_3D import (
        GeometricallyExactShearDeformableBeam3D,
    )
    from pre_processing.element_library.nonlinear.timoshenko.nonlinear_timoshenko_3D import (
        NonlinearTimoshenkoBeamElement3D,
    )

    gd, ed, md, sd = _minimal_linear_ts_mesh()
    td = tempfile.mkdtemp()
    jrd = os.path.join(td, "j")
    os.makedirs(os.path.join(jrd, "element_stiffness_matrices"), exist_ok=True)
    os.makedirs(os.path.join(jrd, "element_force_vectors"), exist_ok=True)
    try:
        ges = GeometricallyExactShearDeformableBeam3D(
            element_id=0,
            element_dictionary={**ed, "types": np.array(["GeometricallyExactShearDeformableBeam3D"])},
            grid_dictionary=gd,
            material_dictionary=md,
            section_dictionary=sd,
            point_load_array=np.empty((0, 9)),
            distributed_load_array=np.empty((0, 9)),
            job_results_dir=jrd,
        )
        u = np.random.randn(12) * 1e-5
        K1 = ges.tangent_stiffness_matrix(u)
        K2 = NonlinearTimoshenkoBeamElement3D(
            element_id=0,
            element_dictionary=ed,
            grid_dictionary=gd,
            material_dictionary=md,
            section_dictionary=sd,
            point_load_array=np.empty((0, 9)),
            distributed_load_array=np.empty((0, 9)),
            job_results_dir=jrd,
        ).tangent_stiffness_matrix(u)
        np.testing.assert_allclose(K1, K2, rtol=1e-12)
    finally:
        shutil.rmtree(td, ignore_errors=True)


def test_gesdb_native_kernel_small_strain_near_tl():
    """Native engineering axial + linear rows ~ TL for ‖u‖ → 0."""
    import tempfile
    import os
    import shutil
    from pre_processing.element_library.nonlinear.large_rotations.geometrically_exact_shear_deformable_beam.geometrically_exact_shear_deformable_beam_3D import (
        GeometricallyExactShearDeformableBeam3D,
    )
    from pre_processing.element_library.nonlinear.timoshenko.nonlinear_timoshenko_3D import (
        NonlinearTimoshenkoBeamElement3D,
    )

    gd, ed, md, sd = _minimal_linear_ts_mesh()
    td = tempfile.mkdtemp()
    jrd = os.path.join(td, "j")
    os.makedirs(os.path.join(jrd, "element_stiffness_matrices"), exist_ok=True)
    os.makedirs(os.path.join(jrd, "element_force_vectors"), exist_ok=True)
    try:
        ges = GeometricallyExactShearDeformableBeam3D(
            element_id=0,
            element_dictionary={**ed, "types": np.array(["GeometricallyExactShearDeformableBeam3D"])},
            grid_dictionary=gd,
            material_dictionary=md,
            section_dictionary=sd,
            point_load_array=np.empty((0, 9)),
            distributed_load_array=np.empty((0, 9)),
            job_results_dir=jrd,
            gesdb_kernel="native",
        )
        nl = NonlinearTimoshenkoBeamElement3D(
            element_id=0,
            element_dictionary=ed,
            grid_dictionary=gd,
            material_dictionary=md,
            section_dictionary=sd,
            point_load_array=np.empty((0, 9)),
            distributed_load_array=np.empty((0, 9)),
            job_results_dir=jrd,
        )
        np.random.seed(11)
        u = np.random.randn(12) * 1e-8
        xi_g = 0.25
        _, _, _, _, E_tl, _ = nl._gp_tl_kinematics(u, xi_g)
        E_nt = ges._tl_voigt_strain_at_gauss(u, xi_g)
        np.testing.assert_allclose(np.asarray(E_nt).ravel(), np.asarray(E_tl).ravel(), rtol=1e-5, atol=1e-18)
    finally:
        shutil.rmtree(td, ignore_errors=True)


def test_gesdb_native_vs_tl_internal_force_moderate_u():
    """Native kernel differs from TL at finite ‖u‖ (documented regression band)."""
    import tempfile
    import os
    import shutil
    from pre_processing.element_library.nonlinear.large_rotations.geometrically_exact_shear_deformable_beam.geometrically_exact_shear_deformable_beam_3D import (
        GeometricallyExactShearDeformableBeam3D,
    )
    from pre_processing.element_library.nonlinear.timoshenko.nonlinear_timoshenko_3D import (
        NonlinearTimoshenkoBeamElement3D,
    )

    gd, ed, md, sd = _minimal_linear_ts_mesh()
    td = tempfile.mkdtemp()
    jrd = os.path.join(td, "j")
    os.makedirs(os.path.join(jrd, "element_stiffness_matrices"), exist_ok=True)
    os.makedirs(os.path.join(jrd, "element_force_vectors"), exist_ok=True)
    try:
        ges = GeometricallyExactShearDeformableBeam3D(
            element_id=0,
            element_dictionary={**ed, "types": np.array(["GeometricallyExactShearDeformableBeam3D"])},
            grid_dictionary=gd,
            material_dictionary=md,
            section_dictionary=sd,
            point_load_array=np.empty((0, 9)),
            distributed_load_array=np.empty((0, 9)),
            job_results_dir=jrd,
            gesdb_kernel="native",
        )
        nl = NonlinearTimoshenkoBeamElement3D(
            element_id=0,
            element_dictionary=ed,
            grid_dictionary=gd,
            material_dictionary=md,
            section_dictionary=sd,
            point_load_array=np.empty((0, 9)),
            distributed_load_array=np.empty((0, 9)),
            job_results_dir=jrd,
        )
        np.random.seed(13)
        u = np.random.randn(12) * 0.02
        f_n = ges.internal_force_vector(u)
        f_t = nl.internal_force_vector(u)
        assert np.linalg.norm(np.asarray(f_n) - np.asarray(f_t)) > 1e-9 * max(
            1.0, np.linalg.norm(f_t)
        )
    finally:
        shutil.rmtree(td, ignore_errors=True)

"""Regression checks for co-rotational beam kinematics and large-rotation registrations."""

from __future__ import annotations

import os
import shutil
import sys
import tempfile
from pathlib import Path

import numpy as np
import pytest

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))


def _job_dirs(base: str) -> str:
    os.makedirs(os.path.join(base, "element_stiffness_matrices"), exist_ok=True)
    os.makedirs(os.path.join(base, "element_force_vectors"), exist_ok=True)
    os.makedirs(os.path.join(base, "logs"), exist_ok=True)
    return base


@pytest.fixture
def minimal_mesh_dicts():
    grid_dictionary = {
        "ids": np.array([0, 1]),
        "coordinates": np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]], dtype=np.float64),
    }
    element_dictionary = {
        "ids": np.array([0]),
        "connectivity": np.array([[0, 1]]),
        "types": np.array(["CorotationalBeamElement3D"]),
        "integration_orders": {
            "axial": np.array([2]),
            "bending_y": np.array([2]),
            "bending_z": np.array([2]),
            "shear_y": np.array([2]),
            "shear_z": np.array([2]),
            "torsion": np.array([2]),
            "load": np.array([2]),
        },
    }
    material_dictionary = {
        "E": np.array([200e9]),
        "G": np.array([77e9]),
        "nu": np.array([0.3]),
        "rho": np.array([7850.0]),
    }
    section_dictionary = {
        "A": np.array([0.01]),
        "I_x": np.array([0.0]),
        "I_y": np.array([1e-6]),
        "I_z": np.array([1e-6]),
        "J_t": np.array([2e-7]),
    }
    return grid_dictionary, element_dictionary, material_dictionary, section_dictionary


def test_corotational_internal_force_zero_at_zero(minimal_mesh_dicts):
    from pre_processing.element_library.nonlinear.large_rotations.corotational.corotational_3D import (
        CorotationalBeamElement3D,
    )

    grid_dictionary, element_dictionary, material_dictionary, section_dictionary = minimal_mesh_dicts
    tmp = tempfile.mkdtemp()
    job_root = _job_dirs(os.path.join(tmp, "job"))
    try:
        el = CorotationalBeamElement3D(
            element_id=0,
            element_dictionary=element_dictionary,
            grid_dictionary=grid_dictionary,
            section_dictionary=section_dictionary,
            material_dictionary=material_dictionary,
            point_load_array=np.empty((0, 9)),
            distributed_load_array=np.empty((0, 9)),
            job_results_dir=job_root,
            kernel="timoshenko",
        )
        u = np.zeros(12, dtype=np.float64)
        f = el.internal_force_vector(u)
        assert np.allclose(f, 0.0, atol=1e-9)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def test_corotational_tangent_matches_stiffness_at_zero(minimal_mesh_dicts):
    from pre_processing.element_library.nonlinear.large_rotations.corotational.corotational_3D import (
        CorotationalBeamElement3D,
    )

    grid_dictionary, element_dictionary, material_dictionary, section_dictionary = minimal_mesh_dicts
    tmp = tempfile.mkdtemp()
    job_root = _job_dirs(os.path.join(tmp, "job"))
    try:
        el = CorotationalBeamElement3D(
            element_id=0,
            element_dictionary=element_dictionary,
            grid_dictionary=grid_dictionary,
            section_dictionary=section_dictionary,
            material_dictionary=material_dictionary,
            point_load_array=np.empty((0, 9)),
            distributed_load_array=np.empty((0, 9)),
            job_results_dir=job_root,
            kernel="timoshenko",
        )
        u = np.zeros(12, dtype=np.float64)
        Ke = el.element_stiffness_matrix().K_e
        Kt = el.tangent_stiffness_matrix(u)
        assert np.allclose(Kt, Ke, rtol=5e-5, atol=5e-5)
        err_sym = np.linalg.norm(Kt - Kt.T)
        assert err_sym < 1e-8
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def test_corotational_elastic_material_tangent_matches_ke_at_zero(minimal_mesh_dicts):
    from pre_processing.element_library.nonlinear.large_rotations.corotational.corotational_3D import (
        CorotationalBeamElement3D,
    )

    grid_dictionary, element_dictionary, material_dictionary, section_dictionary = minimal_mesh_dicts
    tmp = tempfile.mkdtemp()
    job_root = _job_dirs(os.path.join(tmp, "job"))
    try:
        el = CorotationalBeamElement3D(
            element_id=0,
            element_dictionary=element_dictionary,
            grid_dictionary=grid_dictionary,
            section_dictionary=section_dictionary,
            material_dictionary=material_dictionary,
            point_load_array=np.empty((0, 9)),
            distributed_load_array=np.empty((0, 9)),
            job_results_dir=job_root,
            kernel="timoshenko",
            tangent_stiffness_mode="elastic_material",
        )
        u = np.zeros(12, dtype=np.float64)
        Ke = el.element_stiffness_matrix().K_e
        Kt = el.tangent_stiffness_matrix(u)
        assert np.allclose(Kt, Ke, rtol=1e-10, atol=1e-10)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def test_corotational_elastic_material_symmetric_at_small_u(minimal_mesh_dicts):
    from pre_processing.element_library.nonlinear.large_rotations.corotational.corotational_3D import (
        CorotationalBeamElement3D,
    )

    grid_dictionary, element_dictionary, material_dictionary, section_dictionary = minimal_mesh_dicts
    tmp = tempfile.mkdtemp()
    job_root = _job_dirs(os.path.join(tmp, "job"))
    try:
        el = CorotationalBeamElement3D(
            element_id=0,
            element_dictionary=element_dictionary,
            grid_dictionary=grid_dictionary,
            section_dictionary=section_dictionary,
            material_dictionary=material_dictionary,
            point_load_array=np.empty((0, 9)),
            distributed_load_array=np.empty((0, 9)),
            job_results_dir=job_root,
            kernel="timoshenko",
            tangent_stiffness_mode="elastic_material",
        )
        rng = np.random.default_rng(0)
        u = 1e-3 * rng.standard_normal(12)
        Kt = el.tangent_stiffness_matrix(u)
        assert np.linalg.norm(Kt - Kt.T) < 1e-10
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def test_corotational_fd_and_elastic_material_differ_with_transverse_translation(minimal_mesh_dicts):
    """Finite spin of chord: elastic-material block omits spin stiffness vs FD consistent tangent."""
    from pre_processing.element_library.nonlinear.large_rotations.corotational.corotational_3D import (
        CorotationalBeamElement3D,
    )

    grid_dictionary, element_dictionary, material_dictionary, section_dictionary = minimal_mesh_dicts
    tmp = tempfile.mkdtemp()
    job_root = _job_dirs(os.path.join(tmp, "job"))
    try:
        el_fd = CorotationalBeamElement3D(
            element_id=0,
            element_dictionary=element_dictionary,
            grid_dictionary=grid_dictionary,
            section_dictionary=section_dictionary,
            material_dictionary=material_dictionary,
            point_load_array=np.empty((0, 9)),
            distributed_load_array=np.empty((0, 9)),
            job_results_dir=job_root,
            kernel="timoshenko",
            tangent_stiffness_mode="finite_difference",
        )
        el_em = CorotationalBeamElement3D(
            element_id=0,
            element_dictionary=element_dictionary,
            grid_dictionary=grid_dictionary,
            section_dictionary=section_dictionary,
            material_dictionary=material_dictionary,
            point_load_array=np.empty((0, 9)),
            distributed_load_array=np.empty((0, 9)),
            job_results_dir=job_root,
            kernel="timoshenko",
            tangent_stiffness_mode="elastic_material",
        )
        u = np.zeros(12, dtype=np.float64)
        u[1] = 0.05
        u[7] = -0.05
        Kfd = el_fd.tangent_stiffness_matrix(u)
        Kem = el_em.tangent_stiffness_matrix(u)
        assert np.linalg.norm(Kfd - Kem) > 1e-6
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def test_gesdb_matches_tl_timoshenko_internal_force_small_u(minimal_mesh_dicts):
    from pre_processing.element_library.nonlinear.large_rotations.geometrically_exact_shear_deformable_beam.geometrically_exact_shear_deformable_beam_3D import (
        GeometricallyExactShearDeformableBeam3D,
    )
    from pre_processing.element_library.nonlinear.timoshenko.nonlinear_timoshenko_3D import (
        NonlinearTimoshenkoBeamElement3D,
    )

    grid_dictionary, element_dictionary, material_dictionary, section_dictionary = minimal_mesh_dicts
    element_dictionary = {**element_dictionary, "types": np.array(["GeometricallyExactShearDeformableBeam3D"])}
    tmp = tempfile.mkdtemp()
    job_root = _job_dirs(os.path.join(tmp, "job"))
    try:
        ges = GeometricallyExactShearDeformableBeam3D(
            element_id=0,
            element_dictionary=element_dictionary,
            grid_dictionary=grid_dictionary,
            section_dictionary=section_dictionary,
            material_dictionary=material_dictionary,
            point_load_array=np.empty((0, 9)),
            distributed_load_array=np.empty((0, 9)),
            job_results_dir=job_root,
        )
        tl_types = element_dictionary.copy()
        tl_types["types"] = np.array(["NonlinearTimoshenkoBeamElement3D"])
        tl = NonlinearTimoshenkoBeamElement3D(
            element_id=0,
            element_dictionary=tl_types,
            grid_dictionary=grid_dictionary,
            section_dictionary=section_dictionary,
            material_dictionary=material_dictionary,
            point_load_array=np.empty((0, 9)),
            distributed_load_array=np.empty((0, 9)),
            job_results_dir=job_root,
        )
        rng = np.random.default_rng(42)
        u = 1e-4 * rng.standard_normal(12)
        assert np.allclose(ges.internal_force_vector(u), tl.internal_force_vector(u), rtol=1e-12, atol=1e-12)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

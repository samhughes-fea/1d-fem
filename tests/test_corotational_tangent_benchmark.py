"""Benchmark: corotational FD vs elastic_material tangent matrices at finite displacement."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))


@pytest.fixture
def minimal_mesh_dicts():
    import tempfile
    import os

    td = tempfile.mkdtemp()
    jrd = os.path.join(td, "job")
    os.makedirs(os.path.join(jrd, "element_stiffness_matrices"), exist_ok=True)
    os.makedirs(os.path.join(jrd, "element_force_vectors"), exist_ok=True)
    grid_dictionary = {
        "ids": np.array([0, 1], dtype=np.int32),
        "coordinates": np.array([[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]], dtype=np.float64),
    }
    element_dictionary = {
        "ids": np.array([0], dtype=np.int32),
        "connectivity": np.array([[0, 1]], dtype=np.int32),
        "types": np.array(["CorotationalBeamElement3D"]),
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
        "E": np.array([200e9]),
        "G": np.array([77e9]),
        "nu": np.array([0.3]),
        "rho": np.array([7850.0]),
    }
    section_dictionary = {
        "A": np.array([0.02]),
        "I_x": np.array([1e-9]),
        "I_y": np.array([2e-6]),
        "I_z": np.array([2e-6]),
        "J_t": np.array([3e-9]),
        "kappa": np.array([5.0 / 6.0]),
        "alpha": np.array([0.0]),
    }
    yield grid_dictionary, element_dictionary, material_dictionary, section_dictionary, jrd, td
    import shutil

    shutil.rmtree(td, ignore_errors=True)


def test_corotational_fd_vs_elastic_material_frobenius_ratio_bounded(minimal_mesh_dicts):
    from pre_processing.element_library.nonlinear.large_rotations.corotational.corotational_3D import (
        CorotationalBeamElement3D,
    )

    gd, ed, md, sd, jrd, _ = minimal_mesh_dicts
    el_fd = CorotationalBeamElement3D(
        element_id=0,
        element_dictionary=ed,
        grid_dictionary=gd,
        material_dictionary=md,
        section_dictionary=sd,
        point_load_array=np.empty((0, 9)),
        distributed_load_array=np.empty((0, 9)),
        job_results_dir=jrd,
        tangent_stiffness_mode="finite_difference",
    )
    el_em = CorotationalBeamElement3D(
        element_id=0,
        element_dictionary=ed,
        grid_dictionary=gd,
        material_dictionary=md,
        section_dictionary=sd,
        point_load_array=np.empty((0, 9)),
        distributed_load_array=np.empty((0, 9)),
        job_results_dir=jrd,
        tangent_stiffness_mode="elastic_material",
    )
    np.random.seed(42)
    u = np.zeros(12, dtype=np.float64)
    u[7] = 0.05
    u[2] = 0.02
    Kfd = el_fd.tangent_stiffness_matrix(u)
    Kem = el_em.tangent_stiffness_matrix(u)
    nf = np.linalg.norm(Kfd, ord="fro")
    ne = np.linalg.norm(Kem, ord="fro")
    assert nf > 0 and ne > 0
    # Modes differ; both finite symmetric-ish stabilisations — ratio sanity check only.
    ratio = nf / ne
    assert 0.3 < ratio < 3.0

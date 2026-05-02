"""ElementFactory passes simulation_settings nonlinear corotational_tangent_mode."""

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


def _mesh():
    grid_dictionary = {
        "ids": np.array([0, 1], dtype=np.int32),
        "coordinates": np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]], dtype=np.float64),
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
        "A": np.array([0.01]),
        "I_x": np.array([1e-9]),
        "I_y": np.array([1e-6]),
        "I_z": np.array([1e-6]),
        "J_t": np.array([2e-9]),
        "kappa": np.array([5.0 / 6.0]),
        "alpha": np.array([0.0]),
    }
    return grid_dictionary, element_dictionary, material_dictionary, section_dictionary


def test_factory_sets_elastic_material_from_simulation_settings():
    from pre_processing.element_library.element_factory import ElementFactory

    gd, ed, md, sd = _mesh()
    td = tempfile.mkdtemp()
    jrd = os.path.join(td, "j")
    os.makedirs(os.path.join(jrd, "element_stiffness_matrices"), exist_ok=True)
    os.makedirs(os.path.join(jrd, "element_force_vectors"), exist_ok=True)
    os.makedirs(os.path.join(jrd, "logs"), exist_ok=True)
    try:
        factory = ElementFactory(job_results_dir=jrd)
        elements = factory.create_elements_batch(
            element_ids=np.array([0]),
            element_dictionary=ed,
            grid_dictionary=gd,
            material_dictionary=md,
            section_dictionary=sd,
            point_load_array=np.empty((0, 9)),
            distributed_load_array=np.empty((0, 9)),
            simulation_settings={"nonlinear": {"corotational_tangent_mode": "elastic_material"}},
        )
        assert elements[0].tangent_stiffness_mode == "elastic_material"
    finally:
        shutil.rmtree(td, ignore_errors=True)


def test_factory_invalid_corotational_mode_raises():
    from pre_processing.element_library.element_factory import ElementFactory

    gd, ed, md, sd = _mesh()
    td = tempfile.mkdtemp()
    jrd = os.path.join(td, "j")
    os.makedirs(os.path.join(jrd, "element_stiffness_matrices"), exist_ok=True)
    os.makedirs(os.path.join(jrd, "element_force_vectors"), exist_ok=True)
    os.makedirs(os.path.join(jrd, "logs"), exist_ok=True)
    try:
        factory = ElementFactory(job_results_dir=jrd)
        with pytest.raises(ValueError, match="corotational_tangent_mode"):
            factory.create_elements_batch(
                element_ids=np.array([0]),
                element_dictionary=ed,
                grid_dictionary=gd,
                material_dictionary=md,
                section_dictionary=sd,
                point_load_array=np.empty((0, 9)),
                distributed_load_array=np.empty((0, 9)),
                simulation_settings={"nonlinear": {"corotational_tangent_mode": "bogus"}},
            )
    finally:
        shutil.rmtree(td, ignore_errors=True)


def test_factory_sets_gesdb_tl_fallback_from_simulation_settings():
    from pre_processing.element_library.element_factory import ElementFactory

    grid_dictionary = {
        "ids": np.array([0, 1], dtype=np.int32),
        "coordinates": np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]], dtype=np.float64),
    }
    element_dictionary = {
        "ids": np.array([0], dtype=np.int32),
        "connectivity": np.array([[0, 1]], dtype=np.int32),
        "types": np.array(["GeometricallyExactShearDeformableBeam3D"]),
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
        "A": np.array([0.01]),
        "I_x": np.array([1e-9]),
        "I_y": np.array([1e-6]),
        "I_z": np.array([1e-6]),
        "J_t": np.array([2e-9]),
        "kappa": np.array([5.0 / 6.0]),
        "alpha": np.array([0.0]),
    }
    td = tempfile.mkdtemp()
    jrd = os.path.join(td, "j")
    os.makedirs(os.path.join(jrd, "element_stiffness_matrices"), exist_ok=True)
    os.makedirs(os.path.join(jrd, "element_force_vectors"), exist_ok=True)
    os.makedirs(os.path.join(jrd, "logs"), exist_ok=True)
    try:
        factory = ElementFactory(job_results_dir=jrd)
        elements = factory.create_elements_batch(
            element_ids=np.array([0]),
            element_dictionary=element_dictionary,
            grid_dictionary=grid_dictionary,
            material_dictionary=material_dictionary,
            section_dictionary=section_dictionary,
            point_load_array=np.empty((0, 9)),
            distributed_load_array=np.empty((0, 9)),
            simulation_settings={"nonlinear": {"gesdb_tl_fallback": True}},
        )
        assert elements[0]._gesdb_tl_fallback is True
    finally:
        shutil.rmtree(td, ignore_errors=True)


def test_factory_sets_gesdb_kernel_native():
    from pre_processing.element_library.element_factory import ElementFactory

    grid_dictionary = {
        "ids": np.array([0, 1], dtype=np.int32),
        "coordinates": np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]], dtype=np.float64),
    }
    element_dictionary = {
        "ids": np.array([0], dtype=np.int32),
        "connectivity": np.array([[0, 1]], dtype=np.int32),
        "types": np.array(["GeometricallyExactShearDeformableBeam3D"]),
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
        "A": np.array([0.01]),
        "I_x": np.array([1e-9]),
        "I_y": np.array([1e-6]),
        "I_z": np.array([1e-6]),
        "J_t": np.array([2e-9]),
        "kappa": np.array([5.0 / 6.0]),
        "alpha": np.array([0.0]),
    }
    td = tempfile.mkdtemp()
    jrd = os.path.join(td, "j")
    os.makedirs(os.path.join(jrd, "element_stiffness_matrices"), exist_ok=True)
    os.makedirs(os.path.join(jrd, "element_force_vectors"), exist_ok=True)
    os.makedirs(os.path.join(jrd, "logs"), exist_ok=True)
    try:
        factory = ElementFactory(job_results_dir=jrd)
        elements = factory.create_elements_batch(
            element_ids=np.array([0]),
            element_dictionary=element_dictionary,
            grid_dictionary=grid_dictionary,
            material_dictionary=material_dictionary,
            section_dictionary=section_dictionary,
            point_load_array=np.empty((0, 9)),
            distributed_load_array=np.empty((0, 9)),
            simulation_settings={"nonlinear": {"gesdb_kernel": "native"}},
        )
        assert elements[0]._gesdb_kernel == "native"
    finally:
        shutil.rmtree(td, ignore_errors=True)

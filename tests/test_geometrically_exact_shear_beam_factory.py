"""Factory smoke test for classical GEBT stub registration."""

import os
import shutil
import sys
import tempfile
from pathlib import Path

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))


def test_geometrically_exact_shear_beam_factory_instantiation():
    from pre_processing.element_library.element_factory import ElementFactory
    from pre_processing.element_library.nonlinear.geometrically_exact_shear_deformable_beam.geometrically_exact_shear_deformable_beam_3D import (
        GeometricallyExactShearDeformableBeam3D,
    )

    grid_dictionary = {"ids": np.array([0, 1]), "coordinates": np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])}
    element_dictionary = {
        "ids": np.array([0]),
        "connectivity": np.array([[0, 1]]),
        "types": np.array(["GeometricallyExactShearDeformableBeam3D"]),
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
    material_dictionary = {"E": np.array([200e9]), "G": np.array([77e9]), "nu": np.array([0.3]), "rho": np.array([7850.0])}
    section_dictionary = {
        "A": np.array([0.01]),
        "I_x": np.array([0.0]),
        "I_y": np.array([1e-6]),
        "I_z": np.array([1e-6]),
        "J_t": np.array([2e-7]),
    }
    temp_dir = tempfile.mkdtemp()
    job_results_dir = os.path.join(temp_dir, "ge_exact_stub")
    os.makedirs(job_results_dir, exist_ok=True)
    os.makedirs(os.path.join(job_results_dir, "element_stiffness_matrices"), exist_ok=True)
    os.makedirs(os.path.join(job_results_dir, "element_force_vectors"), exist_ok=True)
    os.makedirs(os.path.join(job_results_dir, "logs"), exist_ok=True)
    try:
        factory = ElementFactory(job_results_dir=job_results_dir)
        elements = factory.create_elements_batch(
            element_ids=np.array([0]),
            element_dictionary=element_dictionary,
            grid_dictionary=grid_dictionary,
            material_dictionary=material_dictionary,
            section_dictionary=section_dictionary,
            point_load_array=np.empty((0, 9)),
            distributed_load_array=np.empty((0, 9)),
        )
        assert len(elements) == 1
        assert isinstance(elements[0], GeometricallyExactShearDeformableBeam3D)
    finally:
        shutil.rmtree(temp_dir, ignore_errors=True)

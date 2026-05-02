"""Removed legacy element type strings: ValueError + DeprecationWarning from factory guard."""

import os
import sys
import tempfile
from pathlib import Path

import numpy as np
import pytest

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from pre_processing.element_library.element_factory import ElementFactory
from pre_processing.element_library.removed_element_types import (
    REMOVED_ELEMENT_TYPES,
    ensure_element_type_allowed,
)


@pytest.mark.parametrize(
    "removed_type",
    sorted(REMOVED_ELEMENT_TYPES.keys()),
)
def test_ensure_element_type_allowed_warns_and_raises(removed_type):
    with pytest.warns(DeprecationWarning):
        with pytest.raises(ValueError, match="was removed|not registered"):
            ensure_element_type_allowed(removed_type)


def test_factory_rejects_removed_linear_warping_timoshenko():
    tmp = tempfile.mkdtemp()
    job_results_dir = os.path.join(tmp, "job")
    os.makedirs(os.path.join(job_results_dir, "logs"), exist_ok=True)
    try:
        factory = ElementFactory(job_results_dir=job_results_dir)
        element_dictionary = {
            "ids": np.array([0]),
            "connectivity": np.array([[0, 1]]),
            "types": np.array(["LinearWarpingTimoshenkoBeamElement3D"]),
            "warping": np.array([1], dtype=np.int8),
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
        grid_dictionary = {"coordinates": np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])}
        material_dictionary = {
            "E": np.array([2.1e11]),
            "G": np.array([8.1e10]),
            "nu": np.array([0.3]),
            "rho": np.array([7850.0]),
        }
        section_dictionary = {
            "A": np.array([0.001]),
            "I_x": np.array([0.0]),
            "I_y": np.array([1e-8]),
            "I_z": np.array([1e-8]),
            "J_t": np.array([1e-9]),
            "kappa": np.array([5.0 / 6.0]),
            "alpha": np.array([0.0]),
            "y_sc": np.array([0.0]),
            "z_sc": np.array([0.0]),
            "Gamma": np.array([1e-8]),
        }
        pl = np.empty((0, 9))
        dl = np.empty((0, 9))
        with pytest.warns(DeprecationWarning):
            with pytest.raises(ValueError, match="was removed"):
                factory.create_elements_batch(
                    element_ids=np.array([0]),
                    element_dictionary=element_dictionary,
                    grid_dictionary=grid_dictionary,
                    material_dictionary=material_dictionary,
                    section_dictionary=section_dictionary,
                    point_load_array=pl,
                    distributed_load_array=dl,
                )
    finally:
        import shutil

        shutil.rmtree(tmp, ignore_errors=True)

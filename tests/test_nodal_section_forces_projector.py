"""
Unit tests for NodalSectionForcesProjector.

Asserts that constant section forces at Gauss points yield nodal values
equal to that constant (no boundary spike), and that boundary nodes
receive the element mean.
"""

import os
import sys
import tempfile
from pathlib import Path

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))


def _minimal_formulation_cache_and_mesh():
    """One Euler-Bernoulli element, formulation cache with shape_functions."""
    from pre_processing.element_library.euler_bernoulli.euler_bernoulli_3D import (
        EulerBernoulliBeamElement3D,
    )
    from processing_OOP.static.results.containers import FormulationResultSet

    L = 1.0
    element_dictionary = {
        "ids": np.array([0]),
        "connectivity": np.array([[0, 1]]),
        "types": np.array(["EulerBernoulliBeamElement3D"]),
        "integration_orders": {
            "axial": np.array([3]),
            "bending_y": np.array([3]),
            "bending_z": np.array([3]),
            "shear_y": np.array([2]),
            "shear_z": np.array([2]),
            "torsion": np.array([3]),
            "load": np.array([2]),
        },
    }
    grid_dictionary = {
        "coordinates": np.array([[0.0, 0.0, 0.0], [L, 0.0, 0.0]]),
    }
    material_dictionary = {
        "E": np.array([2.1e11]),
        "G": np.array([8.1e10]),
        "nu": np.array([0.3]),
        "rho": np.array([7850.0]),
    }
    section_dictionary = {
        "A": np.array([1.0e-3]),
        "I_x": np.array([1.0e-9]),
        "I_y": np.array([1.0e-8]),
        "I_z": np.array([1.0e-8]),
        "J_t": np.array([1.0e-9]),
    }
    point_load_array = np.empty((0, 9))
    distributed_load_array = np.array([
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [L, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    ])
    temp_dir = tempfile.mkdtemp()
    job_results_dir = os.path.join(temp_dir, "test_results")
    os.makedirs(job_results_dir, exist_ok=True)
    os.makedirs(os.path.join(job_results_dir, "element_stiffness_matrices"), exist_ok=True)
    os.makedirs(os.path.join(job_results_dir, "element_force_vectors"), exist_ok=True)
    element = EulerBernoulliBeamElement3D(
        element_id=0,
        element_dictionary=element_dictionary,
        grid_dictionary=grid_dictionary,
        material_dictionary=material_dictionary,
        section_dictionary=section_dictionary,
        point_load_array=point_load_array,
        distributed_load_array=distributed_load_array,
        job_results_dir=job_results_dir,
    )
    elem_obj = element.element_stiffness_matrix()
    formulation_cache = FormulationResultSet(
        element_objects=[elem_obj],
        force_objects=[element.element_force_vector()],
    )
    return formulation_cache, element_dictionary, grid_dictionary


def test_nodal_section_forces_projector_constant_vy():
    """Constant Vy at all GPs yields nodal Vy equal to that value (no boundary spike)."""
    from processing_OOP.static.results.compute_tertiary.nodal_section_forces_projector import (
        NodalSectionForcesProjector,
    )

    formulation_cache, element_dictionary, grid_dictionary = _minimal_formulation_cache_and_mesh()
    elem_obj = formulation_cache.element_objects[0]
    n_gp = len(elem_obj.gauss_data)

    # Section forces: [N, Vy, Vz, T, My, Mz]; constant Vy = -500 N
    constant_vy = -500.0
    elem_section_forces = [
        np.array([0.0, constant_vy, 0.0, 0.0, 0.0, 0.0]) for _ in range(n_gp)
    ]
    section_forces_gauss = [elem_section_forces]

    projector = NodalSectionForcesProjector(
        section_forces_gauss=section_forces_gauss,
        formulation_cache=formulation_cache,
        element_dictionary=element_dictionary,
        grid_dictionary=grid_dictionary,
    )
    nodal = projector.project()

    assert nodal.shape == (2, 6)
    # Vy is column index 1
    np.testing.assert_allclose(nodal[:, 1], [constant_vy, constant_vy], rtol=1e-5)
    assert np.all(np.isfinite(nodal))


def test_nodal_section_forces_projector_constant_all_components():
    """Constant section forces (all components) at GPs yield same nodal values."""
    from processing_OOP.static.results.compute_tertiary.nodal_section_forces_projector import (
        NodalSectionForcesProjector,
    )

    formulation_cache, element_dictionary, grid_dictionary = _minimal_formulation_cache_and_mesh()
    elem_obj = formulation_cache.element_objects[0]
    n_gp = len(elem_obj.gauss_data)

    # All components constant: e.g. [100, -500, -300, 10, 20, 30]
    constant = np.array([100.0, -500.0, -300.0, 10.0, 20.0, 30.0])
    elem_section_forces = [constant.copy() for _ in range(n_gp)]
    section_forces_gauss = [elem_section_forces]

    projector = NodalSectionForcesProjector(
        section_forces_gauss=section_forces_gauss,
        formulation_cache=formulation_cache,
        element_dictionary=element_dictionary,
        grid_dictionary=grid_dictionary,
    )
    nodal = projector.project()

    assert nodal.shape == (2, 6)
    np.testing.assert_allclose(nodal[0], constant, rtol=1e-5)
    np.testing.assert_allclose(nodal[1], constant, rtol=1e-5)

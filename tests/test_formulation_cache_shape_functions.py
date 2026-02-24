"""
Unit tests for formulation cache shape functions contract.

Asserts that ElementObject and ForceObject produced by the element library
have shape_functions (and shape_derivatives for stiffness) populated at
every Gauss point, as required by the results pipeline.
"""

import os
import sys
import tempfile
from pathlib import Path

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from pre_processing.element_library.euler_bernoulli.euler_bernoulli_3D import (
    EulerBernoulliBeamElement3D,
)


def _minimal_eb_element_and_job_dir():
    """Create minimal Euler-Bernoulli element and job_results_dir with required subdirs."""
    L = 1.0
    E = 2.1e11
    G = 8.1e10
    A = 1.0e-3
    I_y = 1.0e-8
    I_z = 1.0e-8
    J_t = 1.0e-9
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
        "E": np.array([E]),
        "G": np.array([G]),
        "nu": np.array([0.3]),
        "rho": np.array([7850.0]),
    }
    section_dictionary = {
        "A": np.array([A]),
        "I_x": np.array([J_t]),
        "I_y": np.array([I_y]),
        "I_z": np.array([I_z]),
        "J_t": np.array([J_t]),
    }
    point_load_array = np.empty((0, 9))
    # Two x-positions so interpolator has a valid range; columns [x, y, z, Fx, Fy, Fz, Mx, My, Mz]
    distributed_load_array = np.array([
        [0.0, 0.0, 0.0, 0.0, -10.0, 0.0, 0.0, 0.0, 0.0],
        [L, 0.0, 0.0, 0.0, -10.0, 0.0, 0.0, 0.0, 0.0],
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
    return element, job_results_dir


def test_euler_bernoulli_element_object_shape_functions_populated():
    """ElementObject from Euler-Bernoulli 3D has shape_functions and shape_derivatives at every Gauss point."""
    element, job_results_dir = _minimal_eb_element_and_job_dir()
    try:
        obj = element.element_stiffness_matrix()
        assert obj.gauss_data, "expected at least one Gauss point"
        for i, gp in enumerate(obj.gauss_data):
            assert gp.shape_functions is not None, (
                f"ElementObject element_id={obj.element_id} gauss_data[{i}]: shape_functions must be set"
            )
            assert gp.shape_derivatives is not None, (
                f"ElementObject element_id={obj.element_id} gauss_data[{i}]: shape_derivatives must be set"
            )
    finally:
        import shutil
        if os.path.isdir(job_results_dir):
            shutil.rmtree(os.path.dirname(job_results_dir), ignore_errors=True)


def test_euler_bernoulli_element_object_b1_evaluate_shape_functions():
    """ElementObject from Euler-Bernoulli 3D has evaluate_shape_functions (B1) and it matches operator."""
    element, job_results_dir = _minimal_eb_element_and_job_dir()
    try:
        obj = element.element_stiffness_matrix()
        assert obj.evaluate_shape_functions is not None, "Euler-Bernoulli ElementObject must set evaluate_shape_functions (B1)"
        xi_test = np.array([-0.5, 0.0, 0.5])
        N_cache, dN_cache, d2N_cache = obj.evaluate_shape_functions(xi_test)
        N_op, dN_op, d2N_op = element.shape_function_operator.natural_coordinate_form(xi_test)
        np.testing.assert_allclose(N_cache, N_op, err_msg="B1 callable N(ξ) should match operator")
        np.testing.assert_allclose(dN_cache, dN_op, err_msg="B1 callable dN/dξ should match operator")
        np.testing.assert_allclose(d2N_cache, d2N_op, err_msg="B1 callable d2N/dξ2 should match operator")
    finally:
        import shutil
        if os.path.isdir(job_results_dir):
            shutil.rmtree(os.path.dirname(job_results_dir), ignore_errors=True)


def test_euler_bernoulli_force_object_shape_functions_populated():
    """ForceObject from Euler-Bernoulli 3D has shape_functions at every Gauss point."""
    element, job_results_dir = _minimal_eb_element_and_job_dir()
    try:
        obj = element.element_force_vector()
        for i, gp in enumerate(obj.gauss_data):
            assert gp.shape_functions is not None, (
                f"ForceObject element_id={obj.element_id} gauss_data[{i}]: shape_functions must be set"
            )
    finally:
        import shutil
        if os.path.isdir(job_results_dir):
            shutil.rmtree(os.path.dirname(job_results_dir), ignore_errors=True)


def test_nodal_result_projector_with_shape_functions_cache():
    """NodalResultProjector runs with formulation cache that has shape_functions and returns correct shapes."""
    import shutil
    from processing.static.results.containers import (
        GaussianResults,
        FormulationResultSet,
    )
    from processing.static.results.compute_secondary.nodal_result_projector import (
        NodalResultProjector,
    )

    element, job_results_dir = _minimal_eb_element_and_job_dir()
    try:
        elem_obj = element.element_stiffness_matrix()
        force_obj = element.element_force_vector()
        formulation_cache = FormulationResultSet(
            element_objects=[elem_obj],
            force_objects=[force_obj],
        )
        n_gauss = len(elem_obj.gauss_data)
        # Minimal Gaussian results: one element, list of per-GP values
        elem_strains = [np.ones(6) * (g + 0.1) for g in range(n_gauss)]
        elem_stresses = [np.ones(6) * (g + 0.2) for g in range(n_gauss)]
        elem_energies = [float(g + 0.3) for g in range(n_gauss)]
        gaussian_results = GaussianResults(
            strain=[elem_strains],
            stress=[elem_stresses],
            internal_energy_density=[elem_energies],
        )
        element_dictionary = {
            "ids": np.array([0]),
            "connectivity": np.array([[0, 1]]),
        }
        grid_dictionary = {"coordinates": np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])}
        elements = [element]

        projector = NodalResultProjector(
            elements=elements,
            element_dictionary=element_dictionary,
            grid_dictionary=grid_dictionary,
            gaussian_results=gaussian_results,
            formulation_cache=formulation_cache,
        )
        nodal_results = projector.project()

        assert nodal_results.strain is not None
        assert nodal_results.stress is not None
        assert nodal_results.strain_energy_density is not None
        assert nodal_results.strain.shape == (2, 6)
        assert nodal_results.stress.shape == (2, 6)
        assert nodal_results.strain_energy_density.shape == (2,)
        assert np.all(np.isfinite(nodal_results.strain))
        assert np.all(np.isfinite(nodal_results.stress))
        assert np.all(np.isfinite(nodal_results.strain_energy_density))
    finally:
        if os.path.isdir(job_results_dir):
            shutil.rmtree(os.path.dirname(job_results_dir), ignore_errors=True)

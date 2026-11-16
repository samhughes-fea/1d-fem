"""Element error estimate API from existing Gauss results. No adaptivity."""

from __future__ import annotations

import os
import shutil
import tempfile

import numpy as np
import pytest

from onedfem.processing.static.results.compute_tertiary.element_error_estimate import (
    ZZ_SOURCE,
    element_error_indicators,
)
from onedfem.processing.static.results.containers import GaussianResults, NodalResults
from onedfem.simulation_runner.static.linear_static_simulation import LinearStaticSimulationRunner


ELEMENT = "LinearEulerBernoulliBeamElement3D"


def _cantilever(n_elem: int, L: float, E: float, I: float, rho: float, tip_force):
    from onedfem.pre_processing.element_library.linear.beam.zosdt.euler_bernoulli.linear_euler_bernoulli_3D import (
        LinearEulerBernoulliBeamElement3D,
    )

    n_nodes = n_elem + 1
    dx = L / n_elem
    coordinates = np.array([[i * dx, 0.0, 0.0] for i in range(n_nodes)], dtype=np.float64)
    ids = np.arange(n_nodes, dtype=np.int32)
    section = {
        "A": np.full(n_elem, 0.01),
        "I_x": np.full(n_elem, 1.0e-9),
        "I_y": np.full(n_elem, I),
        "I_z": np.full(n_elem, I),
        "J_t": np.full(n_elem, 1.0e-9),
    }
    material = {
        "E": np.full(n_elem, E),
        "G": np.full(n_elem, E / 2.6),
        "nu": np.full(n_elem, 0.3),
        "rho": np.full(n_elem, rho),
    }
    element_dictionary = {
        "ids": np.arange(n_elem, dtype=np.int32),
        "connectivity": np.array([[i, i + 1] for i in range(n_elem)], dtype=np.int32),
        "types": np.array([ELEMENT] * n_elem),
        "integration_orders": {
            name: np.full(n_elem, 4, dtype=np.int64)
            for name in ("axial", "bending_y", "bending_z", "shear_y", "shear_z", "torsion", "load")
        },
    }
    grid = {"coordinates": coordinates, "ids": ids}
    Fx, Fy, Fz = tip_force
    point_load = np.array([[L, 0.0, 0.0, Fx, Fy, Fz, 0.0, 0.0, 0.0]], dtype=np.float64)
    distributed = np.empty((0, 9))
    temp = tempfile.mkdtemp()
    job = os.path.join(temp, "job")
    os.makedirs(os.path.join(job, "element_stiffness_matrices"), exist_ok=True)
    os.makedirs(os.path.join(job, "element_force_vectors"), exist_ok=True)
    elements, element_objects, force_objects = [], [], []
    for eid in range(n_elem):
        element = LinearEulerBernoulliBeamElement3D(
            element_id=eid,
            element_dictionary=element_dictionary,
            grid_dictionary=grid,
            material_dictionary=material,
            section_dictionary=section,
            point_load_array=point_load,
            distributed_load_array=distributed,
            job_results_dir=job,
        )
        elements.append(element)
        element_objects.append(element.element_stiffness_matrix())
        force_objects.append(element.element_force_vector())
    return {
        "elements": elements,
        "element_objects": element_objects,
        "force_objects": force_objects,
        "temp": temp,
        "job": job,
        "grid": grid,
        "element_dictionary": element_dictionary,
        "material": material,
        "section": section,
        "point_load": point_load,
        "distributed": distributed,
        "n_elem": n_elem,
        "n_nodes": n_nodes,
    }


def _static_runner(built, job_name: str) -> LinearStaticSimulationRunner:
    return LinearStaticSimulationRunner(
        elements=built["elements"],
        grid_dictionary=built["grid"],
        element_dictionary=built["element_dictionary"],
        material_dictionary=built["material"],
        section_dictionary=built["section"],
        point_load_array=built["point_load"],
        distributed_load_array=built["distributed"],
        element_objects=built["element_objects"],
        force_objects=built["force_objects"],
        job_name=job_name,
        job_results_dir=built["job"],
        simulation_settings={"solver": {"type": "direct", "tolerance": 1e-12, "max_iterations": 500}},
    )


def test_missing_gauss_stress_returns_none():
    assert element_error_indicators(None, None, None, {}) is None
    gauss = GaussianResults(stress=None)
    nodal = NodalResults(stress=np.zeros((2, 6)))
    assert element_error_indicators(gauss, object(), nodal, {"connectivity": [[0, 1]]}) is None


def test_zero_fields_give_zero_indicators():
    built = _cantilever(2, 1.0, 210e9, 1e-6, 7850.0, (0.0, 0.0, 1.0))
    try:
        result = _static_runner(built, "error_zero").run()
        gauss = result.secondary.gaussian_results
        n_elem = built["n_elem"]
        n_gp = len(gauss.stress[0])
        n_comp = len(gauss.stress[0][0])
        zero_gauss = GaussianResults(
            stress=[[np.zeros(n_comp) for _ in range(n_gp)] for _ in range(n_elem)],
            internal_energy_density=[[0.0] * n_gp for _ in range(n_elem)],
        )
        zero_nodal = NodalResults(stress=np.zeros_like(result.secondary.nodal_results.stress))
        estimate = element_error_indicators(
            zero_gauss, _formulation(built), zero_nodal, built["element_dictionary"]
        )
        assert estimate is not None
        assert np.allclose(estimate.eta_e, 0.0)
        assert estimate.eta_global == pytest.approx(0.0)
        assert "Zienkiewicz" in ZZ_SOURCE
        assert "No adaptivity" in ZZ_SOURCE
    finally:
        shutil.rmtree(built["temp"], ignore_errors=True)


def _formulation(built):
    from onedfem.processing.static.results.containers.formulation_results import FormulationResultSet

    return FormulationResultSet(element_objects=built["element_objects"], force_objects=built["force_objects"])


def test_static_gauss_results_attach_a_nonnegative_indicator():
    L, E, I, P = 2.0, 210.0e9, 1.0e-6, 1000.0
    built = _cantilever(4, L, E, I, 7850.0, (0.0, 0.0, P))
    try:
        result = _static_runner(built, "error_estimate").run()
        assert result.tertiary.element_error_eta is not None
        eta = np.asarray(result.tertiary.element_error_eta, dtype=np.float64)
        assert eta.shape == (built["n_elem"],)
        assert np.all(np.isfinite(eta))
        assert np.all(eta >= 0.0)
        assert result.tertiary.element_error_eta_global == pytest.approx(float(np.linalg.norm(eta)))
        path = os.path.join(built["job"], "tertiary_results", "elemental", "element_error_estimate.csv")
        assert os.path.isfile(path)
        loaded = np.loadtxt(path, delimiter=",", skiprows=1)
        np.testing.assert_allclose(loaded[:, 1], eta)
        cache = _formulation(built)
        estimate = element_error_indicators(
            result.secondary.gaussian_results,
            cache,
            result.secondary.nodal_results,
            built["element_dictionary"],
        )
        np.testing.assert_allclose(estimate.eta_e, eta, rtol=1e-10)
    finally:
        shutil.rmtree(built["temp"], ignore_errors=True)

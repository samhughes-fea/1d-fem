"""Gauss coordinates, static result object, model totals, and modal element energy.

The mesh is the Euler-Bernoulli cantilever already used for limit-state checks.
Nothing here invents a published number: energy is compared with ``0.5 Uᵀ F``,
Gauss locations with the element axis, and modal element energy with
``0.5 φᵀ K φ`` assembled from the same ``K_e``.
"""

from __future__ import annotations

import os
import shutil
import tempfile

import numpy as np
import pytest

from onedfem.post_processing.stations import axis_point
from onedfem.processing.static.results.containers import ModelTotals, StaticResult
from onedfem.simulation_runner.eigen.eigen_simulation import EigenSimulationRunner
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
    elements, stiffness, mass, element_objects, force_objects = [], [], [], [], []
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
        stiffness.append(np.asarray(element_objects[-1].K_e, dtype=np.float64))
        mass.append(np.asarray(element.element_mass_matrix().M_e, dtype=np.float64))
    return {
        "elements": elements,
        "coordinates": coordinates,
        "stiffness": stiffness,
        "mass": mass,
        "element_objects": element_objects,
        "force_objects": force_objects,
        "temp": temp,
        "job": job,
        "grid": grid,
        "mesh": {"node_ids": ids, "coordinates": coordinates},
        "element_dictionary": element_dictionary,
        "material": material,
        "section": section,
        "point_load": point_load,
        "distributed": distributed,
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


def _eigen_runner(built, job_name: str, n_modes: int) -> EigenSimulationRunner:
    return EigenSimulationRunner(
        settings={
            "elements": np.asarray(built["elements"], dtype=object),
            "mesh_dictionary": built["mesh"],
            "grid_dictionary": built["grid"],
            "element_dictionary": built["element_dictionary"],
            "material_dictionary": built["material"],
            "section_dictionary": built["section"],
            "point_load_array": built["point_load"],
            "distributed_load_array": built["distributed"],
            "job_results_dir": built["job"],
            "simulation_settings": {"eigen": {"num_modes": n_modes}},
            "element_stiffness_matrices": np.stack(built["stiffness"]),
            "element_mass_matrices": np.stack(built["mass"]),
            "element_objects": np.asarray(built["element_objects"], dtype=object),
            "force_objects": np.asarray(built["force_objects"], dtype=object),
            "prescribed_displacement_dict": None,
        },
        job_name=job_name,
    )


def test_gaussian_results_fill_physical_coords_and_weights():
    """Each Gauss point sits on its element axis; weights are the quadrature ``w_g``."""
    L, E, I, P = 2.0, 210.0e9, 1.0e-6, 1000.0
    built = _cantilever(4, L, E, I, 7850.0, (0.0, 0.0, P))
    try:
        result = _static_runner(built, "gauss_coords").run()
        gauss = result.secondary.gaussian_results
        assert gauss.gauss_coords is not None
        assert gauss.gauss_weights is not None
        assert len(gauss.gauss_coords) == len(built["elements"])
        for element, coords, weights, obj in zip(
            built["elements"], gauss.gauss_coords, gauss.gauss_weights, built["element_objects"]
        ):
            assert coords.shape == (len(obj.gauss_data), 3)
            assert weights.shape == (len(obj.gauss_data),)
            for k, gp in enumerate(obj.gauss_data):
                np.testing.assert_allclose(coords[k], axis_point(element, gp.xi), atol=1e-12)
                assert weights[k] == pytest.approx(float(gp.weight))
                assert coords[k, 1] == pytest.approx(0.0, abs=1e-12)
                assert coords[k, 2] == pytest.approx(0.0, abs=1e-12)
                assert 0.0 <= coords[k, 0] <= L
        path = os.path.join(built["job"], "secondary_results", "gaussian", "gauss_points.csv")
        header = open(path, encoding="utf-8").readline()
        for name in ("element_id", "gauss_point", "xi", "weight", "x", "y", "z"):
            assert name in header
    finally:
        shutil.rmtree(built["temp"], ignore_errors=True)


def test_static_run_returns_the_result_object_and_still_writes_the_tree():
    L, E, I, P = 2.0, 210.0e9, 1.0e-6, 1000.0
    built = _cantilever(4, L, E, I, 7850.0, (0.0, 0.0, P))
    try:
        runner = _static_runner(built, "static_result")
        result = runner.run()
        assert isinstance(result, StaticResult)
        assert result.displacements.shape == (built["coordinates"].shape[0] * 6,)
        np.testing.assert_allclose(result.displacements, runner.U_global.ravel())
        tip = P * L**3 / (3.0 * E * I)
        np.testing.assert_allclose(result.displacements[-4], tip, rtol=1e-8, atol=1e-12)
        assert result.primary.global_results.U_global is not None
        assert result.secondary.gaussian_results.strain
        assert result.tertiary.total_strain_energy
        assert os.path.isfile(os.path.join(built["job"], "primary_results", "totals.csv"))
        assert os.path.isdir(os.path.join(built["job"], "secondary_results", "gaussian"))
        assert os.path.isdir(os.path.join(built["job"], "tertiary_results", "elemental"))
    finally:
        shutil.rmtree(built["temp"], ignore_errors=True)


def test_model_totals_match_element_energy_external_work_and_free_residual():
    L, E, I, P = 2.0, 210.0e9, 1.0e-6, 1000.0
    built = _cantilever(4, L, E, I, 7850.0, (0.0, 0.0, P))
    try:
        result = _static_runner(built, "model_totals").run()
        totals = result.totals
        assert isinstance(totals, ModelTotals)
        element_sum = float(np.sum(result.tertiary.total_strain_energy))
        U = result.displacements
        F = np.asarray(result.primary.global_results.F_global, dtype=np.float64).ravel()
        assert totals.strain_energy == pytest.approx(element_sum)
        assert totals.external_work == pytest.approx(0.5 * float(U @ F))
        assert totals.strain_energy == pytest.approx(totals.external_work, rel=1e-8)
        analytic = P**2 * L**3 / (6.0 * E * I)
        assert totals.strain_energy == pytest.approx(analytic, rel=1e-8)
        assert totals.equilibrium_residual < 1e-6
        rows = {}
        with open(os.path.join(built["job"], "primary_results", "totals.csv"), encoding="utf-8") as handle:
            next(handle)
            for line in handle:
                name, value = line.strip().split(",")
                rows[name] = float(value)
        assert rows["strain_energy"] == pytest.approx(totals.strain_energy)
        assert rows["external_work"] == pytest.approx(totals.external_work)
        assert rows["equilibrium_residual"] == pytest.approx(totals.equilibrium_residual)
    finally:
        shutil.rmtree(built["temp"], ignore_errors=True)


def test_per_element_modal_strain_energy_sums_to_the_mode_energy():
    L, E, I, rho = 2.0, 70.0e9, 1.0e-5, 2700.0
    built = _cantilever(8, L, E, I, rho, (0.0, 0.0, 0.0))
    try:
        runner = _eigen_runner(built, "modal_energy", n_modes=2)
        runner.run()
        energy = np.asarray(runner.modal_element_strain_energy, dtype=np.float64)
        phi = np.asarray(runner.vibration_mode_shapes, dtype=np.float64)
        assert energy.shape == (len(built["elements"]), 2)
        ndof = built["coordinates"].shape[0] * 6
        K = np.zeros((ndof, ndof))
        for element, Ke in zip(built["elements"], built["stiffness"]):
            dof = np.asarray(element.assemble_global_dof_indices(), dtype=int)
            K[np.ix_(dof, dof)] += Ke
        for j in range(2):
            assembled = 0.5 * float(phi[:, j] @ K @ phi[:, j])
            assert float(energy[:, j].sum()) == pytest.approx(assembled, rel=1e-10, abs=0.0)
        path = os.path.join(
            built["job"],
            "primary_results",
            "modal_results",
            "modal_energy_modal_element_strain_energy.csv",
        )
        assert os.path.isfile(path)
        loaded = np.loadtxt(path, delimiter=",", skiprows=1)
        np.testing.assert_allclose(loaded[:, 1:], energy)
        gm_path = os.path.join(
            built["job"], "primary_results", "modal_results", "modal_energy_modal_generalized_mass.txt"
        )
        assert os.path.isfile(gm_path)
    finally:
        shutil.rmtree(built["temp"], ignore_errors=True)

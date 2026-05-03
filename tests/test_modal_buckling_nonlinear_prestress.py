"""Linear buckling with nonlinear_static prestress (nonlinear equilibrium path then K_sigma)."""

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


def _build_nl_ts_cantilever_buckling(n_elem: int, L: float, E: float, P_ref: float, prestress: str):
    from pre_processing.element_library.nonlinear.timoshenko.nonlinear_timoshenko_3D import (
        NonlinearTimoshenkoBeamElement3D,
    )

    num_nodes = n_elem + 1
    dx = L / n_elem
    coordinates = np.array([[i * dx, 0.0, 0.0] for i in range(num_nodes)], dtype=np.float64)
    ids = np.arange(num_nodes, dtype=np.int32)
    n_el = n_elem
    A = 0.02
    I_val = 2.0e-6
    section_dictionary = {
        "A": np.full(n_el, A),
        "I_x": np.full(n_el, 1.0e-9),
        "I_y": np.full(n_el, I_val),
        "I_z": np.full(n_el, I_val),
        "J_t": np.full(n_el, 2.0e-9),
        "kappa": np.full(n_el, 5.0 / 6.0),
        "alpha": np.full(n_el, 0.0),
    }
    material_dictionary = {
        "E": np.full(n_el, E),
        "G": np.full(n_el, E / 2.6),
        "nu": np.full(n_el, 0.3),
        "rho": np.full(n_el, 7850.0),
    }
    element_dictionary = {
        "ids": np.arange(n_el, dtype=np.int32),
        "connectivity": np.array([[i, i + 1] for i in range(n_el)], dtype=np.int32),
        "types": np.array(["NonlinearTimoshenkoBeamElement3D"] * n_el),
        "integration_orders": {
            "axial": np.full(n_el, 3, dtype=np.int64),
            "bending_y": np.full(n_el, 3, dtype=np.int64),
            "bending_z": np.full(n_el, 3, dtype=np.int64),
            "shear_y": np.full(n_el, 2, dtype=np.int64),
            "shear_z": np.full(n_el, 2, dtype=np.int64),
            "torsion": np.full(n_el, 3, dtype=np.int64),
            "load": np.full(n_el, 2, dtype=np.int64),
        },
    }
    grid_dictionary = {"coordinates": coordinates, "ids": ids}
    point_load_array = np.array(
        [[L, 0.0, 0.0, -float(P_ref), 0.0, 0.0, 0.0, 0.0, 0.0]],
        dtype=np.float64,
    )
    distributed_load_array = np.empty((0, 9))
    temp_dir = tempfile.mkdtemp()
    job_results_dir = os.path.join(temp_dir, "buckling_nl_pre")
    os.makedirs(job_results_dir, exist_ok=True)
    os.makedirs(os.path.join(job_results_dir, "element_stiffness_matrices"), exist_ok=True)
    os.makedirs(os.path.join(job_results_dir, "element_force_vectors"), exist_ok=True)
    elements = []
    stiffness_objs = []
    mass_objs = []
    force_objs = []
    for eid in range(n_el):
        elem = NonlinearTimoshenkoBeamElement3D(
            element_id=eid,
            element_dictionary=element_dictionary,
            grid_dictionary=grid_dictionary,
            material_dictionary=material_dictionary,
            section_dictionary=section_dictionary,
            point_load_array=point_load_array,
            distributed_load_array=distributed_load_array,
            job_results_dir=job_results_dir,
        )
        elements.append(elem)
        stiffness_objs.append(elem.element_stiffness_matrix())
        mass_objs.append(elem.element_mass_matrix())
        force_objs.append(elem.element_force_vector())
    mesh_dictionary = {"node_ids": ids, "coordinates": coordinates}
    simulation_settings = {
        "modal": {
            "num_modes": 3,
            "analysis": "buckling",
            "buckling_prestress": prestress,
            "buckling_load_factor": 1.0,
        },
        "newton": {
            "tolerance": 1e-6,
            "max_iterations": 80,
            "tolerance_delta_u": 1e-9,
            "relative_tolerance": None,
            "relative_reference": "first_residual",
        },
        "nonlinear": {"num_increments": 1},
    }
    settings = {
        "elements": np.asarray(elements, dtype=object),
        "mesh_dictionary": mesh_dictionary,
        "grid_dictionary": grid_dictionary,
        "element_dictionary": element_dictionary,
        "material_dictionary": material_dictionary,
        "section_dictionary": section_dictionary,
        "point_load_array": point_load_array,
        "distributed_load_array": distributed_load_array,
        "job_results_dir": job_results_dir,
        "simulation_settings": simulation_settings,
        "element_stiffness_matrices": np.stack(
            [np.asarray(o.K_e).astype(np.float64) for o in stiffness_objs]
        ),
        "element_mass_matrices": np.stack([np.asarray(o.M_e).astype(np.float64) for o in mass_objs]),
        "element_objects": np.asarray(stiffness_objs, dtype=object),
        "force_objects": np.asarray(force_objs, dtype=object),
        "prescribed_displacement_dict": None,
    }
    return settings, temp_dir


def test_modal_buckling_nonlinear_static_prestress_smoke():
    from simulation_runner.buckling.buckling_simulation import LinearBucklingSimulationRunner

    L = 2.0
    E = 200.0e9
    P_ref = 1.0
    settings, td = _build_nl_ts_cantilever_buckling(10, L, E, P_ref, "nonlinear_static")
    try:
        runner = LinearBucklingSimulationRunner(settings=settings, job_name="nl_prestress_smoke")
        runner.run()
        lam = np.asarray(runner.secondary_results["global"]["buckling_load_factors"], dtype=np.float64).ravel()
        assert np.all(np.isfinite(lam))
        assert np.any(lam > 0)
    finally:
        shutil.rmtree(td, ignore_errors=True)


def test_modal_buckling_linear_element_rejects_nonlinear_prestress():
    from pre_processing.element_library.linear.beam.first_order_shear_deformation_theory.timoshenko.linear_timoshenko_3D import (
        LinearTimoshenkoBeamElement3D,
    )
    from simulation_runner.buckling.buckling_simulation import LinearBucklingSimulationRunner

    L = 1.0
    E = 200.0e9
    P_ref = 1.0
    num_nodes = 3
    n_el = 2
    coordinates = np.array([[0.0, 0.0, 0.0], [L / 2, 0.0, 0.0], [L, 0.0, 0.0]], dtype=np.float64)
    ids = np.arange(num_nodes, dtype=np.int32)
    section_dictionary = {
        "A": np.full(n_el, 0.02),
        "I_x": np.full(n_el, 1e-9),
        "I_y": np.full(n_el, 1e-6),
        "I_z": np.full(n_el, 1e-6),
        "J_t": np.full(n_el, 2e-9),
        "kappa": np.full(n_el, 5.0 / 6.0),
        "alpha": np.full(n_el, 0.0),
    }
    material_dictionary = {
        "E": np.full(n_el, E),
        "G": np.full(n_el, E / 2.6),
        "nu": np.full(n_el, 0.3),
        "rho": np.full(n_el, 7850.0),
    }
    element_dictionary = {
        "ids": np.arange(n_el, dtype=np.int32),
        "connectivity": np.array([[0, 1], [1, 2]], dtype=np.int32),
        "types": np.array(["LinearTimoshenkoBeamElement3D"] * n_el),
        "integration_orders": {
            "axial": np.full(n_el, 3, dtype=np.int64),
            "bending_y": np.full(n_el, 3, dtype=np.int64),
            "bending_z": np.full(n_el, 3, dtype=np.int64),
            "shear_y": np.full(n_el, 2, dtype=np.int64),
            "shear_z": np.full(n_el, 2, dtype=np.int64),
            "torsion": np.full(n_el, 3, dtype=np.int64),
            "load": np.full(n_el, 2, dtype=np.int64),
        },
    }
    grid_dictionary = {"coordinates": coordinates, "ids": ids}
    point_load_array = np.array([[L, 0.0, 0.0, -float(P_ref), 0.0, 0.0, 0.0, 0.0, 0.0]], dtype=np.float64)
    distributed_load_array = np.empty((0, 9))
    temp_dir = tempfile.mkdtemp()
    job_results_dir = os.path.join(temp_dir, "buckling_lin")
    os.makedirs(job_results_dir, exist_ok=True)
    os.makedirs(os.path.join(job_results_dir, "element_stiffness_matrices"), exist_ok=True)
    os.makedirs(os.path.join(job_results_dir, "element_force_vectors"), exist_ok=True)
    stiffness_objs = []
    mass_objs = []
    force_objs = []
    elements = []
    for eid in range(n_el):
        elem = LinearTimoshenkoBeamElement3D(
            element_id=eid,
            element_dictionary=element_dictionary,
            grid_dictionary=grid_dictionary,
            material_dictionary=material_dictionary,
            section_dictionary=section_dictionary,
            point_load_array=point_load_array,
            distributed_load_array=distributed_load_array,
            job_results_dir=job_results_dir,
        )
        elements.append(elem)
        stiffness_objs.append(elem.element_stiffness_matrix())
        mass_objs.append(elem.element_mass_matrix())
        force_objs.append(elem.element_force_vector())
    mesh_dictionary = {"node_ids": ids, "coordinates": coordinates}
    simulation_settings = {
        "modal": {
            "num_modes": 2,
            "analysis": "buckling",
            "buckling_prestress": "nonlinear_static",
            "buckling_load_factor": 1.0,
        },
    }
    settings = {
        "elements": np.asarray(elements, dtype=object),
        "mesh_dictionary": mesh_dictionary,
        "grid_dictionary": grid_dictionary,
        "element_dictionary": element_dictionary,
        "material_dictionary": material_dictionary,
        "section_dictionary": section_dictionary,
        "point_load_array": point_load_array,
        "distributed_load_array": distributed_load_array,
        "job_results_dir": job_results_dir,
        "simulation_settings": simulation_settings,
        "element_stiffness_matrices": np.stack(
            [np.asarray(o.K_e).astype(np.float64) for o in stiffness_objs]
        ),
        "element_mass_matrices": np.stack([np.asarray(o.M_e).astype(np.float64) for o in mass_objs]),
        "element_objects": np.asarray(stiffness_objs, dtype=object),
        "force_objects": np.asarray(force_objs, dtype=object),
        "prescribed_displacement_dict": None,
    }
    try:
        runner = LinearBucklingSimulationRunner(settings=settings, job_name="reject_nl")
        with pytest.raises(RuntimeError, match="Buckling simulation aborted"):
            runner.run()
    finally:
        shutil.rmtree(temp_dir, ignore_errors=True)


def test_modal_buckling_nonlinear_static_with_linear_twins_smoke():
    """nonlinear_static prestress using modal.buckling_nonlinear_prestress_twins on a linear Timoshenko mesh."""
    from pre_processing.element_library.linear.beam.first_order_shear_deformation_theory.timoshenko.linear_timoshenko_3D import (
        LinearTimoshenkoBeamElement3D,
    )
    from simulation_runner.buckling.buckling_simulation import LinearBucklingSimulationRunner

    L = 2.0
    E = 200.0e9
    P_ref = 1.0
    n_elem = 8
    num_nodes = n_elem + 1
    dx = L / n_elem
    coordinates = np.array([[i * dx, 0.0, 0.0] for i in range(num_nodes)], dtype=np.float64)
    ids = np.arange(num_nodes, dtype=np.int32)
    n_el = n_elem
    A = 0.02
    I_val = 2.0e-6
    section_dictionary = {
        "A": np.full(n_el, A),
        "I_x": np.full(n_el, 1.0e-9),
        "I_y": np.full(n_el, I_val),
        "I_z": np.full(n_el, I_val),
        "J_t": np.full(n_el, 2.0e-9),
        "kappa": np.full(n_el, 5.0 / 6.0),
        "alpha": np.full(n_el, 0.0),
    }
    material_dictionary = {
        "E": np.full(n_el, E),
        "G": np.full(n_el, E / 2.6),
        "nu": np.full(n_el, 0.3),
        "rho": np.full(n_el, 7850.0),
    }
    element_dictionary = {
        "ids": np.arange(n_el, dtype=np.int32),
        "connectivity": np.array([[i, i + 1] for i in range(n_el)], dtype=np.int32),
        "types": np.array(["LinearTimoshenkoBeamElement3D"] * n_el),
        "integration_orders": {
            "axial": np.full(n_el, 3, dtype=np.int64),
            "bending_y": np.full(n_el, 3, dtype=np.int64),
            "bending_z": np.full(n_el, 3, dtype=np.int64),
            "shear_y": np.full(n_el, 2, dtype=np.int64),
            "shear_z": np.full(n_el, 2, dtype=np.int64),
            "torsion": np.full(n_el, 3, dtype=np.int64),
            "load": np.full(n_el, 2, dtype=np.int64),
        },
    }
    grid_dictionary = {"coordinates": coordinates, "ids": ids}
    point_load_array = np.array(
        [[L, 0.0, 0.0, -float(P_ref), 0.0, 0.0, 0.0, 0.0, 0.0]],
        dtype=np.float64,
    )
    distributed_load_array = np.empty((0, 9))
    temp_dir = tempfile.mkdtemp()
    job_results_dir = os.path.join(temp_dir, "buckling_nl_twins")
    os.makedirs(job_results_dir, exist_ok=True)
    os.makedirs(os.path.join(job_results_dir, "element_stiffness_matrices"), exist_ok=True)
    os.makedirs(os.path.join(job_results_dir, "element_force_vectors"), exist_ok=True)
    elements = []
    stiffness_objs = []
    mass_objs = []
    force_objs = []
    for eid in range(n_el):
        elem = LinearTimoshenkoBeamElement3D(
            element_id=eid,
            element_dictionary=element_dictionary,
            grid_dictionary=grid_dictionary,
            material_dictionary=material_dictionary,
            section_dictionary=section_dictionary,
            point_load_array=point_load_array,
            distributed_load_array=distributed_load_array,
            job_results_dir=job_results_dir,
        )
        elements.append(elem)
        stiffness_objs.append(elem.element_stiffness_matrix())
        mass_objs.append(elem.element_mass_matrix())
        force_objs.append(elem.element_force_vector())
    mesh_dictionary = {"node_ids": ids, "coordinates": coordinates}
    simulation_settings = {
        "modal": {
            "num_modes": 2,
            "analysis": "buckling",
            "buckling_prestress": "nonlinear_static",
            "buckling_load_factor": 1.0,
            "buckling_nonlinear_prestress_twins": True,
        },
        "newton": {
            "tolerance": 1e-6,
            "max_iterations": 80,
            "tolerance_delta_u": 1e-9,
            "relative_tolerance": None,
            "relative_reference": "first_residual",
        },
        "nonlinear": {"num_increments": 1},
    }
    settings = {
        "elements": np.asarray(elements, dtype=object),
        "mesh_dictionary": mesh_dictionary,
        "grid_dictionary": grid_dictionary,
        "element_dictionary": element_dictionary,
        "material_dictionary": material_dictionary,
        "section_dictionary": section_dictionary,
        "point_load_array": point_load_array,
        "distributed_load_array": distributed_load_array,
        "job_results_dir": job_results_dir,
        "simulation_settings": simulation_settings,
        "element_stiffness_matrices": np.stack(
            [np.asarray(o.K_e).astype(np.float64) for o in stiffness_objs]
        ),
        "element_mass_matrices": np.stack([np.asarray(o.M_e).astype(np.float64) for o in mass_objs]),
        "element_objects": np.asarray(stiffness_objs, dtype=object),
        "force_objects": np.asarray(force_objs, dtype=object),
        "prescribed_displacement_dict": None,
    }
    try:
        runner = LinearBucklingSimulationRunner(settings=settings, job_name="nl_twins_smoke")
        runner.run()
        lam = np.asarray(runner.secondary_results["global"]["buckling_load_factors"], dtype=np.float64).ravel()
        assert np.all(np.isfinite(lam))
        assert np.any(lam > 0)
    finally:
        shutil.rmtree(temp_dir, ignore_errors=True)

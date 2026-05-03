"""
Regression: linear buckling critical load vs Euler cantilever formula.

Uses ``LinearBucklingSimulationRunner`` with canonical or legacy buckling settings, EB beams without warping,
cantilever base fixed by default modal BCs (first six global DOFs), compressive tip load.
"""

import os
import shutil
import sys
import tempfile
from pathlib import Path

import numpy as np
import pytest

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))


def _euler_cantilever_P_cr(E: float, I: float, L: float) -> float:
    """Critical axial load (Euler), fixed-free, same bending stiffness in both lateral directions."""
    return (np.pi ** 2) * E * I / (4.0 * L ** 2)


def _build_cantilever_modal_case(n_elem: int, L: float, E: float, P_ref: float):
    """Return (settings dict for LinearBucklingSimulationRunner, temp_dir to rmtree)."""
    from pre_processing.element_library.linear.beam.zero_order_shear_deformation_theory.euler_bernoulli.linear_euler_bernoulli_3D import (
        LinearEulerBernoulliBeamElement3D,
    )

    num_nodes = n_elem + 1
    dx = L / n_elem
    coordinates = np.array([[i * dx, 0.0, 0.0] for i in range(num_nodes)], dtype=np.float64)
    ids = np.arange(num_nodes, dtype=np.int32)

    A = 0.01
    I_val = 1.0e-6
    n_el = n_elem
    # Section/material arrays are indexed by element row index (same length as element_dictionary ids).
    section_dictionary = {
        "A": np.full(n_el, A),
        "I_x": np.full(n_el, 1.0e-9),
        "I_y": np.full(n_el, I_val),
        "I_z": np.full(n_el, I_val),
        "J_t": np.full(n_el, 1.0e-9),
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
        "types": np.array(["LinearEulerBernoulliBeamElement3D"] * n_el),
        "integration_orders": {
            "axial": np.full(n_el, 4, dtype=np.int64),
            "bending_y": np.full(n_el, 4, dtype=np.int64),
            "bending_z": np.full(n_el, 4, dtype=np.int64),
            "shear_y": np.full(n_el, 3, dtype=np.int64),
            "shear_z": np.full(n_el, 3, dtype=np.int64),
            "torsion": np.full(n_el, 4, dtype=np.int64),
            "load": np.full(n_el, 3, dtype=np.int64),
        },
    }
    grid_dictionary = {"coordinates": coordinates, "ids": ids}

    # Tip compression at x = L (global); columns [x, y, z, Fx, Fy, Fz, Mx, My, Mz]
    point_load_array = np.array(
        [[L, 0.0, 0.0, -float(P_ref), 0.0, 0.0, 0.0, 0.0, 0.0]],
        dtype=np.float64,
    )
    distributed_load_array = np.empty((0, 9))

    temp_dir = tempfile.mkdtemp()
    job_results_dir = os.path.join(temp_dir, "buckling_job")
    os.makedirs(job_results_dir, exist_ok=True)
    os.makedirs(os.path.join(job_results_dir, "element_stiffness_matrices"), exist_ok=True)
    os.makedirs(os.path.join(job_results_dir, "element_force_vectors"), exist_ok=True)

    elements = []
    stiffness_objs = []
    mass_objs = []
    force_objs = []
    for eid in range(n_el):
        elem = LinearEulerBernoulliBeamElement3D(
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

    mesh_dictionary = {
        "node_ids": ids,
        "coordinates": coordinates,
    }

    simulation_settings = {
        "modal": {
            "num_modes": 4,
            "analysis": "buckling",
            "buckling_prestress": "linear_static",
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
        # np.asarray(..., dtype=object) stacks equal-shaped matrices into (n,12,12) with dtype object;
        # modal assembly iterates the first axis — use float np.stack instead.
        "element_stiffness_matrices": np.stack(
            [np.asarray(o.K_e).astype(np.float64) for o in stiffness_objs]
        ),
        "element_mass_matrices": np.stack(
            [np.asarray(o.M_e).astype(np.float64) for o in mass_objs]
        ),
        "element_objects": np.asarray(stiffness_objs, dtype=object),
        "force_objects": np.asarray(force_objs, dtype=object),
        "prescribed_displacement_dict": None,
    }
    return settings, temp_dir, I_val


@pytest.mark.parametrize("n_elem", [16, 32])
def test_modal_buckling_cantilever_matches_euler(n_elem: int):
    """Lowest buckling load factor × reference tip compressive load vs Euler cantilever P_cr."""
    from simulation_runner.buckling.buckling_simulation import LinearBucklingSimulationRunner

    L = 2.5
    E = 210.0e9
    P_ref = 1.0
    settings, temp_dir, I_val = _build_cantilever_modal_case(n_elem, L, E, P_ref)
    try:
        runner = LinearBucklingSimulationRunner(settings=settings, job_name="buckling_euler_test")
        runner.run()
        lambdas = np.asarray(
            runner.secondary_results["global"]["buckling_load_factors"], dtype=np.float64
        ).ravel()
        assert lambdas.size >= 2
        P_theory = _euler_cantilever_P_cr(E, I_val, L)
        P_modes = lambdas * P_ref
        # Isotropic I_y = I_z: two orthogonal lateral buckling planes share the same critical load;
        # the generalized eigen solver may split degenerate pairs — match theory on any of the first modes.
        rtol = 0.06 if n_elem == 16 else 0.03
        err = np.min(np.abs(P_modes[:8] - P_theory))
        np.testing.assert_allclose(err, 0.0, atol=rtol * P_theory)
    finally:
        shutil.rmtree(temp_dir, ignore_errors=True)


def test_buckling_load_factor_scales_eigenvalues():
    """Doubling ``buckling_load_factor`` halves λ: same physical critical load λ·(factor·P_ref)."""
    from simulation_runner.buckling.buckling_simulation import LinearBucklingSimulationRunner

    L = 2.5
    E = 210.0e9
    P_ref = 1.0
    n_elem = 24
    settings1, td1, _ = _build_cantilever_modal_case(n_elem, L, E, P_ref)
    settings2, td2, _ = _build_cantilever_modal_case(n_elem, L, E, P_ref)
    settings2["simulation_settings"]["modal"]["buckling_load_factor"] = 2.0
    try:
        r1 = LinearBucklingSimulationRunner(settings=settings1, job_name="bf1")
        r1.run()
        r2 = LinearBucklingSimulationRunner(settings=settings2, job_name="bf2")
        r2.run()
        lam1 = np.asarray(r1.secondary_results["global"]["buckling_load_factors"], dtype=np.float64).ravel()
        lam2 = np.asarray(r2.secondary_results["global"]["buckling_load_factors"], dtype=np.float64).ravel()
        assert lam1.size >= 4 and lam2.size >= 4
        s1 = np.sort(lam1[:6])
        s2 = np.sort(lam2[:6])
        np.testing.assert_allclose(s2, 0.5 * s1, rtol=0.05)
    finally:
        shutil.rmtree(td1, ignore_errors=True)
        shutil.rmtree(td2, ignore_errors=True)


def test_buckling_prestress_none_raises():
    from simulation_runner.buckling.buckling_simulation import LinearBucklingSimulationRunner

    L = 2.5
    E = 210.0e9
    P_ref = 1.0
    settings, td, _ = _build_cantilever_modal_case(8, L, E, P_ref)
    settings["simulation_settings"]["modal"]["buckling_prestress"] = "none"
    try:
        runner = LinearBucklingSimulationRunner(settings=settings, job_name="none_prestress")
        with pytest.raises(RuntimeError, match="Buckling simulation aborted") as exc_info:
            runner.run()
        assert isinstance(exc_info.value.__cause__, ValueError)
        assert "buckling_prestress" in str(exc_info.value.__cause__)
    finally:
        shutil.rmtree(td, ignore_errors=True)


def test_modal_buckling_timoshenko_smoke():
    """Linear Timoshenko mesh supports ``linear_geometric_stiffness_matrix`` — buckling run completes."""
    from pre_processing.element_library.linear.beam.first_order_shear_deformation_theory.timoshenko.linear_timoshenko_3D import (
        LinearTimoshenkoBeamElement3D,
    )
    from simulation_runner.buckling.buckling_simulation import LinearBucklingSimulationRunner

    L = 2.0
    E = 200.0e9
    P_ref = 1.0
    n_elem = 12
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
    job_results_dir = os.path.join(temp_dir, "buckling_ts")
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
            "num_modes": 3,
            "analysis": "buckling",
            "buckling_prestress": "linear_static",
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
        "element_mass_matrices": np.stack(
            [np.asarray(o.M_e).astype(np.float64) for o in mass_objs]
        ),
        "element_objects": np.asarray(stiffness_objs, dtype=object),
        "force_objects": np.asarray(force_objs, dtype=object),
        "prescribed_displacement_dict": None,
    }
    try:
        runner = LinearBucklingSimulationRunner(settings=settings, job_name="ts_buckling")
        runner.run()
        lam = np.asarray(runner.secondary_results["global"]["buckling_load_factors"], dtype=np.float64).ravel()
        assert np.all(np.isfinite(lam))
        assert np.any(lam > 0)
    finally:
        shutil.rmtree(temp_dir, ignore_errors=True)


def test_modal_buckling_euler_warping_mesh_smoke():
    """7 DOF/node EB + Γ: linear buckling completes (embedded 12×12 K_σ)."""
    from pre_processing.element_library.linear.beam.zero_order_shear_deformation_theory.euler_bernoulli.linear_euler_bernoulli_3D import (
        LinearEulerBernoulliBeamElement3D,
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
    A = 0.015
    I_val = 2.0e-6
    gamma = 5.0e-8
    section_dictionary = {
        "A": np.full(n_el, A),
        "I_x": np.full(n_el, 1.0e-9),
        "I_y": np.full(n_el, I_val),
        "I_z": np.full(n_el, I_val),
        "J_t": np.full(n_el, 2.0e-9),
        "kappa": np.full(n_el, 5.0 / 6.0),
        "alpha": np.full(n_el, 0.0),
        "y_sc": np.full(n_el, 0.0),
        "z_sc": np.full(n_el, 0.0),
        "Gamma": np.full(n_el, gamma),
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
        "types": np.array(["LinearEulerBernoulliBeamElement3D"] * n_el),
        "warping": np.ones(n_el, dtype=np.int8),
        "integration_orders": {
            "axial": np.full(n_el, 4, dtype=np.int64),
            "bending_y": np.full(n_el, 4, dtype=np.int64),
            "bending_z": np.full(n_el, 4, dtype=np.int64),
            "shear_y": np.full(n_el, 3, dtype=np.int64),
            "shear_z": np.full(n_el, 3, dtype=np.int64),
            "torsion": np.full(n_el, 4, dtype=np.int64),
            "load": np.full(n_el, 3, dtype=np.int64),
        },
    }
    grid_dictionary = {"coordinates": coordinates, "ids": ids}
    point_load_array = np.array(
        [[L, 0.0, 0.0, -float(P_ref), 0.0, 0.0, 0.0, 0.0, 0.0]],
        dtype=np.float64,
    )
    distributed_load_array = np.empty((0, 9))
    temp_dir = tempfile.mkdtemp()
    job_results_dir = os.path.join(temp_dir, "buckling_warp")
    os.makedirs(job_results_dir, exist_ok=True)
    os.makedirs(os.path.join(job_results_dir, "element_stiffness_matrices"), exist_ok=True)
    os.makedirs(os.path.join(job_results_dir, "element_force_vectors"), exist_ok=True)
    elements = []
    stiffness_objs = []
    mass_objs = []
    force_objs = []
    for eid in range(n_el):
        elem = LinearEulerBernoulliBeamElement3D(
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
            "buckling_prestress": "linear_static",
            "buckling_load_factor": 1.0,
        },
    }
    # Clamp root (all 7 DOFs) and tip warping χ — avoids singular warping mechanisms (see JOB_INPUT_BEAM_WARPING).
    tip_chi = int((num_nodes - 1) * 7 + 6)
    prescribed_displacement_dict = {
        "global_dof": np.concatenate([np.arange(7, dtype=np.int32), np.array([tip_chi], dtype=np.int32)]),
        "value": np.zeros(8, dtype=np.float64),
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
        "element_mass_matrices": np.stack(
            [np.asarray(o.M_e).astype(np.float64) for o in mass_objs]
        ),
        "element_objects": np.asarray(stiffness_objs, dtype=object),
        "force_objects": np.asarray(force_objs, dtype=object),
        "prescribed_displacement_dict": prescribed_displacement_dict,
    }
    try:
        runner = LinearBucklingSimulationRunner(settings=settings, job_name="eb_warp_buckling")
        runner.run()
        lam = np.asarray(runner.secondary_results["global"]["buckling_load_factors"], dtype=np.float64).ravel()
        assert np.all(np.isfinite(lam))
        assert np.any(lam > 0)
    finally:
        shutil.rmtree(temp_dir, ignore_errors=True)


def test_processing_modal_buckling_unit_eigenpair():
    """Smoke: K + λ K_g singular pair yields positive λ from generalized eigen."""
    from scipy.sparse import csr_matrix

    from processing.buckling import solve_linear_buckling_eigenpairs

    K = csr_matrix(np.diag([1.0, 2.0, 3.0]))
    Kg = csr_matrix(np.diag([-0.5, -1.0, -1.5]))
    lam, vecs = solve_linear_buckling_eigenpairs(K, Kg, num_modes=2)
    assert lam[0] > 0
    np.testing.assert_allclose(lam[0], 2.0, rtol=1e-10)

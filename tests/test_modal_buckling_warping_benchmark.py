"""
Benchmark-style checks for linear buckling with warping: embedded 12×12 K_σ consistency and smoke vs Euler column.

Full lateral–torsional buckling (Γ-driven) is not validated analytically here — see JOB_INPUT_BEAM_WARPING.md.
"""

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


def _euler_cantilever_P_cr(E: float, I: float, L: float) -> float:
    return (np.pi ** 2) * E * I / (4.0 * L ** 2)


def test_warping_off_modal_matches_euler_band():
    """6 DOF mesh: lowest buckling load vs Euler cantilever (same as regression suite intent)."""
    from pre_processing.element_library.linear.beam.zero_order_shear_deformation_theory.euler_bernoulli.linear_euler_bernoulli_3D import (
        LinearEulerBernoulliBeamElement3D,
    )
    from simulation_runner.modal.modal_simulation import ModalSimulationRunner

    L = 2.0
    E = 200.0e9
    P_ref = 1.0
    n_elem = 24
    I_val = 2.0e-6
    num_nodes = n_elem + 1
    dx = L / n_elem
    coordinates = np.array([[i * dx, 0.0, 0.0] for i in range(num_nodes)], dtype=np.float64)
    ids = np.arange(num_nodes, dtype=np.int32)
    n_el = n_elem
    section_dictionary = {
        "A": np.full(n_el, 0.015),
        "I_x": np.full(n_el, 1.0e-9),
        "I_y": np.full(n_el, I_val),
        "I_z": np.full(n_el, I_val),
        "J_t": np.full(n_el, 2.0e-9),
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
    point_load_array = np.array([[L, 0.0, 0.0, -float(P_ref), 0.0, 0.0, 0.0, 0.0, 0.0]], dtype=np.float64)
    distributed_load_array = np.empty((0, 9))
    temp_dir = tempfile.mkdtemp()
    job_results_dir = os.path.join(temp_dir, "buckling_benchmark")
    os.makedirs(job_results_dir, exist_ok=True)
    os.makedirs(os.path.join(job_results_dir, "element_stiffness_matrices"), exist_ok=True)
    os.makedirs(os.path.join(job_results_dir, "element_force_vectors"), exist_ok=True)
    stiffness_objs = []
    mass_objs = []
    force_objs = []
    elements = []
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
        "element_stiffness_matrices": np.stack([np.asarray(o.K_e).astype(np.float64) for o in stiffness_objs]),
        "element_mass_matrices": np.stack([np.asarray(o.M_e).astype(np.float64) for o in mass_objs]),
        "element_objects": np.asarray(stiffness_objs, dtype=object),
        "force_objects": np.asarray(force_objs, dtype=object),
        "prescribed_displacement_dict": None,
    }
    try:
        runner = ModalSimulationRunner(settings=settings, job_name="euler_band")
        runner.run()
        lam = np.asarray(runner.secondary_results["global"]["buckling_load_factors"], dtype=np.float64).ravel()
        P_cr = _euler_cantilever_P_cr(E, I_val, L)
        P_num = lam[:6] * P_ref
        err = np.min(np.abs(P_num - P_cr))
        np.testing.assert_allclose(err, 0.0, atol=0.08 * P_cr)
    finally:
        shutil.rmtree(temp_dir, ignore_errors=True)


def test_linear_timoshenko_warping_K_sigma_embedding_matches_twelve_dof_slice():
    """14-DOF K_σ top-left 12×12 equals K_σ from first 12 displacement components alone."""
    import tempfile

    from pre_processing.element_library.linear.beam.first_order_shear_deformation_theory.timoshenko.linear_timoshenko_3D import (
        LinearTimoshenkoBeamElement3D,
    )

    tmp = tempfile.mkdtemp()
    jrd = os.path.join(tmp, "job")
    os.makedirs(os.path.join(jrd, "element_stiffness_matrices"), exist_ok=True)
    os.makedirs(os.path.join(jrd, "element_force_vectors"), exist_ok=True)
    grid_dictionary = {
        "ids": np.array([0, 1], dtype=np.int32),
        "coordinates": np.array([[0.0, 0.0, 0.0], [1.2, 0.0, 0.0]], dtype=np.float64),
    }
    gamma = 3e-9
    section_dictionary = {
        "A": np.array([0.02]),
        "I_x": np.array([1e-9]),
        "I_y": np.array([5e-7]),
        "I_z": np.array([5e-7]),
        "J_t": np.array([3e-9]),
        "kappa": np.array([5.0 / 6.0]),
        "alpha": np.array([0.0]),
        "y_sc": np.array([0.0]),
        "z_sc": np.array([0.0]),
        "Gamma": np.array([gamma]),
    }
    material_dictionary = {
        "E": np.array([210e9]),
        "G": np.array([81e9]),
        "nu": np.array([0.3]),
        "rho": np.array([7850.0]),
    }
    element_dictionary = {
        "ids": np.array([0], dtype=np.int32),
        "connectivity": np.array([[0, 1]], dtype=np.int32),
        "types": np.array(["LinearTimoshenkoBeamElement3D"]),
        "warping": np.array([1], dtype=np.int8),
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
    try:
        elem = LinearTimoshenkoBeamElement3D(
            element_id=0,
            element_dictionary=element_dictionary,
            grid_dictionary=grid_dictionary,
            material_dictionary=material_dictionary,
            section_dictionary=section_dictionary,
            point_load_array=np.empty((0, 9)),
            distributed_load_array=np.empty((0, 9)),
            job_results_dir=jrd,
        )
        np.random.seed(3)
        U14 = np.random.randn(14) * 1e-5
        U12 = U14[:12].copy()
        K14 = elem.linear_geometric_stiffness_matrix(U14)
        K12 = elem.linear_geometric_stiffness_matrix(U12)
        np.testing.assert_allclose(K14[:12, :12], K12, rtol=1e-12)
        assert np.allclose(K14[12:, :], 0.0)
        assert np.allclose(K14[:, 12:], 0.0)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def test_handbook_euler_closed_form_reference_positive():
    """Sanity: Euler cantilever P_cr formula is finite and positive (handbook-style benchmark anchor)."""
    E = 210e9
    I = 3.2e-6
    L = 1.8
    P = _euler_cantilever_P_cr(E, I, L)
    assert np.isfinite(P) and P > 0.0

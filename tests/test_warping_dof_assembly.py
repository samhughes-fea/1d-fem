"""
Tests for 7-DOF warping beam elements: DOF assembly and global index mapping.

Asserts that warping elements return 14 global DOF indices (7 per node) and that
total_dof is set from max dof_per_node when using the static runner.
"""

import os
import sys
from pathlib import Path

import numpy as np
import pytest

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from pre_processing.element_library.linear.beam.first_order_shear_deformation_theory.timoshenko.linear_timoshenko_3D import (
    LinearTimoshenkoBeamElement3D,
)


def _minimal_section_dict_with_gamma():
    """Section dict with 10 entries (A, I_x, I_y, I_z, J_t, kappa, alpha, y_sc, z_sc, Gamma)."""
    return {
        "A": np.array([0.001]),
        "I_x": np.array([0.0]),
        "I_y": np.array([1e-8]),
        "I_z": np.array([1e-8]),
        "J_t": np.array([1e-9]),
        "kappa": np.array([5.0 / 6.0]),
        "alpha": np.array([0.0]),
        "y_sc": np.array([0.0]),
        "z_sc": np.array([0.0]),
        "Gamma": np.array([1.0e-8]),
    }


def test_warping_dof_assembly():
    """7-DOF warping element returns 14 global DOF indices in element order [n1 0..5, n2 0..5, chi1, chi2]."""
    L = 1.0
    element_dictionary = {
        "ids": np.array([0]),
        "connectivity": np.array([[0, 1]]),
        "types": np.array(["LinearTimoshenkoBeamElement3D"]),
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
    grid_dictionary = {"coordinates": np.array([[0.0, 0.0, 0.0], [L, 0.0, 0.0]])}
    material_dictionary = {
        "E": np.array([2.1e11]),
        "G": np.array([8.1e10]),
        "nu": np.array([0.3]),
        "rho": np.array([7850.0]),
    }
    section_dictionary = _minimal_section_dict_with_gamma()
    point_load_array = np.empty((0, 9))
    distributed_load_array = np.empty((0, 9))

    import tempfile
    temp_dir = tempfile.mkdtemp()
    job_results_dir = os.path.join(temp_dir, "warping_test")
    os.makedirs(job_results_dir, exist_ok=True)
    os.makedirs(os.path.join(job_results_dir, "element_stiffness_matrices"), exist_ok=True)
    os.makedirs(os.path.join(job_results_dir, "element_force_vectors"), exist_ok=True)

    try:
        element = LinearTimoshenkoBeamElement3D(
            element_id=0,
            element_dictionary=element_dictionary,
            grid_dictionary=grid_dictionary,
            material_dictionary=material_dictionary,
            section_dictionary=section_dictionary,
            point_load_array=point_load_array,
            distributed_load_array=distributed_load_array,
            job_results_dir=job_results_dir,
        )
        dof_indices = element.assemble_global_dof_indices()
        assert len(dof_indices) == 14, f"Expected 14 DOF indices, got {len(dof_indices)}"
        assert element.dof_per_node == 7
        # element order: six standard DOFs of each node, then the two warping DOFs
        expected = np.concatenate([np.arange(0, 6), np.arange(7, 13), [6, 13]])
        np.testing.assert_array_equal(dof_indices, expected)
    finally:
        import shutil
        if os.path.isdir(temp_dir):
            shutil.rmtree(temp_dir, ignore_errors=True)


def test_process_job_warping_timoshenko_cantilever_is_assembled_in_element_order(tmp_path):
    """A 7-DOF mesh gives the same tip response as the closed form (regression for the DOF map).

    Element matrices are built as [12 standard, chi1, chi2]; mapping them node-major
    (chi1 at local 6) mixed u_x of node 2 with chi of node 1.
    """
    import csv

    from workflow_orchestrator.run_job import process_job, setup_job_results_directory

    L, n_el, E, A, Iz = 2.0, 4, 2.1e11, 1e-3, 2e-6
    Fx, Mx, G, Jt = 1e4, 50.0, 8.1e10, 1e-7
    job = tmp_path / "job_ts_warp_dof_map"
    job.mkdir()
    nodes = "\n".join(f"{i} {L * i / n_el:.6f} 0.0 0.0" for i in range(n_el + 1))
    (job / "grid.txt").write_text(f"[Grid]\n[node_id] [x] [y] [z]\n{nodes}\n", encoding="utf-8")
    rows = "\n".join(f"{e} {e} {e + 1} LinearTimoshenkoBeamElement3D 2 2 2 2 2 2 2 1" for e in range(n_el))
    (job / "element.txt").write_text(
        "[Element]\n[element_id] [node1] [node2] [element_type] [axial_order] [bending_y_order] "
        f"[bending_z_order] [shear_y_order] [shear_z_order] [torsion_order] [load_order] [warping]\n{rows}\n",
        encoding="utf-8",
    )
    (job / "material.txt").write_text(
        "[Material]\n[element_id] [E] [G] [nu] [rho]\n" + "\n".join(f"{e} {E} {G} 0.3 0" for e in range(n_el)) + "\n",
        encoding="utf-8",
    )
    (job / "section.txt").write_text(
        "[Section]\n[element_id] [A] [I_x] [I_y] [I_z] [J_t] [kappa] [alpha] [y_sc] [z_sc] [Gamma]\n"
        + "\n".join(f"{e} {A} 0.0 1e-6 {Iz} {Jt} 0.8333 0.0 0.0 0.0 0.0" for e in range(n_el)) + "\n",
        encoding="utf-8",
    )
    dofs = ["UX", "UY", "UZ", "RX", "RY", "RZ", "W"]
    (job / "prescribed_displacement.txt").write_text(
        "[Prescribed Displacement]\n[id] [node_id] [dof] [value] [type] [comment]\n"
        + "\n".join(f"{i} 0 {d} 0.0 displacement" for i, d in enumerate(dofs)) + "\n",
        encoding="utf-8",
    )
    (job / "point_load.txt").write_text(
        "[Point load]\n[x] [y] [z] [F_x] [F_y] [F_z] [M_x] [M_y] [M_z]\n"
        f"{L} 0.0 0.0 {Fx} 0.0 0.0 {Mx} 0.0 0.0\n",
        encoding="utf-8",
    )
    (job / "simulation_settings.txt").write_text(
        "[Simulation]\n[Type]\nstatic\n\n[Solver]\ntype = direct\ntolerance = 1e-12\nmax_iterations = 500\n",
        encoding="utf-8",
    )
    res = setup_job_results_directory("test_ts_warp_dof_map")
    process_job(str(job), res, {}, {}, force_serial=True, max_processes_per_job=1)
    with open(os.path.join(res, "primary_results", "global", "U_global.csv"), newline="", encoding="utf-8") as f:
        rows = list(csv.reader(f))
    U = np.array([float(r[1]) for r in rows[1:]])
    assert U.size == 7 * (n_el + 1)
    tip = U[7 * n_el:7 * n_el + 7]
    assert tip[0] == pytest.approx(Fx * L / (E * A), rel=1e-9)
    # Gamma = 0: free warping, Saint-Venant twist
    assert tip[3] == pytest.approx(Mx * L / (G * Jt), rel=1e-9)
    assert np.all(np.abs(tip[[1, 2, 4, 5]]) < 1e-12)

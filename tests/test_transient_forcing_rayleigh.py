"""Regression: transient ``F(t)`` from file/analytic and Rayleigh ``C`` when element damping absent."""

from __future__ import annotations

import shutil
import sys
from pathlib import Path

import numpy as np
from scipy.sparse import coo_matrix

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from processing.dynamic.transient_forcing import build_transient_force_func
from simulation_runner.eigen.eigen_simulation import EigenSimulationRunner
from simulation_runner.transient.dynamic_simulation import TransientSimulationRunner
from tests.test_modal_buckling_euler_column import _build_cantilever_modal_case


def _element_matrix_to_coo(m):
    if hasattr(m, "tocoo"):
        return m.tocoo()
    return coo_matrix(np.asarray(m, dtype=np.float64))


def test_transient_force_time_series_scale():
    settings, job_tmp, _ = _build_cantilever_modal_case(4, 1.0, 210e9, 1.0)
    try:
        ts = Path(job_tmp) / "fscale.txt"
        ts.write_text("0.0 0.5\n0.05 1.0\n", encoding="utf-8")
        eos = list(np.asarray(settings["element_objects"], dtype=object).ravel())
        elements = list(np.asarray(settings["elements"], dtype=object).ravel())
        element_stiffness_matrices_dyn = np.asarray(
            [_element_matrix_to_coo(o.K_e) for o in eos], dtype=object
        )
        mass_objs = [e.element_mass_matrix() for e in elements]
        element_mass_matrices_dyn = np.asarray(
            [_element_matrix_to_coo(mo.M_e) for mo in mass_objs], dtype=object
        )
        dyn_settings = {
            "elements": settings["elements"],
            "mesh_dictionary": settings["mesh_dictionary"],
            "grid_dictionary": settings["grid_dictionary"],
            "element_dictionary": settings["element_dictionary"],
            "material_dictionary": settings["material_dictionary"],
            "section_dictionary": settings["section_dictionary"],
            "point_load_array": settings["point_load_array"],
            "distributed_load_array": settings["distributed_load_array"],
            "element_stiffness_matrices": element_stiffness_matrices_dyn,
            "element_mass_matrices": element_mass_matrices_dyn,
            "element_objects": settings["element_objects"],
            "force_objects": settings["force_objects"],
            "job_results_dir": settings["job_results_dir"],
            "job_dir": str(job_tmp),
            "simulation_settings": {
                "transient": {
                    "time_step": 0.01,
                    "end_time": 0.04,
                    "load_scale": 1.0,
                    "force_time_series_file": "fscale.txt",
                }
            },
        }
        TransientSimulationRunner(settings=dyn_settings, job_name="ft_series").run()
        root = Path(settings["job_results_dir"])
        u = np.loadtxt(root / "primary_results" / "dynamic_results" / "ft_series_displacements.txt")
        assert u.shape[0] >= 2
    finally:
        shutil.rmtree(job_tmp, ignore_errors=True)


def test_transient_rayleigh_adds_damping_matrix():
    settings, job_tmp, _ = _build_cantilever_modal_case(4, 1.0, 210e9, 1.0)
    try:
        eos = list(np.asarray(settings["element_objects"], dtype=object).ravel())
        elements = list(np.asarray(settings["elements"], dtype=object).ravel())
        element_stiffness_matrices_dyn = np.asarray(
            [_element_matrix_to_coo(o.K_e) for o in eos], dtype=object
        )
        mass_objs = [e.element_mass_matrix() for e in elements]
        element_mass_matrices_dyn = np.asarray(
            [_element_matrix_to_coo(mo.M_e) for mo in mass_objs], dtype=object
        )
        dyn_settings = {
            "elements": settings["elements"],
            "mesh_dictionary": settings["mesh_dictionary"],
            "grid_dictionary": settings["grid_dictionary"],
            "element_dictionary": settings["element_dictionary"],
            "material_dictionary": settings["material_dictionary"],
            "section_dictionary": settings["section_dictionary"],
            "point_load_array": settings["point_load_array"],
            "distributed_load_array": settings["distributed_load_array"],
            "element_stiffness_matrices": element_stiffness_matrices_dyn,
            "element_mass_matrices": element_mass_matrices_dyn,
            "element_damping_matrices": None,
            "element_objects": settings["element_objects"],
            "force_objects": settings["force_objects"],
            "job_results_dir": settings["job_results_dir"],
            "simulation_settings": {
                "transient": {
                    "time_step": 0.01,
                    "end_time": 0.03,
                    "rayleigh_alpha": 0.02,
                    "rayleigh_beta": 1.0e-6,
                }
            },
        }
        TransientSimulationRunner(settings=dyn_settings, job_name="rayleigh_dyn").run()
        root = Path(settings["job_results_dir"])
        assert (root / "primary_results" / "dynamic_results" / "rayleigh_dyn_displacements.txt").is_file()
    finally:
        shutil.rmtree(job_tmp, ignore_errors=True)


def test_eigen_secondary_generalized_mass_file():
    settings, job_tmp, _ = _build_cantilever_modal_case(6, 1.0, 210e9, 1.0)
    settings["simulation_settings"] = {"eigen": {"num_modes": 3}}
    try:
        EigenSimulationRunner(settings=settings, job_name="eigen_gm").run()
        root_modal = Path(settings["job_results_dir"]) / "primary_results" / "modal_results"
        p = root_modal / "eigen_gm_modal_generalized_mass.txt"
        assert p.is_file()
        gm = np.loadtxt(str(p))
        assert gm.shape == (3,)
        assert np.all(gm > 0.0)
        pp = root_modal / "eigen_gm_modal_load_participation.txt"
        assert pp.is_file()
        part = np.loadtxt(str(pp))
        assert part.shape == (3,)
    finally:
        shutil.rmtree(job_tmp, ignore_errors=True)


def test_build_transient_force_func_sin_burst():
    F_ref = np.ones(4, dtype=np.float64)
    cfg = {
        "load_scale": 2.0,
        "load_ramp": False,
        "end_time": 1.0,
        "force_analytic": "sin_burst",
        "force_analytic_amplitude": 0.5,
        "force_analytic_frequency_hz": 2.0,
        "force_analytic_phase_rad": 0.0,
        "force_analytic_t_start": 0.1,
        "force_analytic_t_end": 0.2,
    }
    fn = build_transient_force_func(
        F_ref, cfg, total_dof=4, job_dir=None, end_time=1.0
    )
    f_out = fn(0.125)
    assert f_out.shape == (4,)
    # Inside burst: 2*pi*2*0.125 = pi/2 => sin = 1 => envelope 1.5 * base 2.0
    np.testing.assert_allclose(f_out, 3.0 * F_ref, rtol=1e-9)
    f_flat = fn(0.0)
    np.testing.assert_allclose(f_flat, 2.0 * F_ref, rtol=1e-9)

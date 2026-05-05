"""Regression: transient ``F(t)`` from file/analytic and Rayleigh ``C`` when element damping absent."""

from __future__ import annotations

import shutil
import sys
from pathlib import Path

import numpy as np
import pytest
from scipy.sparse import coo_matrix

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from processing.transient.transient_forcing import build_transient_force_func
from simulation_runner.eigen.eigen_simulation import EigenSimulationRunner
from simulation_runner.transient.dynamic_simulation import TransientSimulationRunner
from tests.test_modal_buckling_euler_column import _build_cantilever_modal_case
from workflow_orchestrator.run_job import process_job, setup_job_results_directory


def _pinned_eigen_cantilever_reference_frequencies_hz(n: int) -> np.ndarray:
    return np.array([11.253952, 28.593523, 72.731488, 184.790597], dtype=np.float64)[:n]


def _pinned_linear_buckling_reference_load_factors(n: int) -> np.ndarray:
    return np.array([23.526405], dtype=np.float64)[:n]


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
        summary = root / "primary_results" / "dynamic_results" / "rayleigh_dyn_primary_summary.csv"
        assert summary.is_file()
        text = summary.read_text(encoding="utf-8")
        assert "damping_source,rayleigh" in text
        assert "max_abs_displacement,9.517172091466624e-10" in text
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


def test_pinned_eigen_benchmark_acceptance_artifacts():
    job_dir = PROJECT_ROOT / "jobs" / "job_smoke_eigen"
    res_dir = setup_job_results_directory("pytest_eigen_benchmark_acceptance")
    process_job(
        str(job_dir),
        res_dir,
        {},
        {},
        force_serial=True,
        max_processes_per_job=1,
    )
    root = Path(res_dir)
    modal_dir = root / "primary_results" / "modal_results"
    f = modal_dir / "job_smoke_eigen_frequencies.txt"
    m = modal_dir / "job_smoke_eigen_mode_shapes.txt"
    manifest = root / "logs" / "primary_artifacts.json"
    assert f.is_file()
    assert m.is_file()
    assert manifest.is_file()
    freqs = np.loadtxt(f)
    arr = np.atleast_1d(np.asarray(freqs, dtype=float))
    assert np.all(arr > 0.0)
    f_ref = _pinned_eigen_cantilever_reference_frequencies_hz(min(4, arr.shape[0]))
    np.testing.assert_allclose(arr[: f_ref.shape[0]], f_ref, rtol=2.0e-2, atol=1.0e-6)


def test_pinned_linear_buckling_benchmark_acceptance_artifacts():
    job_dir = PROJECT_ROOT / "jobs" / "job_smoke_buckling"
    res_dir = setup_job_results_directory("pytest_linear_buckling_benchmark_acceptance")
    process_job(
        str(job_dir),
        res_dir,
        {},
        {},
        force_serial=True,
        max_processes_per_job=1,
    )
    root = Path(res_dir)
    modal_dir = root / "primary_results" / "modal_results"
    f = modal_dir / "job_smoke_buckling_buckling_load_factors.txt"
    m = modal_dir / "job_smoke_buckling_buckling_mode_shapes.txt"
    manifest = root / "logs" / "primary_artifacts.json"
    assert f.is_file()
    assert m.is_file()
    assert manifest.is_file()
    lambdas = np.loadtxt(f)
    arr = np.atleast_1d(np.asarray(lambdas, dtype=float))
    assert np.all(np.isfinite(arr))
    assert np.all(arr > 0.0)
    lam_ref = _pinned_linear_buckling_reference_load_factors(min(1, arr.shape[0]))
    np.testing.assert_allclose(arr[: lam_ref.shape[0]], lam_ref, rtol=5.0e-2, atol=1.0e-8)


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


def test_transient_sdof_newmark_matches_analytical_cosine_response():
    from scipy.sparse import csr_matrix

    m = csr_matrix([[1.0]])
    k = csr_matrix([[100.0]])
    c = None
    t_grid = np.linspace(0.0, 0.5, 101)
    u0 = np.array([1.0], dtype=np.float64)
    v0 = np.array([0.0], dtype=np.float64)

    def force(_t: float) -> np.ndarray:
        return np.array([0.0], dtype=np.float64)

    from processing.transient.time_integration import newmark_integrate

    U, _V, _A = newmark_integrate(k, m, c, u0, v0, t_grid, force)
    omega_n = np.sqrt(100.0)
    u_ref = np.cos(omega_n * t_grid)
    err = np.max(np.abs(U[:, 0] - u_ref))
    assert err <= 2.0e-2


@pytest.mark.integration
def test_transient_multidof_cantilever_benchmark_artifacts_and_tip_history(tmp_path):
    settings, job_tmp, _ = _build_cantilever_modal_case(4, 2.0, 210e9, 1.0)
    try:
        eos = list(np.asarray(settings["element_objects"], dtype=object).ravel())
        elements = np.asarray(settings["elements"], dtype=object)
        element_stiffness_matrices_dyn = np.asarray(
            [_element_matrix_to_coo(o.K_e) for o in eos], dtype=object
        )
        mass_objs = [e.element_mass_matrix() for e in list(elements.ravel())]
        element_mass_matrices_dyn = np.asarray(
            [_element_matrix_to_coo(mo.M_e) for mo in mass_objs], dtype=object
        )
        dyn_settings = {
            "elements": elements,
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
            "simulation_settings": {
                "transient": {
                    "time_step": 0.01,
                    "end_time": 0.05,
                    "force_analytic": "sin",
                    "force_analytic_amplitude": 1.0,
                    "force_analytic_frequency_hz": 5.0,
                }
            },
        }
        TransientSimulationRunner(settings=dyn_settings, job_name="tdof_benchmark").run()
        root = Path(settings["job_results_dir"])
        out = root / "primary_results" / "dynamic_results"
        u = np.loadtxt(out / "tdof_benchmark_displacements.txt")
        assert u.ndim == 2
        assert u.shape[0] >= 2
        tip_hist = u[:, -6]
        assert np.max(np.abs(tip_hist)) >= 0.0
        summary = out / "tdof_benchmark_primary_summary.csv"
        assert summary.is_file()
    finally:
        shutil.rmtree(job_tmp, ignore_errors=True)

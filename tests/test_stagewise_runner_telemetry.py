"""Assert RuntimeMonitorTelemetry stage markers for stagewise runners."""

from __future__ import annotations

import shutil
import sys
from pathlib import Path

import numpy as np
import pytest
from scipy.sparse import coo_matrix

PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from tests.test_modal_buckling_euler_column import _build_cantilever_modal_case
from workflow_orchestrator.run_job import process_job, setup_job_results_directory


def _element_matrix_to_coo(m):
    if hasattr(m, "tocoo"):
        return m.tocoo()
    return coo_matrix(np.asarray(m, dtype=np.float64))


def test_dynamic_runner_runtime_monitor_telemetry_stages():
    from simulation_runner.transient.dynamic_simulation import DynamicSimulationRunner

    settings, tmp, _ = _build_cantilever_modal_case(8, 2.5, 210e9, 1.0)
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
        "simulation_settings": {"dynamic": {"time_step": 0.01, "end_time": 0.03}},
    }
    try:
        DynamicSimulationRunner(settings=dyn_settings, job_name="telemetry_dyn").run()
        root = Path(settings["job_results_dir"])
        log_path = root / "diagnostics" / "RuntimeMonitorTelemetry.log"
        assert log_path.is_file(), f"missing telemetry log: {log_path}"
        text = log_path.read_text(encoding="utf-8", errors="replace")
        assert "AssembleDynamicGlobalSystem BEGIN" in text
        assert "AssembleDynamicGlobalSystem END" in text
        assert "ModifyDynamicGlobalSystem BEGIN" in text
        assert "IntegrateTransientSystem END" in text
        tdiag = root / "logs" / "transient_run_diagnostic.log"
        assert tdiag.is_file(), f"missing transient diagnostic log: {tdiag}"
        integ = root / "logs" / "IntegrateTransientSystem.log"
        assert integ.is_file()
        assert "Transient stability snapshot" in integ.read_text(encoding="utf-8", errors="replace")
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def test_dynamic_runner_stage_log_files_under_logs():
    from simulation_runner.transient.dynamic_simulation import DynamicSimulationRunner

    settings, tmp, _ = _build_cantilever_modal_case(8, 2.5, 210e9, 1.0)
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
        "simulation_settings": {"dynamic": {"time_step": 0.01, "end_time": 0.02}},
    }
    try:
        DynamicSimulationRunner(settings=dyn_settings, job_name="stage_logs_dyn").run()
        root = Path(settings["job_results_dir"])
        logs = root / "logs"
        assert (logs / "AssembleDynamicGlobalSystem.log").is_file()
        assert (logs / "ModifyDynamicGlobalSystem.log").is_file()
        assert (logs / "IntegrateTransientSystem.log").is_file()
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


@pytest.mark.integration
def test_eigen_process_job_runtime_monitor_telemetry_stages():
    job_dir = PROJECT_ROOT / "jobs" / "job_smoke_eigen"
    if not job_dir.is_dir():
        pytest.skip(f"missing job directory: {job_dir}")
    res_dir = setup_job_results_directory("pytest_eigen_telemetry_stages")
    jt: dict = {}
    je: dict = {}
    try:
        process_job(
            str(job_dir),
            res_dir,
            jt,
            je,
            force_serial=True,
            max_processes_per_job=1,
        )
        log_path = Path(res_dir) / "diagnostics" / "RuntimeMonitorTelemetry.log"
        assert log_path.is_file(), f"missing telemetry log: {log_path}"
        text = log_path.read_text(encoding="utf-8", errors="replace")
        assert "PrepareSpectralLocalMatrices BEGIN" in text
        assert "AssembleSpectralGlobalSystem END" in text
        assert "SolveGeneralizedEigenproblem" in text
        logs = Path(res_dir) / "logs"
        assert (logs / "PrepareSpectralLocalMatrices.log").is_file()
        assert (logs / "ModifySpectralGlobalSystem.log").is_file()
    finally:
        shutil.rmtree(res_dir, ignore_errors=True)


@pytest.mark.integration
def test_buckling_process_job_runtime_monitor_telemetry_stages():
    job_dir = PROJECT_ROOT / "jobs" / "job_smoke_buckling"
    if not job_dir.is_dir():
        pytest.skip(f"missing job directory: {job_dir}")
    res_dir = setup_job_results_directory("pytest_buckling_telemetry_stages")
    jt: dict = {}
    je: dict = {}
    try:
        process_job(
            str(job_dir),
            res_dir,
            jt,
            je,
            force_serial=True,
            max_processes_per_job=1,
        )
        log_path = Path(res_dir) / "diagnostics" / "RuntimeMonitorTelemetry.log"
        assert log_path.is_file(), f"missing telemetry log: {log_path}"
        text = log_path.read_text(encoding="utf-8", errors="replace")
        assert "BucklingPrestress BEGIN" in text
        assert "LinearStaticPrestress BEGIN" in text
        assert "AssembleBucklingGeometricStiffness" in text
    finally:
        shutil.rmtree(res_dir, ignore_errors=True)

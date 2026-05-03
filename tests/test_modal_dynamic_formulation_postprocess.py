"""
Integration: modal vibration and dynamic runners optional formulation-cache secondary/tertiary path.

Runs with FEM_FORMULATION_CACHE_STRICT_SHAPE enforced via pytest monkeypatch.
"""

from __future__ import annotations

import shutil
import sys
from pathlib import Path

import numpy as np
import pytest
from scipy.sparse import coo_matrix

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from tests.test_modal_buckling_euler_column import _build_cantilever_modal_case


def _element_matrix_to_coo(m):
    if hasattr(m, "tocoo"):
        return m.tocoo()
    return coo_matrix(np.asarray(m, dtype=np.float64))


@pytest.fixture
def strict_shape_env(monkeypatch):
    monkeypatch.setenv("FEM_FORMULATION_CACHE_STRICT_SHAPE", "1")


def test_modal_buckling_prestress_secondary_tertiary(strict_shape_env):
    """Modal buckling branch + prestress snapshot U (avoids vibration ARPACK on tiny meshes)."""
    from simulation_runner.buckling.buckling_simulation import BucklingSimulationRunner

    settings, tmp, _ = _build_cantilever_modal_case(8, 2.5, 210e9, 1.0)
    settings["simulation_settings"] = {
        "modal": {
            "num_modes": 2,
            "analysis": "buckling",
            "buckling_prestress": "linear_static",
            "buckling_load_factor": 1.0,
        },
        "post_processing": {
            "run_secondary_tertiary_modal": True,
            "buckling_displacement": "prestress",
            "modal_amplitude": 1.0,
        },
    }
    try:
        BucklingSimulationRunner(settings=settings, job_name="modal_pp_test").run()
        root = Path(settings["job_results_dir"])
        assert (root / "secondary_results" / "secondary_summary.csv").is_file()
        assert (root / "tertiary_results" / "tertiary_summary.csv").is_file()
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def test_dynamic_secondary_tertiary(strict_shape_env):
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
        "simulation_settings": {
            "dynamic": {"time_step": 0.01, "end_time": 0.05},
            "post_processing": {
                "run_secondary_tertiary_dynamic": True,
                "dynamic_time_index": -1,
            },
        },
    }
    try:
        DynamicSimulationRunner(settings=dyn_settings, job_name="dynamic_pp_test").run()
        root = Path(settings["job_results_dir"])
        assert (root / "secondary_results" / "secondary_summary.csv").is_file()
        assert (root / "tertiary_results" / "tertiary_summary.csv").is_file()
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def test_parse_post_processing_section(tmp_path):
    from pre_processing.parsing.simulation_settings_parser import parse_simulation_settings

    p = tmp_path / "simulation_settings.txt"
    p.write_text(
        "\n".join(
            [
                "[Simulation]",
                "[Type]",
                "modal",
                "[PostProcessing]",
                "run_secondary_tertiary_modal = true",
                "modal_mode_index = 1",
                "buckling_displacement = prestress",
                "dynamic_time_index = -2",
            ]
        ),
        encoding="utf-8",
    )
    cfg = parse_simulation_settings(str(p))
    assert cfg["type"] == "eigen"
    pp = cfg["post_processing"]
    assert pp["run_secondary_tertiary_modal"] is True
    assert pp["modal_mode_index"] == 1
    assert pp["buckling_displacement"] == "prestress"
    assert pp["dynamic_time_index"] == -2

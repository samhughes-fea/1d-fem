"""Lock modular runner orchestration: staged ``self.*`` calls must appear in pipeline order."""

from __future__ import annotations

import ast
import sys
from pathlib import Path

import pytest

PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))


def _function_source_lines(path: Path, class_name: str, method_name: str) -> str:
    tree = ast.parse(path.read_text(encoding="utf-8"))
    for node in tree.body:
        if isinstance(node, ast.ClassDef) and node.name == class_name:
            for item in node.body:
                if isinstance(item, ast.FunctionDef) and item.name == method_name:
                    end = getattr(item, "end_lineno", None)
                    if end is None:
                        pytest.fail("AST FunctionDef missing end_lineno; use Python 3.8+")
                    lines = path.read_text(encoding="utf-8").splitlines()
                    return "\n".join(lines[item.lineno - 1 : end])
    raise AssertionError(f"{class_name}.{method_name} not found in {path}")


def _assert_needles_monotonic(source: str, needles: list[str]) -> None:
    last = -1
    for needle in needles:
        pos = source.find(needle)
        assert pos != -1, f"expected orchestration to reference {needle!r}"
        assert pos > last, (
            f"expected {needle!r} after prior stage (got pos {pos} <= last {last}); "
            "staged runner order may have regressed"
        )
        last = pos


def test_harmonic_runner_run_self_stage_order() -> None:
    path = PROJECT_ROOT / "simulation_runner" / "harmonic" / "harmonic_simulation.py"
    src = _function_source_lines(path, "HarmonicSimulationRunner", "run")
    _assert_needles_monotonic(
        src,
        [
            "self.prepare_harmonic_job_tree()",
            "self.collect_harmonic_assembly_inputs()",
            "self.assemble_harmonic_structural_matrices(",
            "self.assemble_harmonic_load_vector(",
            "self.modify_harmonic_structural_matrices(",
            "self.build_harmonic_damping_matrix(",
            "self.log_harmonic_structural_diagnostics(",
            "self.build_harmonic_complex_load_vector(",
            "self._harmonic_prescribed_partition(",
            "self.solve_harmonic_frequency_sweep(",
            "self.save_harmonic_primary_results(",
            "self.write_harmonic_primary_artifact_manifest()",
            "self.run_harmonic_post_processing_if_enabled(",
        ],
    )


def test_dynamic_runner_run_self_stage_order() -> None:
    path = PROJECT_ROOT / "simulation_runner" / "transient" / "dynamic_simulation.py"
    src = _function_source_lines(path, "TransientSimulationRunner", "run")
    _assert_needles_monotonic(
        src,
        [
            "self.setup_simulation()",
            "self.assemble_dynamic_global_system(",
            "self.apply_rayleigh_damping_if_needed(",
            "self.modify_dynamic_global_system(",
            "self.integrate_transient_system(",
            "self.save_transient_primary_results(",
            "self.run_dynamic_post_processing_if_enabled(",
            "self.write_transient_primary_artifact_manifest()",
        ],
    )


def test_vibration_buckling_backend_vibration_pipeline_order() -> None:
    path = PROJECT_ROOT / "simulation_runner" / "spectral" / "vibration_buckling_backend.py"
    src = _function_source_lines(path, "VibrationBucklingBackend", "_run_vibration_analysis")
    _assert_needles_monotonic(
        src,
        [
            "self.prepare_spectral_local_matrices(",
            "self.assemble_spectral_global_system(",
            "self.modify_spectral_global_system(",
            "self.solve_vibration_eigenpairs(",
            "self.save_vibration_primary_results(",
            "self.finalize_vibration_secondary_or_post(",
        ],
    )


def test_vibration_buckling_backend_buckling_pipeline_order() -> None:
    path = PROJECT_ROOT / "simulation_runner" / "spectral" / "vibration_buckling_backend.py"
    src = _function_source_lines(path, "VibrationBucklingBackend", "_run_buckling_analysis")
    _assert_needles_monotonic(
        src,
        [
            "self.prepare_spectral_local_matrices(",
            "self.assemble_spectral_global_system(",
            "self.solve_buckling_linearized_problem(",
            "self.save_buckling_primary_results(",
            "self.finalize_buckling_optional_post(",
        ],
    )

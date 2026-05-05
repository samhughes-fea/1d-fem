"""
Tests for post_processing validation_visualisers (FEM vs Abaqus).

- job_to_abaqus_script: generated script exists and contains expected Abaqus API strings.
- Comparison scripts: run without error when abaqus_results is empty or missing (no Abaqus required).
"""

import sys
from pathlib import Path

import pytest

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))


def test_job_to_abaqus_script_generates_file_with_expected_content():
    """Generate Abaqus script for job_0000_n8; assert file exists and contains WirePolyLine, BeamSection, job name."""
    from post_processing.validation_visualisers.abaqus.job_to_abaqus_script import (
        _parse_job,
        _generate_script_content,
    )
    from post_processing.validation_visualisers.abaqus.config import (
        ABAQUS_GENERATED_DIR,
        ABAQUS_RESULTS_DIR,
    )

    job_dir = PROJECT_ROOT / "jobs" / "job_0000_n8"
    if not job_dir.is_dir():
        pytest.skip("jobs/job_0000_n8 not found")
    data = _parse_job(job_dir)
    out_csv_dir = str(ABAQUS_RESULTS_DIR / data["job_name"])
    content = _generate_script_content(data, out_csv_dir)
    assert "WirePolyLine" in content
    assert "BeamSection" in content
    assert "job_0000_n8" in content or "job_0000_n8_abaqus" in content
    assert "OUT_CSV_DIR" in content
    assert "U_global.csv" in content
    # abqpy: driverUtils, executeOnCaeStartup, IN_ABAQUS, saveAs, exit after saveAs
    assert "driverUtils" in content
    assert "executeOnCaeStartup" in content
    assert "IN_ABAQUS" in content
    assert "mdb.saveAs(" in content
    assert "sys.exit(0)" in content


def test_job_to_abaqus_script_parses_transient_benchmark_contract() -> None:
    from post_processing.validation_visualisers.abaqus.job_to_abaqus_script import _parse_job

    job_dir = PROJECT_ROOT / "jobs" / "job_benchmark_transient_cantilever_multidof"
    if not job_dir.is_dir():
        pytest.skip("transient benchmark job not found")
    data = _parse_job(job_dir)
    assert data["job_name"] == "job_benchmark_transient_cantilever_multidof"
    assert data["simulation_settings_path"].endswith("simulation_settings.txt")
    assert data["simulation_type"] == "transient"


def test_job_to_abaqus_script_emits_transient_step_metadata() -> None:
    from post_processing.validation_visualisers.abaqus.job_to_abaqus_script import _parse_job, _generate_script_content
    from post_processing.validation_visualisers.abaqus.config import ABAQUS_RESULTS_DIR

    job_dir = PROJECT_ROOT / "jobs" / "job_benchmark_transient_cantilever_multidof"
    if not job_dir.is_dir():
        pytest.skip("transient benchmark job not found")
    data = _parse_job(job_dir)
    out_csv_dir = str(ABAQUS_RESULTS_DIR / data["job_name"])
    content = _generate_script_content(data, out_csv_dir)
    assert "ImplicitDynamicsStep" in content
    assert "transient_reference_contract.txt" in content


def test_job_to_abaqus_script_emits_nonlinear_static_tip_history_contract() -> None:
    from post_processing.validation_visualisers.abaqus.job_to_abaqus_script import _parse_job, _generate_script_content
    from post_processing.validation_visualisers.abaqus.config import ABAQUS_RESULTS_DIR

    job_dir = PROJECT_ROOT / "jobs" / "job_benchmark_nl_static_cantilever_tip"
    if not job_dir.is_dir():
        pytest.skip("nonlinear static benchmark job not found")
    data = _parse_job(job_dir)
    out_csv_dir = str(ABAQUS_RESULTS_DIR / data["job_name"])
    content = _generate_script_content(data, out_csv_dir)
    assert "export_tip_history=True" in content
    assert "extract_odb_to_csv" in content


def test_tip_load_fine_reference_pair_docs_exist() -> None:
    plan = PROJECT_ROOT / "docs" / "conventions" / "NONLINEAR_STATIC_ABAQUS_VALIDATION_PLAN.md"
    text = plan.read_text(encoding="utf-8")
    assert "job_benchmark_nl_static_cantilever_tip_n64" in text
    assert "job_benchmark_nl_static_cantilever_tip_n500" in text


def test_wave1_pyabaqus_contract_docs_exist() -> None:
    docs = {
        "LINEAR_STATIC_PYABAQUS_REFERENCE_CONTRACT.md": [
            "Linear static pyAbaqus reference contract",
            "U_global.csv",
            "rotation_source.txt",
        ],
        "HARMONIC_PYABAQUS_REFERENCE_CONTRACT.md": [
            "Harmonic pyAbaqus reference contract",
            "frequency_response.csv",
            "frequency_hz",
        ],
        "EIGEN_PYABAQUS_REFERENCE_CONTRACT.md": [
            "Eigen pyAbaqus reference contract",
            "eigen_frequencies.csv",
            "mode_shapes.csv",
        ],
        "LINEAR_BUCKLING_PYABAQUS_REFERENCE_CONTRACT.md": [
            "Linear buckling pyAbaqus reference contract",
            "buckling_load_factors.csv",
            "buckling_mode_shapes.csv",
        ],
    }
    base = PROJECT_ROOT / "docs" / "conventions"
    for name, required_snippets in docs.items():
        path = base / name
        assert path.is_file(), f"missing contract doc: {name}"
        text = path.read_text(encoding="utf-8")
        for snippet in required_snippets:
            assert snippet in text, f"expected snippet {snippet!r} in {name}"


def test_simulation_type_dispatch_payload_reads_canonical_type() -> None:
    from post_processing.validation_visualisers.abaqus.simulation_type_dispatch import build_validation_dispatch_payload

    job_dir = PROJECT_ROOT / "jobs" / "job_smoke_eigen"
    payload = build_validation_dispatch_payload(str(job_dir))
    assert payload["simulation_type"] == "eigen"
    assert isinstance(payload["simulation_settings"], dict)
    assert payload["artifact_contract"]["contract_name"] == "eigen_reference"
    assert "eigen_frequencies.csv" in payload["artifact_contract"]["expected_files"]


def test_simulation_type_dispatch_wave1_contracts_cover_core_families() -> None:
    from post_processing.validation_visualisers.abaqus.simulation_type_dispatch import build_validation_dispatch_payload

    cases = {
        "job_0000_n8": ("static_reference", "U_global.csv"),
        "job_benchmark_transient_cantilever_multidof": ("transient_reference", "transient_reference_contract.txt"),
        "job_smoke_eigen": ("eigen_reference", "eigen_frequencies.csv"),
        "job_smoke_buckling": ("linear_buckling_reference", "buckling_load_factors.csv"),
    }
    for job_name, (contract_name, expected_file) in cases.items():
        job_dir = PROJECT_ROOT / "jobs" / job_name
        if not job_dir.is_dir():
            pytest.skip(f"{job_name} not found")
        payload = build_validation_dispatch_payload(str(job_dir))
        assert payload["artifact_contract"]["contract_name"] == contract_name
        assert expected_file in payload["artifact_contract"]["expected_files"]


def test_job_to_abaqus_script_emits_wave1_contract_artifact_metadata() -> None:
    from post_processing.validation_visualisers.abaqus.job_to_abaqus_script import _parse_job, _generate_script_content
    from post_processing.validation_visualisers.abaqus.config import ABAQUS_RESULTS_DIR

    cases = {
        "job_0000_n8": ["ARTIFACT_CONTRACT_NAME = \"static_reference\"", "artifact_contract.txt"],
        "job_benchmark_transient_cantilever_multidof": ["ARTIFACT_CONTRACT_NAME = \"transient_reference\"", "transient_reference_contract.txt"],
        "job_smoke_eigen": ["ARTIFACT_CONTRACT_NAME = \"eigen_reference\"", "eigen_frequencies.csv", "mode_shapes.csv"],
        "job_smoke_buckling": ["ARTIFACT_CONTRACT_NAME = \"linear_buckling_reference\"", "buckling_load_factors.csv", "buckling_mode_shapes.csv"],
    }
    for job_name, snippets in cases.items():
        job_dir = PROJECT_ROOT / "jobs" / job_name
        if not job_dir.is_dir():
            pytest.skip(f"{job_name} not found")
        data = _parse_job(job_dir)
        out_csv_dir = str(ABAQUS_RESULTS_DIR / data["job_name"])
        content = _generate_script_content(data, out_csv_dir)
        for snippet in snippets:
            assert snippet in content


def test_job_to_abaqus_script_has_first_helper_slice() -> None:
    p = PROJECT_ROOT / "post_processing" / "validation_visualisers" / "abaqus" / "job_to_abaqus_script.py"
    text = p.read_text(encoding="utf-8")
    assert "def _build_script_preamble(" in text


def test_job_to_abaqus_script_has_second_helper_slice() -> None:
    p = PROJECT_ROOT / "post_processing" / "validation_visualisers" / "abaqus" / "job_to_abaqus_script.py"
    text = p.read_text(encoding="utf-8")
    assert "def _build_step_and_model_block(" in text


def test_job_to_abaqus_script_has_third_helper_slice() -> None:
    p = PROJECT_ROOT / "post_processing" / "validation_visualisers" / "abaqus" / "job_to_abaqus_script.py"
    text = p.read_text(encoding="utf-8")
    assert "def _build_loads_and_job_block(" in text


def test_job_to_abaqus_script_has_results_export_prologue_helper_slice() -> None:
    p = PROJECT_ROOT / "post_processing" / "validation_visualisers" / "abaqus" / "job_to_abaqus_script.py"
    text = p.read_text(encoding="utf-8")
    assert "def _build_results_export_prologue_block(" in text


def test_job_to_abaqus_script_has_results_export_displacement_helper_slice() -> None:
    p = PROJECT_ROOT / "post_processing" / "validation_visualisers" / "abaqus" / "job_to_abaqus_script.py"
    text = p.read_text(encoding="utf-8")
    assert "def _build_results_export_displacement_block(" in text


def test_extract_odb_results_mentions_tip_history_flag() -> None:
    p = PROJECT_ROOT / "post_processing" / "validation_visualisers" / "abaqus" / "extract_odb_results.py"
    text = p.read_text(encoding="utf-8")
    assert "export_tip_history" in text
    assert "tip_load_history.csv" in text


def test_extract_odb_results_mentions_wave1_placeholder_exports() -> None:
    p = PROJECT_ROOT / "post_processing" / "validation_visualisers" / "abaqus" / "extract_odb_results.py"
    text = p.read_text(encoding="utf-8")
    assert "frequency_response.csv" in text
    assert "eigen_frequencies.csv" in text
    assert "buckling_load_factors.csv" in text


def test_extract_odb_results_tip_history_helpers_pick_last_node_u2() -> None:
    from post_processing.validation_visualisers.abaqus.extract_odb_results import (
        _tip_node_label_from_u_field,
        _tip_u2_from_field_values,
    )

    class _Value:
        def __init__(self, node_label, data):
            self.nodeLabel = node_label
            self.data = data

    values = [
        _Value(1, (0.0, -1.0e-4, 0.0)),
        _Value(3, (0.0, -3.0e-4, 0.0)),
        _Value(2, (0.0, -2.0e-4, 0.0)),
    ]

    tip_node = _tip_node_label_from_u_field(values)
    assert tip_node == 3
    assert _tip_u2_from_field_values(values, tip_node) == -3.0e-4


def test_extract_odb_results_tip_history_helpers_tolerate_missing_tip() -> None:
    from post_processing.validation_visualisers.abaqus.extract_odb_results import (
        _tip_node_label_from_u_field,
        _tip_u2_from_field_values,
    )

    assert _tip_node_label_from_u_field([]) is None
    assert _tip_u2_from_field_values([], None) == 0.0


def test_deformation_comparison_runs():
    """Run deformation comparison (FEM vs Abaqus); must not raise when abaqus_results empty."""
    from post_processing.validation_visualisers.deformation.deformation_comparison import (
        run_deformation_comparison,
    )
    run_deformation_comparison()


def test_wave1_family_validation_runner_entrypoints_exist(tmp_path) -> None:
    from post_processing.validation_visualisers.static.run_linear_static_validation import run_linear_static_validation
    from post_processing.validation_visualisers.harmonic.run_harmonic_validation import run_harmonic_validation
    from post_processing.validation_visualisers.eigen.run_eigen_validation import run_eigen_validation
    from post_processing.validation_visualisers.buckling.run_buckling_validation import run_buckling_validation

    runners = [
        run_linear_static_validation,
        run_harmonic_validation,
        run_eigen_validation,
        run_buckling_validation,
    ]
    for runner in runners:
        out = runner(tmp_path)
        assert isinstance(out, dict)
        assert "family" in out
        assert "summary_file" in out


def test_wave1_family_validation_runners_report_reference_artifact_presence(tmp_path) -> None:
    from post_processing.validation_visualisers.harmonic.run_harmonic_validation import run_harmonic_validation
    from post_processing.validation_visualisers.eigen.run_eigen_validation import run_eigen_validation
    from post_processing.validation_visualisers.buckling.run_buckling_validation import run_buckling_validation

    for runner in [run_harmonic_validation, run_eigen_validation, run_buckling_validation]:
        out = runner(tmp_path)
        assert "reference_dir" in out
        assert "present_files" in out
        summary_text = Path(out["summary_file"]).read_text(encoding="utf-8")
        assert "expected_files=" in summary_text
        assert "present_files=" in summary_text


def test_gitignore_mentions_local_abaqus_and_validation_hygiene_patterns() -> None:
    p = PROJECT_ROOT / ".gitignore"
    text = p.read_text(encoding="utf-8")
    assert "abaqus.rpy*" in text
    assert "abaqus*.rec" in text
    assert "*_abaqus.odb" in text
    assert "post_processing/validation_visualisers/abaqus/generated/run_*.py" in text
    assert "post_processing/validation_visualisers/abaqus_results/*/" in text
    assert "post_processing/results/validation/" in text


def test_section_forces_comparison_runs():
    """Run section forces comparison; must not raise."""
    from post_processing.validation_visualisers.section_forces.section_forces_comparison import (
        run_section_forces_comparison,
    )
    run_section_forces_comparison()


def test_discover_job_names_from_results_returns_expected_names(tmp_path):
    """discover_job_names_from_results returns sorted unique job_XXXX_nN from mock result dirs."""
    from post_processing.validation_visualisers.job_discovery import (
        discover_job_names_from_results,
    )

    # Create mock result dirs matching timestamped pattern
    (tmp_path / "job_0001_n8_2026-02-26_16-19-42-602110_pid18316_b582c064").mkdir()
    (tmp_path / "job_0000_n16_2026-02-26_16-19-30-571529_pid18316_690d4b79").mkdir()
    (tmp_path / "job_0001_n8_2026-02-27_10-00-00-000000_pid999_abcdef01").mkdir()  # duplicate job
    (tmp_path / "job_0005_n32_2026-02-26_16-17-48-353941_pid18316_8d7c1563").mkdir()
    (tmp_path / "not_a_job_dir").mkdir()
    (tmp_path / "job_0000_n16_other_suffix").mkdir()  # no pid/hash, should not match

    names = discover_job_names_from_results(tmp_path)
    assert names == ["job_0000_n16", "job_0001_n8", "job_0005_n32"]

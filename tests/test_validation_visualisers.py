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


def test_simulation_type_dispatch_payload_reads_canonical_type() -> None:
    from post_processing.validation_visualisers.abaqus.simulation_type_dispatch import build_validation_dispatch_payload

    job_dir = PROJECT_ROOT / "jobs" / "job_smoke_eigen"
    payload = build_validation_dispatch_payload(str(job_dir))
    assert payload["simulation_type"] == "eigen"
    assert isinstance(payload["simulation_settings"], dict)


def test_extract_odb_results_mentions_tip_history_flag() -> None:
    p = PROJECT_ROOT / "post_processing" / "validation_visualisers" / "abaqus" / "extract_odb_results.py"
    text = p.read_text(encoding="utf-8")
    assert "export_tip_history" in text
    assert "tip_load_history.csv" in text


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

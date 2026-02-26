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


def test_deformation_comparison_runs():
    """Run deformation comparison (FEM vs Abaqus); must not raise when abaqus_results empty."""
    from post_processing.validation_visualisers.deflection_tables.deformation_comparison import (
        run_deformation_comparison,
    )
    run_deformation_comparison()


def test_section_forces_comparison_runs():
    """Run section forces comparison; must not raise."""
    from post_processing.validation_visualisers.section_forces.section_forces_comparison import (
        run_section_forces_comparison,
    )
    run_section_forces_comparison()

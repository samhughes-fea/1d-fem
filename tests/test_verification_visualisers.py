"""
Tests for post_processing verification_visualisers scripts.

Runs verification entry points and asserts they complete without exception.
When FEM results exist (post_processing/results/job_*), optionally asserts
that output CSVs exist and error columns are within tolerance.
"""

import sys
from pathlib import Path

import pytest

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))


def test_roark_verification_runs():
    """Run Roark displacement/rotation verification; must not raise."""
    from post_processing.verification_visualisers.roark.roark_verification import (
        run_roark_verification,
    )
    run_roark_verification()


def test_roark_section_forces_verification_runs():
    """Run Roark section forces (V, M) verification; must not raise."""
    from post_processing.verification_visualisers.roark.roark_section_forces_verification import (
        run_section_forces_verification,
    )
    run_section_forces_verification()


def test_deformation_convergence_runs():
    """Run deformation convergence; must not raise."""
    from post_processing.verification_visualisers.roark.deformation_convergence import (
        VisualiseDeformationConvergence,
    )
    VisualiseDeformationConvergence().process_convergence_plot()


def test_distributed_load_convergence_runs():
    """Run distributed-load convergence script; must not raise."""
    from post_processing.verification_visualisers.roark.distributed_load_convergence import (
        run_distributed_convergence,
    )
    run_distributed_convergence()


@pytest.mark.skipif(
    not (PROJECT_ROOT / "post_processing" / "verification_visualisers" / "roark" / "shear_deformable_verification.py").is_file(),
    reason="shear_deformable_verification.py not present",
)
def test_shear_deformable_verification_runs():
    """Run Timoshenko/Levinson verification; must not raise."""
    from post_processing.verification_visualisers.roark.shear_deformable_verification import (
        run_shear_deformable_verification,
    )
    run_shear_deformable_verification()


@pytest.mark.skipif(
    not (PROJECT_ROOT / "post_processing" / "results").is_dir(),
    reason="No post_processing/results dir",
)
def test_roark_verification_output_when_results_exist():
    """If any job results exist, Roark verification should produce deformation_plots dir (and optionally CSV)."""
    from post_processing.verification_visualisers.roark.roark_verification import (
        run_roark_verification,
    )
    run_roark_verification()
    out_dir = PROJECT_ROOT / "post_processing" / "verification_visualisers" / "roark" / "deformation_plots"
    assert out_dir.is_dir()

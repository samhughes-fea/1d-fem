"""
End-to-end: process_job with NonlinearTimoshenkoBeamElement3D and [warping]=1.

Uses a tiny 1-element, 2-node job in a temp directory to keep runtime small.
"""

from __future__ import annotations

import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from workflow_orchestrator.run_job import process_job, setup_job_results_directory


def _write_minimal_nl_ts_warp_job(job_dir: Path) -> None:
    job_dir.mkdir(parents=True, exist_ok=True)
    (job_dir / "grid.txt").write_text(
        "[Grid]\n"
        "[node_id]   [x]         [y]       [z]\n"
        "0           0.000000    0.0       0.0\n"
        "1           2.000000    0.0       0.0\n",
        encoding="utf-8",
    )
    (job_dir / "element.txt").write_text(
        "[Element]\n"
        "[element_id]  [node1]  [node2]  [element_type]                [axial_order]     "
        "[bending_y_order] [bending_z_order] [shear_y_order]   [shear_z_order]   [torsion_order]   [load_order]  [warping]\n"
        "0             0        1        NonlinearTimoshenkoBeamElement3D   3                 3                 3                 2                 2                 3                 2                 1\n",
        encoding="utf-8",
    )
    (job_dir / "material.txt").write_text(
        "[Material]\n[element_id]  [E]         [G]          [nu]   [rho]\n"
        "0             2.1e+11     8.1000e+10   0.3    7850\n",
        encoding="utf-8",
    )
    (job_dir / "section.txt").write_text(
        "[Section]\n"
        "[element_id]  [A]  [I_x]  [I_y]  [I_z]  [J_t]  [kappa]  [alpha]  [y_sc]  [z_sc]  [Gamma]\n"
        "0  0.00131  0.0  3.23400e-07  2.08769e-06  2.60673e-08  0.83333333333333337  0.0015936564885496184  0.0  0.0  1.0e-8\n",
        encoding="utf-8",
    )
    (job_dir / "precurvature.txt").write_text(
        "# t\n[Precurvature]\n[element_id] [k_x0] [k_y0] [k_z0]\n0 0.0 0.0 0.0\n",
        encoding="utf-8",
    )
    (job_dir / "prescribed_displacement.txt").write_text(
        "[Prescribed Displacement]\n"
        "[id]     [node_id]  [dof]   [value]     [type]          [comment]\n"
        "0        0          UX      0.0         displacement\n"
        "1        0          UY      0.0         displacement\n"
        "2        0          UZ      0.0         displacement\n"
        "3        0          RX      0.0         displacement\n"
        "4        0          RY      0.0         displacement\n"
        "5        0          RZ      0.0         displacement\n"
        "6        0          CHI     0.0         displacement\n"
        "7        1          CHI     0.0         displacement\n",
        encoding="utf-8",
    )
    (job_dir / "point_load.txt").write_text(
        "[Point load]\n"
        "         [x]          [y]          [z]        [F_x]        [F_y]        [F_z]        [M_x]        [M_y]        [M_z]\n"
        "    2.000000     0.000000     0.000000     0.000000  -100.000000     0.000000     0.000000     0.000000     0.000000\n",
        encoding="utf-8",
    )
    (job_dir / "simulation_settings.txt").write_text(
        "[Simulation]\n[Type]\nstatic_nonlinear\n\n[Newton]\n"
        "tolerance = 1e-4\n"
        "relative_tolerance = 1e-6\n"
        "relative_reference = first_residual\n"
        "max_iterations = 40\n"
        "tolerance_delta_u = 1e-9\n\n[Solver]\n"
        "tolerance = 1e-10\n"
        "max_iterations = 2000\n\n[Nonlinear]\n"
        "num_increments = 1\n"
        "line_search = false\n",
        encoding="utf-8",
    )


def test_process_job_minimal_nl_ts_warp_converges(tmp_path: Path) -> None:
    job_dir = tmp_path / "job_minimal_nl_ts_warp"
    _write_minimal_nl_ts_warp_job(job_dir)
    res_dir = setup_job_results_directory("test_minimal_nl_ts_warp")
    jt: dict = {}
    je: dict = {}
    process_job(
        str(job_dir),
        res_dir,
        jt,
        je,
        force_serial=True,
        max_processes_per_job=1,
    )
    manifest = Path(res_dir) / "logs" / "run_manifest.json"
    assert manifest.is_file()
    assert (Path(res_dir) / "primary_results").is_dir()

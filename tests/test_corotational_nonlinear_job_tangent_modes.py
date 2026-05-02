"""
Nonlinear static smoke: same mesh/load with CorotationalBeamElement3D and FD vs elastic_material tangents.

See ``docs/element_library/large_rotation_vs_total_lagrangian.md`` and plan epic 4.
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

import pytest

PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from workflow_orchestrator.run_job import process_job, setup_job_results_directory


def _write_corotational_nl_job(job_dir: Path, *, tangent_mode: str) -> None:
    job_dir.mkdir(parents=True, exist_ok=True)
    (job_dir / "grid.txt").write_text(
        "[Grid]\n[node_id]   [x]         [y]       [z]\n"
        "0           0.000000    0.0       0.0\n"
        "1           1.000000    0.0       0.0\n",
        encoding="utf-8",
    )
    (job_dir / "element.txt").write_text(
        "[Element]\n"
        "[element_id]  [node1]  [node2]  [element_type]                [axial_order]     "
        "[bending_y_order] [bending_z_order] [shear_y_order]   [shear_z_order]   [torsion_order]   [load_order]\n"
        "0             0        1        CorotationalBeamElement3D   3                 3                 "
        "3                 2                 2                 3                 2\n",
        encoding="utf-8",
    )
    (job_dir / "material.txt").write_text(
        "[Material]\n[element_id]  [E]         [G]          [nu]   [rho]\n"
        "0             2.1e+11     8.1e+10      0.3    7850\n",
        encoding="utf-8",
    )
    (job_dir / "section.txt").write_text(
        "[Section]\n"
        "[element_id]  [A]      [I_x]    [I_y]      [I_z]      [J_t]      [kappa]  [alpha]\n"
        "0             0.01     1e-9     1e-6       1e-6       2e-9       0.833333  0.0\n",
        encoding="utf-8",
    )
    (job_dir / "precurvature.txt").write_text(
        "# t\n[Precurvature]\n[element_id] [k_x0] [k_y0] [k_z0]\n0 0.0 0.0 0.0\n",
        encoding="utf-8",
    )
    (job_dir / "prescribed_displacement.txt").write_text(
        "[Prescribed Displacement]\n"
        "[id]     [node_id]  [dof]   [value]     [type]\n"
        "0        0          UX      0.0         displacement\n"
        "1        0          UY      0.0         displacement\n"
        "2        0          UZ      0.0         displacement\n"
        "3        0          RX      0.0         displacement\n"
        "4        0          RY      0.0         displacement\n"
        "5        0          RZ      0.0         displacement\n",
        encoding="utf-8",
    )
    (job_dir / "point_load.txt").write_text(
        "[Point load]\n"
        "         [x]          [y]          [z]        [F_x]        [F_y]        [F_z]        [M_x]        [M_y]        [M_z]\n"
        "    1.000000     0.000000     0.000000     0.000000  -500.000000     0.000000     0.000000     0.000000     0.000000\n",
        encoding="utf-8",
    )
    (job_dir / "simulation_settings.txt").write_text(
        "[Simulation]\n[Type]\nstatic_nonlinear\n\n[Newton]\n"
        "tolerance = 1e-5\n"
        "relative_tolerance = 1e-5\n"
        "relative_reference = first_residual\n"
        "max_iterations = 35\n"
        "tolerance_delta_u = 1e-9\n\n[Solver]\n"
        "tolerance = 1e-12\n"
        "max_iterations = 2000\n\n[Nonlinear]\n"
        "num_increments = 1\n"
        "line_search = false\n"
        f"corotational_tangent_mode = {tangent_mode}\n",
        encoding="utf-8",
    )


def _read_newton_converged(summary_csv: Path) -> bool:
    with summary_csv.open(newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    if not rows:
        return False
    last = rows[-1]
    raw = last.get("newton_converged", "")
    return str(raw).strip() in ("1", "True", "true")


@pytest.mark.parametrize("tangent_mode", ["finite_difference", "elastic_material"])
def test_corotational_nl_job_converges_both_tangent_modes(tmp_path: Path, tangent_mode: str) -> None:
    job_dir = tmp_path / f"job_corot_{tangent_mode}"
    _write_corotational_nl_job(job_dir, tangent_mode=tangent_mode)
    tag = f"test_corot_nl_{tangent_mode}"
    res_dir = setup_job_results_directory(tag)
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
    summary = Path(res_dir) / "primary_results" / "primary_summary.csv"
    assert summary.is_file(), f"missing {summary}"
    assert _read_newton_converged(summary), f"Newton did not converge for corotational_tangent_mode={tangent_mode}"
    nh = Path(res_dir) / "logs" / "newton_history.csv"
    assert nh.is_file()

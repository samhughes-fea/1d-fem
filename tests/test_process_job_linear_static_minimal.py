"""
End-to-end smoke: ``process_job`` with linear static ([Type] static) on a tiny beam.

Exercises ``LinearStaticSimulationRunner`` via ``run_job``. Section table uses 7 entries
(A … J_t, kappa, alpha) so ``section_array`` satisfies linear EB validation without warping.
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

import numpy as np
import pytest

PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from workflow_orchestrator.run_job import process_job, setup_job_results_directory


def _write_minimal_linear_static_job(job_dir: Path) -> None:
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
        "[bending_y_order] [bending_z_order] [shear_y_order]   [shear_z_order]   [torsion_order]   [load_order]\n"
        "0             0        1        LinearEulerBernoulliBeamElement3D   3                 3                 "
        "3                 0                 0                 3                 2\n",
        encoding="utf-8",
    )
    (job_dir / "material.txt").write_text(
        "[Material]\n[element_id]  [E]         [G]          [nu]   [rho]\n"
        "0             2.1e+11     8.1000e+10   0.3    7850\n",
        encoding="utf-8",
    )
    # Seven section scalars (no y_sc/z_sc/Gamma) → ``section_array.size == 7`` for non-warp mesh.
    (job_dir / "section.txt").write_text(
        "[Section]\n"
        "[element_id]  [A]  [I_x]  [I_y]  [I_z]  [J_t]  [kappa]  [alpha]\n"
        "0  0.00131  0.0  3.23400e-07  2.08769e-06  2.60673e-08  0.83333333333333337  0.0015936564885496184\n",
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
        "5        0          RZ      0.0         displacement\n",
        encoding="utf-8",
    )
    (job_dir / "point_load.txt").write_text(
        "[Point load]\n"
        "         [x]          [y]          [z]        [F_x]        [F_y]        [F_z]        [M_x]        [M_y]        [M_z]\n"
        "    2.000000     0.000000     0.000000     0.000000  -1000.000000     0.000000     0.000000     0.000000     0.000000\n",
        encoding="utf-8",
    )
    (job_dir / "simulation_settings.txt").write_text(
        "[Simulation]\n[Type]\nstatic\n\n[Solver]\n"
        "type = direct\n"
        "tolerance = 1e-12\n"
        "max_iterations = 500\n",
        encoding="utf-8",
    )


def test_process_job_minimal_linear_static(tmp_path: Path) -> None:
    job_dir = tmp_path / "job_minimal_linear_static"
    _write_minimal_linear_static_job(job_dir)
    res_dir = setup_job_results_directory("test_minimal_linear_static")
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
    root = Path(res_dir)
    assert (root / "logs" / "run_manifest.json").is_file()
    assert (root / "primary_results").is_dir()
    u_path = root / "primary_results" / "global" / "U_global.csv"
    assert u_path.is_file()
    with open(u_path, newline="", encoding="utf-8") as f:
        rows = list(csv.reader(f))
    assert len(rows) >= 2
    vals = [float(x) for x in rows[1][1:]]
    assert np.all(np.isfinite(vals))

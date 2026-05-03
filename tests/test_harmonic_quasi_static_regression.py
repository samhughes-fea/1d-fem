"""Quasi-static check: harmonic response at very low frequency matches linear static ``U_global``."""

from __future__ import annotations

import csv
import shutil
import sys
from pathlib import Path

import numpy as np
import pytest

PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from workflow_orchestrator.run_job import process_job, setup_job_results_directory


def _read_u_global_csv(path: Path) -> np.ndarray:
    rows: list[float] = []
    with path.open(newline="", encoding="utf-8") as f:
        reader = csv.reader(f)
        next(reader, None)
        for row in reader:
            rows.append(float(row[1]))
    return np.asarray(rows, dtype=np.float64)


@pytest.mark.integration
def test_harmonic_very_low_frequency_matches_linear_static(tmp_path: Path) -> None:
    """Same mesh/load as ``job_smoke_harmonic``; compare static ``U`` to harmonic column at min Hz."""
    src = PROJECT_ROOT / "jobs" / "job_smoke_harmonic"
    j_static = tmp_path / "job_hs_static"
    j_harm = tmp_path / "job_hs_harmonic"
    shutil.copytree(src, j_static)
    shutil.copytree(src, j_harm)

    (j_static / "simulation_settings.txt").write_text(
        "[Simulation]\n[Type]\nstatic\n",
        encoding="utf-8",
    )
    (j_harm / "simulation_settings.txt").write_text(
        "[Simulation]\n[Type]\nharmonic\n\n"
        "[Harmonic]\n"
        "enabled = true\n"
        "frequency_min_hz = 1e-8\n"
        "frequency_max_hz = 2e-8\n"
        "num_frequency_points = 2\n"
        "modal_damping_ratio = 0.0\n",
        encoding="utf-8",
    )

    res_s = setup_job_results_directory("pytest_harmonic_quasi_static_ref")
    res_h = setup_job_results_directory("pytest_harmonic_quasi_static_sw")
    process_job(str(j_static), res_s, {}, {}, force_serial=True)
    process_job(str(j_harm), res_h, {}, {}, force_serial=True)

    u_static = _read_u_global_csv(Path(res_s) / "primary_results" / "global" / "U_global.csv")

    harm_dir = Path(res_h) / "primary_results" / "harmonic_results"
    real_files = sorted(harm_dir.glob("*_displacement_real.txt"))
    assert len(real_files) == 1
    u_mat = np.loadtxt(real_files[0])
    u_h = u_mat[:, 0] if u_mat.ndim == 2 else u_mat

    np.testing.assert_allclose(u_static, u_h, rtol=5e-4, atol=1e-7)

    imag_files = sorted(harm_dir.glob("*_displacement_imag.txt"))
    assert len(imag_files) == 1
    ui = np.loadtxt(imag_files[0])
    ui_col = ui[:, 0] if ui.ndim == 2 else ui
    assert np.linalg.norm(ui_col, ord=np.inf) < 1e-6 * max(1.0, np.linalg.norm(u_static, ord=np.inf))

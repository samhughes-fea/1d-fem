from __future__ import annotations

import csv
import math
import sys
from pathlib import Path

_HERE = Path(__file__).resolve()
_PROJECT_ROOT = _HERE.parents[2]
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))

from validation.abaqus.config import ABAQUS_RESULTS_DIR, PROJECT_ROOT
from validation.common_runner_contract import build_runner_result, write_runner_summary
from validation.transient.analytical_transient import (
    forced_sdof_steady_state_amplitude,
    load_transient_sdof_params,
)

_TRANSIENT_JOB_DIR = PROJECT_ROOT / "jobs" / "job_benchmark_transient_cantilever_multidof"
_ANALYTICAL_TOLERANCE = 0.50  # 50 % relative; SDOF vs multi-DOF FEM is a rough approximation
_MAX_ABS_ERROR_TOL = 1.0e-8   # FEM vs Abaqus comparison tolerance


def _read_transient_history(csv_path: Path) -> list[tuple[float, float]]:
    rows: list[tuple[float, float]] = []
    with open(csv_path, "r", encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            rows.append((float(row["time"]), float(row["tip_displacement"])))
    return rows


def _max_abs_tip(rows: list[tuple[float, float]]) -> float:
    return max((abs(r[1]) for r in rows), default=0.0)


def _interpolate_at_time(history: list[tuple[float, float]], t: float) -> float:
    if not history:
        return 0.0
    history = sorted(history, key=lambda r: r[0])
    if t <= history[0][0]:
        return history[0][1]
    if t >= history[-1][0]:
        return history[-1][1]
    for (t0, u0), (t1, u1) in zip(history, history[1:]):
        if t0 <= t <= t1:
            frac = (t - t0) / (t1 - t0) if not math.isclose(t1, t0, abs_tol=1e-15) else 0.0
            return u0 + frac * (u1 - u0)
    return history[-1][1]


def _compare_transient_histories(
    fem_rows: list[tuple[float, float]],
    abaqus_rows: list[tuple[float, float]],
) -> dict[str, float | int | bool | str]:
    max_abs_error = 0.0
    for t_fem, u_fem in fem_rows:
        u_aba = _interpolate_at_time(abaqus_rows, t_fem)
        max_abs_error = max(max_abs_error, abs(u_fem - u_aba))
    passed = max_abs_error <= _MAX_ABS_ERROR_TOL
    return {
        "ready": True,
        "passed": passed,
        "status": "pass" if passed else "fail",
        "n_rows_fem": len(fem_rows),
        "n_rows_abaqus": len(abaqus_rows),
        "max_abs_error": max_abs_error,
    }


def _analytical_comparison(fem_csv: Path) -> dict[str, float | bool]:
    if not fem_csv.is_file():
        return {}
    try:
        rows = _read_transient_history(fem_csv)
        if not rows:
            return {}
        fem_max = _max_abs_tip(rows)
        params = load_transient_sdof_params(_TRANSIENT_JOB_DIR)
        analytical_amp = forced_sdof_steady_state_amplitude(params)
        if analytical_amp <= 0 or not math.isfinite(analytical_amp):
            return {}
        rel_error = abs(fem_max - analytical_amp) / analytical_amp
        return {
            "analytical_ready": True,
            "analytical_passed": bool(rel_error <= _ANALYTICAL_TOLERANCE),
            "analytical_max_relative_error": rel_error,
        }
    except Exception:
        return {}


def run_transient_validation(results_root: str | Path) -> dict[str, str | bool]:
    """Canonical comparator runner for the transient cantilever benchmark.

    Status vocabulary: pass, fail, missing_artifacts, analytical_only_pass.
    analytical_only_pass is issued when FEM results are present but Abaqus
    reference artifacts are absent; the FEM max tip displacement is compared
    against the first-mode SDOF analytical amplitude within 50 % tolerance.
    """
    results_root = Path(results_root)
    out_dir = results_root / "validation" / "transient"
    out_dir.mkdir(parents=True, exist_ok=True)
    summary = out_dir / "transient_validation_summary.txt"
    abaqus_job_name = "job_benchmark_transient_cantilever_multidof_n500"
    reference_dir = ABAQUS_RESULTS_DIR / abaqus_job_name
    expected = ["transient_reference.csv"]
    present = [name for name in expected if (reference_dir / name).is_file()]
    fem_csv = results_root / "primary_results" / "transient" / "tip_displacement_history.csv"
    abaqus_csv = reference_dir / "transient_reference.csv"
    tolerances = {
        "max_abs_error": _MAX_ABS_ERROR_TOL,
        "analytical_amplitude_tolerance": _ANALYTICAL_TOLERANCE,
    }

    analytical = _analytical_comparison(fem_csv)
    analytical_ready = bool(analytical.get("analytical_ready", False))
    analytical_passed = bool(analytical.get("analytical_passed", False))

    fem_present = fem_csv.is_file()
    abaqus_present = abaqus_csv.is_file()

    if fem_present and abaqus_present:
        result = _compare_transient_histories(
            _read_transient_history(fem_csv),
            _read_transient_history(abaqus_csv),
        )
        ready = True
        passed = bool(result["passed"])
        status = str(result["status"])
        metrics: dict[str, float | int | str | bool] = {
            "n_rows_fem": int(result["n_rows_fem"]),
            "n_rows_abaqus": int(result["n_rows_abaqus"]),
            "max_abs_error": float(result["max_abs_error"]),
            "analytical_passed": analytical_passed,
        }
    elif fem_present:
        ready = analytical_ready
        passed = analytical_passed
        status = "analytical_only_pass" if analytical_passed else "fail"
        metrics = {
            "n_rows_fem": 0,
            "n_rows_abaqus": 0,
            "max_abs_error": float("nan"),
            "analytical_passed": analytical_passed,
            "analytical_max_relative_error": float(
                analytical.get("analytical_max_relative_error", float("nan"))
            ),
        }
    else:
        ready = False
        passed = False
        status = "missing_artifacts"
        metrics = {
            "n_rows_fem": 0,
            "n_rows_abaqus": 0,
            "max_abs_error": float("nan"),
            "analytical_passed": False,
        }

    artifact_paths = {"fem_csv": fem_csv, "abaqus_csv": abaqus_csv}
    write_runner_summary(
        summary,
        family="transient",
        ready=ready,
        passed=passed,
        status=status,
        reference_dir=reference_dir,
        expected_files=expected,
        present_files=present,
        tolerances=tolerances,
        metrics=metrics,
        artifact_paths=artifact_paths,
    )
    return build_runner_result(
        family="transient",
        summary_file=summary,
        reference_dir=reference_dir,
        expected_files=expected,
        present_files=present,
        ready=ready,
        passed=passed,
        status=status,
        tolerances=tolerances,
        metrics=metrics,
        artifact_paths=artifact_paths,
    )

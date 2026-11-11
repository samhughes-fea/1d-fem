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
from validation.nonlinear_buckling.analytical_nl_buckling import (
    load_euler_column_params,
)

_NL_BUCKLING_JOB_DIR = PROJECT_ROOT / "jobs" / "job_benchmark_nl_buckling_imperfect_column"
# Path-comparison tolerance: generous for post-buckling trajectory agreement.
_MAX_PATH_ERROR_TOL = 1.0e-4
# analytical_only_pass: FEM load factors must stay below Euler λ_cr, and
# at least one increment must have converged.
_EULER_MARGIN = 0.99  # FEM max load factor must be < EULER_MARGIN * λ_cr


def _read_continuation_history(csv_path: Path) -> list[dict[str, float | bool]]:
    rows: list[dict[str, float | bool]] = []
    with open(csv_path, "r", encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            rows.append(
                {
                    "load_factor": float(row["load_factor"]),
                    "converged": str(row.get("converged", "True")).strip().lower() in ("true", "1"),
                    "tip_displacement": float(row["tip_displacement"]),
                }
            )
    return rows


def _interpolate_tip_at_load(
    history: list[dict[str, float | bool]], target_lf: float
) -> float:
    sorted_h = sorted(history, key=lambda r: r["load_factor"])
    if target_lf <= sorted_h[0]["load_factor"]:
        return float(sorted_h[0]["tip_displacement"])
    if target_lf >= sorted_h[-1]["load_factor"]:
        return float(sorted_h[-1]["tip_displacement"])
    for r0, r1 in zip(sorted_h, sorted_h[1:]):
        lf0, lf1 = float(r0["load_factor"]), float(r1["load_factor"])
        if lf0 <= target_lf <= lf1:
            frac = (target_lf - lf0) / (lf1 - lf0) if not math.isclose(lf1, lf0, abs_tol=1e-15) else 0.0
            return float(r0["tip_displacement"]) + frac * (float(r1["tip_displacement"]) - float(r0["tip_displacement"]))
    return float(sorted_h[-1]["tip_displacement"])


def _compare_paths(
    fem_rows: list[dict[str, float | bool]],
    ref_rows: list[dict[str, float | bool]],
) -> dict[str, float | int | bool | str]:
    max_path_error = 0.0
    for fem_r in fem_rows:
        ref_tip = _interpolate_tip_at_load(ref_rows, float(fem_r["load_factor"]))
        max_path_error = max(max_path_error, abs(float(fem_r["tip_displacement"]) - ref_tip))
    passed = max_path_error <= _MAX_PATH_ERROR_TOL
    return {
        "ready": True,
        "passed": passed,
        "status": "pass" if passed else "fail",
        "n_rows_fem": len(fem_rows),
        "n_rows_reference": len(ref_rows),
        "max_path_error": max_path_error,
    }


def _analytical_comparison(fem_csv: Path) -> dict[str, float | bool]:
    if not fem_csv.is_file():
        return {}
    try:
        rows = _read_continuation_history(fem_csv)
        if not rows:
            return {}
        params = load_euler_column_params(_NL_BUCKLING_JOB_DIR)
        lambda_cr = params["lambda_cr"]
        max_lf = max(float(r["load_factor"]) for r in rows)
        n_converged = sum(1 for r in rows if r["converged"])
        # Physical sanity: pre-buckling regime and at least one converged step.
        pre_buckling = max_lf < _EULER_MARGIN * lambda_cr
        has_convergence = n_converged >= 1
        passed = pre_buckling and has_convergence
        return {
            "analytical_ready": True,
            "analytical_passed": passed,
            "analytical_lambda_cr": lambda_cr,
            "analytical_max_load_factor": max_lf,
            "analytical_n_converged": n_converged,
        }
    except Exception:
        return {}


def run_nonlinear_buckling_validation(results_root: str | Path) -> dict[str, str | bool]:
    """Canonical comparator runner for the nonlinear buckling imperfect-column benchmark.

    Status vocabulary: pass, fail, missing_artifacts, analytical_only_pass.

    analytical_only_pass requires:
    - at least one converged increment
    - all FEM load factors below 99 % of the Euler critical load factor

    pass requires FEM and Abaqus path agreement within max_path_error ≤ 1e-4.
    """
    results_root = Path(results_root)
    out_dir = results_root / "validation" / "nonlinear_buckling"
    out_dir.mkdir(parents=True, exist_ok=True)
    summary = out_dir / "nonlinear_buckling_validation_summary.txt"
    abaqus_job_name = "job_benchmark_nl_buckling_imperfect_column_n500"
    reference_dir = ABAQUS_RESULTS_DIR / abaqus_job_name
    expected = ["continuation_history.csv"]
    present = [name for name in expected if (reference_dir / name).is_file()]
    fem_csv = (
        results_root
        / "primary_results"
        / "nonlinear_buckling_results"
        / "continuation_history.csv"
    )
    abaqus_csv = reference_dir / "continuation_history.csv"
    tolerances = {
        "max_path_error": _MAX_PATH_ERROR_TOL,
        "euler_pre_buckling_margin": _EULER_MARGIN,
    }

    analytical = _analytical_comparison(fem_csv)
    analytical_ready = bool(analytical.get("analytical_ready", False))
    analytical_passed = bool(analytical.get("analytical_passed", False))

    fem_present = fem_csv.is_file()
    abaqus_present = abaqus_csv.is_file()

    if fem_present and abaqus_present:
        result = _compare_paths(
            _read_continuation_history(fem_csv),
            _read_continuation_history(abaqus_csv),
        )
        ready = True
        passed = bool(result["passed"])
        status = str(result["status"])
        metrics: dict[str, float | int | str | bool] = {
            "n_rows_fem": int(result["n_rows_fem"]),
            "n_rows_reference": int(result["n_rows_reference"]),
            "max_path_error": float(result["max_path_error"]),
            "analytical_passed": analytical_passed,
        }
    elif fem_present:
        ready = analytical_ready
        passed = analytical_passed
        status = "analytical_only_pass" if analytical_passed else "fail"
        metrics = {
            "n_rows_fem": 0,
            "n_rows_reference": 0,
            "max_path_error": float("nan"),
            "analytical_passed": analytical_passed,
            "analytical_lambda_cr": float(analytical.get("analytical_lambda_cr", float("nan"))),
            "analytical_n_converged": int(analytical.get("analytical_n_converged", 0)),
        }
    else:
        ready = False
        passed = False
        status = "missing_artifacts"
        metrics = {
            "n_rows_fem": 0,
            "n_rows_reference": 0,
            "max_path_error": float("nan"),
            "analytical_passed": False,
        }

    artifact_paths = {"fem_csv": fem_csv, "abaqus_csv": abaqus_csv}
    write_runner_summary(
        summary,
        family="nonlinear_buckling",
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
        family="nonlinear_buckling",
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

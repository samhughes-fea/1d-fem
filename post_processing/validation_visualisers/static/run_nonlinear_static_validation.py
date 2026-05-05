from __future__ import annotations

import csv
from pathlib import Path

from post_processing.validation_visualisers.abaqus.config import ABAQUS_RESULTS_DIR
from post_processing.validation_visualisers.static.compare_nonlinear_static_tip_history import (
    run_nonlinear_static_tip_history_validation_from_benchmark_paths,
)


MAX_ABS_ERROR_TOL = 1.0e-5
MAX_LOAD_ALIGNMENT_ERROR_TOL = 1.0e-12


def evaluate_benchmark_pass_fail(metrics: dict[str, float]) -> bool:
    return (
        float(metrics.get("max_abs_error", float("inf"))) <= MAX_ABS_ERROR_TOL
        and float(metrics.get("max_load_factor_alignment_error", float("inf"))) <= MAX_LOAD_ALIGNMENT_ERROR_TOL
    )


def run_benchmark(
    results_root: str | Path,
    *,
    job_name: str = "job_benchmark_nl_static_cantilever_tip",
    abaqus_reference_job_name: str | None = None,
) -> dict[str, str | bool]:
    return run_nonlinear_static_tip_history_validation_from_benchmark_paths(
        results_root=results_root,
        abaqus_root=ABAQUS_RESULTS_DIR,
        job_name=job_name,
        abaqus_reference_job_name=abaqus_reference_job_name,
    )


def discover_suite_job_names() -> list[str]:
    jobs_root = Path(__file__).resolve().parents[3] / "jobs"
    names = []
    for p in sorted(jobs_root.iterdir()):
        if not p.is_dir():
            continue
        nm = p.name
        if nm.startswith("job_benchmark_nl_static_"):
            names.append(nm)
    return names


def write_suite_summary(rows: list[dict[str, str | bool]], out_path: str | Path) -> None:
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow([
            "job_name",
            "ready",
            "passed",
            "status_file",
            "max_load_factor_alignment_error",
            "max_abs_error",
            "mean_abs_error",
        ])
        for r in rows:
            w.writerow([
                r.get("job_name"),
                r.get("ready"),
                r.get("passed", ""),
                r.get("status_file"),
                r.get("max_load_factor_alignment_error", ""),
                r.get("max_abs_error", ""),
                r.get("mean_abs_error", ""),
            ])


def _read_summary_metrics(summary_csv: str | Path) -> dict[str, float]:
    out: dict[str, float] = {}
    with open(summary_csv, "r", encoding="utf-8") as f:
        r = csv.reader(f)
        next(r, None)
        for row in r:
            if len(row) >= 2 and row[0] in (
                "max_load_factor_alignment_error",
                "max_abs_error",
                "mean_abs_error",
            ):
                out[row[0]] = float(row[1])
    return out


def run_suite(results_root: str | Path) -> dict[str, str | list[dict[str, str | bool]]]:
    rows: list[dict[str, str | bool]] = []
    for job_name in discover_suite_job_names():
        out = run_benchmark(results_root, job_name=job_name)
        out["job_name"] = job_name
        if out.get("ready") and out.get("summary_csv"):
            metrics = _read_summary_metrics(out["summary_csv"])
            out.update(metrics)
            out["passed"] = evaluate_benchmark_pass_fail(metrics)
        rows.append(out)
    suite_dir = Path(results_root) / "validation" / "nonlinear_static_suite"
    summary_csv = suite_dir / "suite_summary.csv"
    write_suite_summary(rows, summary_csv)
    return {"rows": rows, "summary_csv": str(summary_csv)}

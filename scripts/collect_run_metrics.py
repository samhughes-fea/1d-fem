"""
Aggregate per-job result folders under ``post_processing/results`` into one CSV.

Reads ``logs/run_manifest.json`` when present and merges ``primary_summary.csv``
and row counts for ``newton_history.csv`` / ``inner_solve_history.csv``.

Usage (from repo root)::

    python scripts/collect_run_metrics.py
    python scripts/collect_run_metrics.py --out post_processing/runs_summary.csv
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import sys
from pathlib import Path

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)


def _count_csv_rows(path: Path, subtract_header: bool = True) -> int | None:
    if not path.is_file():
        return None
    try:
        with open(path, newline="", encoding="utf-8") as f:
            n = sum(1 for _ in f)
        return max(0, n - (1 if subtract_header and n else 0))
    except OSError:
        return None


def _read_primary_summary_row(summary_path: Path) -> dict[str, str]:
    out: dict[str, str] = {}
    if not summary_path.is_file():
        return out
    try:
        with open(summary_path, newline="", encoding="utf-8") as f:
            r = csv.DictReader(f)
            row = next(iter(r), None)
            if row:
                out = {k: (v or "") for k, v in row.items()}
    except OSError:
        pass
    return out


def collect_runs(
    results_root: Path,
) -> list[dict[str, str | int | float | None]]:
    rows: list[dict[str, str | int | float | None]] = []
    if not results_root.is_dir():
        return rows

    for d in sorted(results_root.iterdir(), key=lambda p: p.name):
        if not d.is_dir():
            continue
        manifest_path = d / "logs" / "run_manifest.json"
        job_name = d.name
        wall = None
        git = None
        if manifest_path.is_file():
            try:
                with open(manifest_path, encoding="utf-8") as f:
                    m = json.load(f)
                job_name = str(m.get("job_name", job_name))
                wall = m.get("wall_time_sec")
                git = m.get("git_commit")
            except (OSError, json.JSONDecodeError):
                pass

        nh = d / "logs" / "newton_history.csv"
        ish = d / "logs" / "inner_solve_history.csv"
        ps = d / "primary_results" / "primary_summary.csv"
        prim = _read_primary_summary_row(ps)

        row: dict[str, str | int | float | None] = {
            "results_folder": d.name,
            "job_name": job_name,
            "git_commit": git,
            "wall_time_sec": wall,
            "n_newton_history_rows": _count_csv_rows(nh),
            "n_inner_solve_rows": _count_csv_rows(ish),
        }
        for k in (
            "newton_iterations_total",
            "newton_converged",
            "load_increments_completed",
            "total_dof",
            "n_elements",
        ):
            if k in prim:
                row[f"primary_{k}"] = prim[k]
        rows.append(row)
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description="Summarize job result folders.")
    parser.add_argument(
        "--results-root",
        type=Path,
        default=Path(ROOT) / "post_processing" / "results",
        help="Directory containing per-run result folders",
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=Path(ROOT) / "post_processing" / "runs_summary.csv",
        help="Output CSV path",
    )
    args = parser.parse_args()

    rows = collect_runs(args.results_root)
    if not rows:
        print(f"No subdirectories under {args.results_root}", file=sys.stderr)
        return

    fieldnames: list[str] = []
    for r in rows:
        for k in r:
            if k not in fieldnames:
                fieldnames.append(k)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in fieldnames})

    print(f"Wrote {len(rows)} rows -> {args.out.resolve()}")


if __name__ == "__main__":
    main()

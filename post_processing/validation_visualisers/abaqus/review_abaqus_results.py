# post_processing/validation_visualisers/abaqus/review_abaqus_results.py
"""
Review abaqus_results directories: file completeness and Abaqus run performance.

Scans post_processing/validation_visualisers/abaqus_results/, checks expected files
per job dir, parses .sta (completion status, total time) and .msg (error/warning counts),
writes a CSV report, optional Markdown summary, and an errors/inconsistencies log.
Run from project root:

  python post_processing/validation_visualisers/abaqus/review_abaqus_results.py
  python post_processing/validation_visualisers/abaqus/review_abaqus_results.py --expected --output output/abaqus_results_review.csv
  python post_processing/validation_visualisers/abaqus/review_abaqus_results.py --no-log-errors
"""
from __future__ import annotations

import argparse
import csv
import re
import sys
from datetime import datetime
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
VALIDATION_DIR = SCRIPT_DIR.parent
PROJECT_ROOT = VALIDATION_DIR.parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from post_processing.validation_visualisers.abaqus.config import (
    ABAQUS_RESULTS_DIR,
    JOBS_DIR,
)

# job_XXXX_nN (e.g. job_0000_n8, job_0007_n128)
JOB_DIR_PATTERN = re.compile(r"^job_\d{4}_n\d+$")

DEFAULT_OUTPUT_CSV = VALIDATION_DIR / "output" / "abaqus_results_review.csv"
DEFAULT_OUTPUT_MD = VALIDATION_DIR / "output" / "abaqus_performance_summary.md"
DEFAULT_OUTPUT_LOG = VALIDATION_DIR / "output" / "abaqus_results_errors_and_inconsistencies.log"
MAX_MSG_LINES_PER_JOB = 20


def discover_result_jobs() -> list[str]:
    """Return sorted list of job names (job_XXXX_nN) that exist under ABAQUS_RESULTS_DIR."""
    if not ABAQUS_RESULTS_DIR.is_dir():
        return []
    names = []
    for p in ABAQUS_RESULTS_DIR.iterdir():
        if p.is_dir() and JOB_DIR_PATTERN.match(p.name):
            names.append(p.name)
    return sorted(names)


def discover_expected_jobs() -> list[str]:
    """Return sorted list of job names that exist under JOBS_DIR (expected to have results)."""
    if not JOBS_DIR.is_dir():
        return []
    names = []
    for p in JOBS_DIR.iterdir():
        if p.is_dir() and JOB_DIR_PATTERN.match(p.name):
            names.append(p.name)
    return sorted(names)


def _parse_sta(sta_path: Path) -> tuple[str, float | None]:
    """
    Parse .sta file. Return (status, total_time_sec).
    status: COMPLETED | ABORTED | NO_STA | UNKNOWN.
    total_time_sec: from SUMMARY row TOTAL TIME/FREQ column if present, else None.
    """
    if not sta_path.is_file():
        return ("NO_STA", None)
    try:
        text = sta_path.read_text(encoding="utf-8", errors="replace")
    except OSError:
        return ("NO_STA", None)
    lines = text.splitlines()
    status = "UNKNOWN"
    if any("THE ANALYSIS HAS COMPLETED SUCCESSFULLY" in line for line in lines):
        status = "COMPLETED"
    elif any("THE ANALYSIS WAS NOT COMPLETED" in line or "ABORTED" in line for line in lines):
        status = "ABORTED"

    total_time_sec: float | None = None
    for line in lines:
        # Summary data row: "   1     1   1     0     1     1  1.00       1.00       1.000"
        stripped = line.strip()
        if not stripped or not stripped[0].isdigit():
            continue
        parts = stripped.split()
        floats = []
        for p in parts:
            try:
                floats.append(float(p))
            except ValueError:
                break
        if len(floats) >= 7:
            # Column 7 (0-indexed 6) is TOTAL TIME/FREQ in seconds (Abaqus .sta format)
            total_time_sec = floats[6]
            break

    return (status, total_time_sec)


def _parse_msg(msg_path: Path) -> tuple[int, int]:
    """Count lines containing 'ERROR' and 'WARNING' in .msg file. Return (error_count, warning_count)."""
    if not msg_path.is_file():
        return (0, 0)
    try:
        text = msg_path.read_text(encoding="utf-8", errors="replace")
    except OSError:
        return (0, 0)
    errors = sum(1 for line in text.splitlines() if "ERROR" in line.upper())
    warnings = sum(1 for line in text.splitlines() if "WARNING" in line.upper())
    return (errors, warnings)


def _collect_errors_and_inconsistencies(
    job_dir: Path,
    job_name: str,
    base: str,
    status: str,
    has_u: bool,
    has_rot: bool,
    has_log: bool,
) -> list[tuple[str, str]]:
    """
    Read run_log.txt, rotation_source.txt, .msg, and check U_global.csv for this job.
    Return list of (source, message) for errors and inconsistencies.
    """
    out: list[tuple[str, str]] = []
    run_log_path = job_dir / "run_log.txt"
    rot_path = job_dir / "rotation_source.txt"
    msg_path = job_dir / f"{base}.msg"
    u_path = job_dir / "U_global.csv"

    # Inconsistencies: missing required files
    if not has_u:
        out.append(("[inconsistency]", "missing U_global.csv"))
    if not has_rot:
        out.append(("[inconsistency]", "missing rotation_source.txt"))
    if not has_log:
        out.append(("[inconsistency]", "missing run_log.txt"))
    if status != "COMPLETED":
        out.append(("[inconsistency]", f"analysis status is {status!r}, not COMPLETED"))

    # run_log.txt: error lines and ODB has UR
    odb_has_ur: bool | None = None
    if run_log_path.is_file():
        try:
            lines = run_log_path.read_text(encoding="utf-8", errors="replace").splitlines()
        except OSError:
            lines = []
        for line in lines:
            line_stripped = line.strip()
            if "ODB has UR:" in line_stripped:
                odb_has_ur = "True" in line_stripped
                out.append(("[run_log]", line_stripped))
            elif "failed" in line_stripped.lower() or "error" in line_stripped.lower():
                out.append(("[run_log]", line_stripped))
    else:
        lines = []

    # rotation_source.txt
    rotation_source: str | None = None
    if rot_path.is_file():
        try:
            rot_lines = rot_path.read_text(encoding="utf-8", errors="replace").splitlines()
            rotation_source = rot_lines[0].strip() if rot_lines else ""
        except OSError:
            rotation_source = ""
        if rotation_source not in ("ODB", "derived", "none"):
            out.append(("[inconsistency]", f"rotation_source.txt content is {rotation_source!r}, expected ODB, derived, or none"))
    else:
        rotation_source = None

    # ODB has UR False but rotation_source ODB
    if odb_has_ur is False and rotation_source == "ODB":
        out.append(("[inconsistency]", "rotation_source is \"ODB\" but run_log contains \"ODB has UR: False\""))

    # U_global.csv: header and non-empty
    if u_path.is_file():
        try:
            u_text = u_path.read_text(encoding="utf-8", errors="replace")
            u_lines = u_text.splitlines()
        except OSError:
            out.append(("[inconsistency]", "U_global.csv could not be read"))
        else:
            if not u_lines:
                out.append(("[inconsistency]", "U_global.csv is empty"))
            elif u_lines[0].strip() != "Global DOF,Value":
                out.append(("[inconsistency]", f"U_global.csv header is {u_lines[0][:50]!r}, expected \"Global DOF,Value\""))
            elif len(u_lines) < 2:
                out.append(("[inconsistency]", "U_global.csv has no data rows"))
    # else already reported missing

    # .msg: lines containing ERROR or WARNING (capped)
    if msg_path.is_file():
        try:
            msg_lines = msg_path.read_text(encoding="utf-8", errors="replace").splitlines()
        except OSError:
            msg_lines = []
        collected = []
        for i, line in enumerate(msg_lines, start=1):
            if "ERROR" in line.upper() or "WARNING" in line.upper():
                collected.append((i, line.strip()[:200]))  # truncate long lines
        if collected:
            for idx, (lnum, content) in enumerate(collected):
                if idx >= MAX_MSG_LINES_PER_JOB:
                    remaining = len(collected) - MAX_MSG_LINES_PER_JOB
                    out.append(("[.msg]", f"... {remaining} more ERROR/WARNING lines"))
                    break
                out.append(("[.msg]", f"line {lnum}: {content}"))

    return out


def run_review(
    *,
    include_expected: bool = False,
    output_csv: Path,
    output_md: Path | None = None,
    output_log: Path | None = None,
) -> int:
    """
    Scan abaqus_results, collect per-job file presence and performance, write report(s).
    Optionally write errors/inconsistencies log. Return 0 on success.
    """
    result_jobs = discover_result_jobs()
    expected_jobs = discover_expected_jobs() if include_expected else []

    rows: list[dict[str, str | int | float | None]] = []
    errors_per_job: list[tuple[str, list[tuple[str, str]]]] = []

    for job_name in result_jobs:
        job_dir = ABAQUS_RESULTS_DIR / job_name
        base = f"{job_name}_abaqus"
        has_u = (job_dir / "U_global.csv").is_file()
        has_sf = (job_dir / "section_forces.csv").is_file()
        has_nodal_sf = (job_dir / "nodal_section_forces.csv").is_file()
        has_rot = (job_dir / "rotation_source.txt").is_file()
        has_log = (job_dir / "run_log.txt").is_file()
        has_inp = (job_dir / f"{base}.inp").is_file()
        has_odb = (job_dir / f"{base}.odb").is_file()
        has_sta = (job_dir / f"{base}.sta").is_file()
        has_msg = (job_dir / f"{base}.msg").is_file()

        status = "NO_STA"
        total_time_sec: float | None = None
        msg_errors = 0
        msg_warnings = 0

        if has_sta:
            status, total_time_sec = _parse_sta(job_dir / f"{base}.sta")
        if has_msg:
            msg_errors, msg_warnings = _parse_msg(job_dir / f"{base}.msg")

        rows.append({
            "job_name": job_name,
            "has_U_global": has_u,
            "has_section_forces": has_sf,
            "has_nodal_section_forces": has_nodal_sf,
            "has_rotation_source": has_rot,
            "has_run_log": has_log,
            "has_inp": has_inp,
            "has_odb": has_odb,
            "has_sta": has_sta,
            "has_msg": has_msg,
            "status": status,
            "total_time_sec": total_time_sec,
            "msg_errors": msg_errors,
            "msg_warnings": msg_warnings,
        })

        if output_log is not None:
            entries = _collect_errors_and_inconsistencies(
                job_dir, job_name, base, status, has_u, has_rot, has_log
            )
            errors_per_job.append((job_name, entries))

    # Write CSV
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "job_name", "has_U_global", "has_section_forces", "has_nodal_section_forces", "has_rotation_source",
        "has_run_log", "has_inp", "has_odb", "has_sta", "has_msg",
        "status", "total_time_sec", "msg_errors", "msg_warnings",
    ]
    with open(output_csv, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for r in rows:
            row_export = {k: r[k] for k in fieldnames}
            row_export["total_time_sec"] = r["total_time_sec"] if r["total_time_sec"] is not None else ""
            w.writerow(row_export)

    # Optional Markdown summary
    if output_md is not None:
        output_md.parent.mkdir(parents=True, exist_ok=True)
        completed = sum(1 for r in rows if r["status"] == "COMPLETED")
        aborted = sum(1 for r in rows if r["status"] == "ABORTED")
        no_sta = sum(1 for r in rows if r["status"] == "NO_STA")
        total_time_sum = sum(r["total_time_sec"] or 0 for r in rows)
        missing_dirs = sorted(set(expected_jobs) - set(result_jobs)) if include_expected else []

        lines = [
            "# Abaqus results review summary",
            "",
            f"- **Result directories scanned:** {len(result_jobs)}",
            f"- **Completed successfully:** {completed}",
            f"- **Aborted / not completed:** {aborted}",
            f"- **No .sta file:** {no_sta}",
            f"- **Total solver time (sec):** {total_time_sum:.2f}",
            "",
            "Validation performance (FEM vs Abaqus agreement) is produced by comparison scripts:",
            "- `deformation/deformation_comparison.py` → `deformation/deformation_plots/`",
            "- `section_forces/section_forces_comparison.py` → `section_forces/section_forces_plots/`",
            "- `grid_convergence_study/gci_richardson_abaqus_report.py` → `grid_convergence_study/gci_tables/`",
            "- See those directories for overlay plots and GCI/review CSVs.",
            "",
        ]
        if missing_dirs:
            lines.extend([
                "## Expected jobs with no result directory",
                "",
                "The following jobs exist under `jobs/` but have no `abaqus_results/job_XXXX_nN/` dir:",
                "",
            ])
            for j in missing_dirs:
                lines.append(f"- {j}")
            lines.append("")

        with open(output_md, "w", encoding="utf-8") as f:
            f.write("\n".join(lines))

    # Errors and inconsistencies log
    if output_log is not None:
        output_log.parent.mkdir(parents=True, exist_ok=True)
        jobs_with_issues = [j for j, entries in errors_per_job if entries]
        log_lines = [
            f"Generated by review_abaqus_results.py at {datetime.now().isoformat()}",
            "",
            f"Total jobs: {len(result_jobs)}. Jobs with at least one error or inconsistency: {len(jobs_with_issues)}.",
            f"List: {', '.join(jobs_with_issues) if jobs_with_issues else '(none)'}",
            "",
        ]
        for job_name, entries in errors_per_job:
            log_lines.append(f"## {job_name}")
            if not entries:
                log_lines.append("No errors or inconsistencies.")
            else:
                for source, message in entries:
                    log_lines.append(f"  {source} {message}")
            log_lines.append("")
        with open(output_log, "w", encoding="utf-8") as f:
            f.write("\n".join(log_lines))

    return 0


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Review abaqus_results directories and Abaqus run performance."
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help=f"Output CSV path (default: {DEFAULT_OUTPUT_CSV.relative_to(PROJECT_ROOT)}).",
    )
    parser.add_argument(
        "--md",
        type=Path,
        default=None,
        help="Optional Markdown summary path (default: output/abaqus_performance_summary.md if --expected).",
    )
    parser.add_argument(
        "--expected",
        action="store_true",
        help="Include expected job list from jobs/ and report missing result dirs.",
    )
    parser.add_argument(
        "--no-log-errors",
        action="store_true",
        help="Do not write the errors/inconsistencies log file.",
    )
    args = parser.parse_args()

    output_csv = args.output if args.output is not None else DEFAULT_OUTPUT_CSV
    if not output_csv.is_absolute():
        output_csv = PROJECT_ROOT / output_csv
    output_md = args.md
    if output_md is not None:
        if not output_md.is_absolute():
            output_md = PROJECT_ROOT / output_md
    else:
        output_md = PROJECT_ROOT / DEFAULT_OUTPUT_MD
    output_log = None if args.no_log_errors else (PROJECT_ROOT / DEFAULT_OUTPUT_LOG)

    return run_review(
        include_expected=args.expected,
        output_csv=output_csv,
        output_md=output_md,
        output_log=output_log,
    )


if __name__ == "__main__":
    sys.exit(main())

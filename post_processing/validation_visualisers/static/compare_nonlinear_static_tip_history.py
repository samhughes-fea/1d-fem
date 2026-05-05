from __future__ import annotations

import csv
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def benchmark_paths(
    results_root: str | Path,
    abaqus_root: str | Path,
    job_name: str,
    *,
    abaqus_reference_job_name: str | None = None,
) -> dict[str, str]:
    results_root = Path(results_root)
    abaqus_root = Path(abaqus_root)
    abaqus_job = abaqus_reference_job_name or job_name
    return {
        "fem_csv": str(results_root / "primary_results" / "nonlinear_static_validation" / f"{job_name}_tip_load_history.csv"),
        "abaqus_csv": str(abaqus_root / abaqus_job / "tip_load_history.csv"),
        "output_dir": str(results_root / "validation" / "nonlinear_static" / job_name),
    }


def infer_abaqus_reference_job_name(job_name: str) -> str:
    if job_name.endswith("_n64"):
        return job_name[:-4] + "_n500"
    return job_name


def read_tip_history_csv(path: str | Path) -> list[tuple[float, float]]:
    rows: list[tuple[float, float]] = []
    with open(path, "r", encoding="utf-8") as f:
        r = csv.reader(f)
        header = next(r, None)
        if not header:
            return rows
        normalized = [str(col).strip().lower() for col in header]
        try:
            load_factor_idx = normalized.index("load_factor")
            tip_displacement_idx = normalized.index("tip_displacement")
        except ValueError:
            load_factor_idx = 0
            tip_displacement_idx = 1
        for row in r:
            if not row:
                continue
            if len(row) <= max(load_factor_idx, tip_displacement_idx):
                continue
            rows.append((float(row[load_factor_idx]), float(row[tip_displacement_idx])))
    return rows


def _interpolate_displacement_at_load_factor(
    history: list[tuple[float, float]],
    target_load_factor: float,
) -> tuple[float, float]:
    if not history:
        raise ValueError("Cannot interpolate empty tip-history data")

    history = sorted(history, key=lambda pair: pair[0])
    first_lf, first_u = history[0]
    last_lf, last_u = history[-1]

    if target_load_factor <= first_lf:
        return first_lf, first_u
    if target_load_factor >= last_lf:
        return last_lf, last_u

    for (lf0, u0), (lf1, u1) in zip(history, history[1:]):
        if math.isclose(target_load_factor, lf0, rel_tol=0.0, abs_tol=1e-12):
            return lf0, u0
        if math.isclose(target_load_factor, lf1, rel_tol=0.0, abs_tol=1e-12):
            return lf1, u1
        if lf0 <= target_load_factor <= lf1:
            if math.isclose(lf1, lf0, rel_tol=0.0, abs_tol=1e-15):
                return lf0, 0.5 * (u0 + u1)
            t = (target_load_factor - lf0) / (lf1 - lf0)
            return target_load_factor, u0 + t * (u1 - u0)

    return last_lf, last_u


def compare_tip_histories(
    fem_csv: str | Path,
    abaqus_csv: str | Path,
) -> list[dict[str, float]]:
    fem = read_tip_history_csv(fem_csv)
    aba = read_tip_history_csv(abaqus_csv)
    out: list[dict[str, float]] = []
    for lf_f, u_f in fem:
        lf_a, u_a = _interpolate_displacement_at_load_factor(aba, lf_f)
        out.append(
            {
                "load_factor_fem": lf_f,
                "load_factor_abaqus": lf_a,
                "tip_displacement_fem": u_f,
                "tip_displacement_abaqus": u_a,
                "load_factor_alignment_error": abs(lf_f - lf_a),
                "abs_error": abs(u_f - u_a),
            }
        )
    return out


def write_tip_history_comparison_csv(rows: list[dict[str, float]], out_path: str | Path) -> None:
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow([
            "load_factor_fem",
            "load_factor_abaqus",
            "tip_displacement_fem",
            "tip_displacement_abaqus",
            "load_factor_alignment_error",
            "abs_error",
        ])
        for r in rows:
            w.writerow([
                r["load_factor_fem"],
                r["load_factor_abaqus"],
                r["tip_displacement_fem"],
                r["tip_displacement_abaqus"],
                r["load_factor_alignment_error"],
                r["abs_error"],
            ])


def write_tip_history_summary(rows: list[dict[str, float]], out_path: str | Path) -> None:
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    max_err = max((r["abs_error"] for r in rows), default=0.0)
    mean_err = sum((r["abs_error"] for r in rows), 0.0) / max(1, len(rows))
    max_load_alignment_err = max((r["load_factor_alignment_error"] for r in rows), default=0.0)
    with open(out_path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["metric", "value"])
        w.writerow(["num_rows", len(rows)])
        w.writerow(["max_load_factor_alignment_error", float(max_load_alignment_err)])
        w.writerow(["max_abs_error", float(max_err)])
        w.writerow(["mean_abs_error", float(mean_err)])


def plot_tip_history_comparison(rows: list[dict[str, float]], out_path: str | Path) -> None:
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    lf_f = [r["load_factor_fem"] for r in rows]
    u_f = [r["tip_displacement_fem"] for r in rows]
    lf_a = [r["load_factor_abaqus"] for r in rows]
    u_a = [r["tip_displacement_abaqus"] for r in rows]
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(lf_f, u_f, "o-", label="FEM")
    ax.plot(lf_a, u_a, "s--", label="Abaqus")
    ax.set_xlabel("Load factor")
    ax.set_ylabel("Tip displacement")
    ax.legend()
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def run_nonlinear_static_tip_history_validation(
    *,
    fem_csv: str | Path,
    abaqus_csv: str | Path,
    output_dir: str | Path,
) -> dict[str, str]:
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    rows = compare_tip_histories(fem_csv, abaqus_csv)
    comp_csv = output_dir / "tip_history_comparison.csv"
    summary_csv = output_dir / "tip_history_summary.csv"
    plot_png = output_dir / "tip_history_comparison.png"
    write_tip_history_comparison_csv(rows, comp_csv)
    write_tip_history_summary(rows, summary_csv)
    plot_tip_history_comparison(rows, plot_png)
    return {
        "comparison_csv": str(comp_csv),
        "summary_csv": str(summary_csv),
        "plot_png": str(plot_png),
    }


def run_nonlinear_static_tip_history_validation_from_benchmark_paths(
    *,
    results_root: str | Path,
    abaqus_root: str | Path,
    job_name: str,
    abaqus_reference_job_name: str | None = None,
) -> dict[str, str | bool]:
    abaqus_job = abaqus_reference_job_name or infer_abaqus_reference_job_name(job_name)
    paths = benchmark_paths(results_root, abaqus_root, job_name, abaqus_reference_job_name=abaqus_job)
    fem_csv = Path(paths["fem_csv"])
    abaqus_csv = Path(paths["abaqus_csv"])
    output_dir = Path(paths["output_dir"])
    output_dir.mkdir(parents=True, exist_ok=True)
    status_file = output_dir / "status.txt"
    if not fem_csv.is_file() or not abaqus_csv.is_file():
        status_file.write_text("missing_reference\n", encoding="utf-8")
        return {
            "ready": False,
            "status_file": str(status_file),
            "fem_csv": str(fem_csv),
            "abaqus_csv": str(abaqus_csv),
            "abaqus_reference_job_name": abaqus_job,
        }
    out = run_nonlinear_static_tip_history_validation(
        fem_csv=fem_csv,
        abaqus_csv=abaqus_csv,
        output_dir=output_dir,
    )
    status_file.write_text("comparison_complete\n", encoding="utf-8")
    out["status_file"] = str(status_file)
    out["ready"] = True
    out["abaqus_reference_job_name"] = abaqus_job
    return out

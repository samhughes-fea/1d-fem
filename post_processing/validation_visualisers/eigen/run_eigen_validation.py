from __future__ import annotations

from pathlib import Path

from post_processing.validation_visualisers.abaqus.config import ABAQUS_RESULTS_DIR


def run_eigen_validation(results_root: str | Path) -> dict[str, str | bool]:
    results_root = Path(results_root)
    out_dir = results_root / "validation" / "eigen"
    out_dir.mkdir(parents=True, exist_ok=True)
    summary = out_dir / "eigen_validation_summary.txt"
    reference_dir = ABAQUS_RESULTS_DIR / "job_benchmark_eigen_cantilever"
    expected = ["eigen_frequencies.csv", "mode_shapes.csv"]
    present = [name for name in expected if (reference_dir / name).is_file()]
    summary.write_text(
        "status=scaffold\n"
        "family=eigen\n"
        f"reference_dir={reference_dir}\n"
        f"expected_files={','.join(expected)}\n"
        f"present_files={','.join(present)}\n",
        encoding="utf-8",
    )
    return {
        "family": "eigen",
        "ready": bool(present),
        "summary_file": str(summary),
        "reference_dir": str(reference_dir),
        "present_files": ",".join(present),
    }

from __future__ import annotations

from pathlib import Path

from post_processing.validation_visualisers.deformation.deformation_comparison import run_deformation_comparison
from post_processing.validation_visualisers.section_forces.section_forces_comparison import run_section_forces_comparison


def run_linear_static_validation(results_root: str | Path) -> dict[str, str | bool]:
    results_root = Path(results_root)
    out_dir = results_root / "validation" / "linear_static"
    out_dir.mkdir(parents=True, exist_ok=True)
    run_deformation_comparison()
    run_section_forces_comparison()
    summary = out_dir / "linear_static_validation_summary.txt"
    summary.write_text("status=scaffold\nfamily=linear_static\n", encoding="utf-8")
    return {"family": "linear_static", "ready": False, "summary_file": str(summary)}

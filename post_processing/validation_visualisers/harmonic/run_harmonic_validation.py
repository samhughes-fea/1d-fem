from __future__ import annotations

from pathlib import Path


def run_harmonic_validation(results_root: str | Path) -> dict[str, str | bool]:
    results_root = Path(results_root)
    out_dir = results_root / "validation" / "harmonic"
    out_dir.mkdir(parents=True, exist_ok=True)
    summary = out_dir / "harmonic_validation_summary.txt"
    summary.write_text("status=scaffold\nfamily=harmonic\n", encoding="utf-8")
    return {"family": "harmonic", "ready": False, "summary_file": str(summary)}

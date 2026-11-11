"""Analytical reference for the nonlinear buckling imperfect-column benchmark.

Provides the Euler critical load factor for the cantilever (fixed-free) column,
derived from beam properties in the benchmark job directory. Used for the
analytical_only_pass check in run_nonlinear_buckling_validation when Abaqus
reference artifacts are absent.

For a fixed-free (cantilever) column the Euler critical load is:
    P_cr = π² EI_min / (4 L²)

where I_min is the minimum second moment of area (weak-axis bending).
"""
from __future__ import annotations

import math
from pathlib import Path


def load_euler_column_params(job_dir: str | Path) -> dict[str, float]:
    """Derive cantilever Euler column parameters from the benchmark job directory."""
    job_dir = Path(job_dir)
    E, _ = _parse_material(job_dir / "material.txt")
    I_y = _parse_weak_axis_I(job_dir / "section.txt")
    L = _parse_beam_length(job_dir / "grid.txt")
    P = _parse_axial_load(job_dir / "point_load.txt")
    P_cr = math.pi**2 * E * I_y / (4.0 * L**2)
    return {
        "L": L,
        "E": E,
        "I_y": I_y,
        "P_applied": P,
        "P_cr": P_cr,
        "lambda_cr": P_cr / P if P > 0.0 else float("inf"),
    }


def _parse_material(path: Path) -> tuple[float, float]:
    for line in path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line or line.startswith("[") or line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) >= 5:
            return float(parts[1]), float(parts[4])
    raise ValueError(f"Cannot parse material from {path}")


def _parse_weak_axis_I(path: Path) -> float:
    """Return I_y (weak axis, index 3) from first data row of section.txt."""
    for line in path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line or line.startswith("[") or line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) >= 4:
            return float(parts[3])
    raise ValueError(f"Cannot parse section from {path}")


def _parse_beam_length(path: Path) -> float:
    max_x = 0.0
    for line in path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line or line.startswith("[") or line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) >= 2:
            try:
                max_x = max(max_x, float(parts[1]))
            except ValueError:
                pass
    return max_x


def _parse_axial_load(path: Path) -> float:
    """Return |F_x| from the point load data row."""
    for line in path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line or line.startswith("[") or line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) >= 9:
            return abs(float(parts[3]))
    raise ValueError(f"Cannot parse axial load from {path}")

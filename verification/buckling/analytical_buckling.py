from __future__ import annotations

from dataclasses import dataclass
from math import pi
from pathlib import Path


@dataclass(frozen=True)
class EulerBucklingBenchmark:
    length_m: float
    youngs_modulus_pa: float
    second_moment_m4: float
    effective_length_factor: float
    reference_axial_load_n: float


def _data_lines(path: str | Path) -> list[list[str]]:
    lines: list[list[str]] = []
    for raw_line in Path(path).read_text(encoding="utf-8").splitlines():
        stripped = raw_line.strip()
        if not stripped or stripped.startswith("#") or stripped.startswith("["):
            continue
        lines.append(stripped.split())
    return lines


def load_euler_buckling_benchmark(job_dir: str | Path) -> EulerBucklingBenchmark:
    job_dir = Path(job_dir)
    grid_lines = _data_lines(job_dir / "grid.txt")
    material_lines = _data_lines(job_dir / "material.txt")
    section_lines = _data_lines(job_dir / "section.txt")
    point_load_lines = _data_lines(job_dir / "point_load.txt")
    x_coords = [float(parts[1]) for parts in grid_lines]
    length_m = max(x_coords) - min(x_coords)
    youngs_modulus_pa = float(material_lines[0][1])
    second_moment_m4 = float(section_lines[0][4])
    reference_axial_load_n = abs(float(point_load_lines[0][3]))
    return EulerBucklingBenchmark(
        length_m=length_m,
        youngs_modulus_pa=youngs_modulus_pa,
        second_moment_m4=second_moment_m4,
        effective_length_factor=2.0,
        reference_axial_load_n=reference_axial_load_n,
    )


def euler_critical_load_n(benchmark: EulerBucklingBenchmark) -> float:
    return (pi**2 * benchmark.youngs_modulus_pa * benchmark.second_moment_m4) / (
        (benchmark.effective_length_factor * benchmark.length_m) ** 2
    )


def euler_buckling_load_factor(benchmark: EulerBucklingBenchmark) -> float:
    return euler_critical_load_n(benchmark) / benchmark.reference_axial_load_n

from __future__ import annotations

from dataclasses import dataclass
from math import pi, sqrt
from pathlib import Path


@dataclass(frozen=True)
class CantileverEigenBenchmark:
    length_m: float
    youngs_modulus_pa: float
    density_kg_per_m3: float
    area_m2: float
    second_moment_m4: float


CANTILEVER_BETA_ROOTS = (
    1.875104068711961,
    4.694091132974174,
    7.854757438237612,
    10.995540734875467,
)


def _data_lines(path: str | Path) -> list[list[str]]:
    lines: list[list[str]] = []
    for raw_line in Path(path).read_text(encoding="utf-8").splitlines():
        stripped = raw_line.strip()
        if not stripped or stripped.startswith("#") or stripped.startswith("["):
            continue
        lines.append(stripped.split())
    return lines


def load_cantilever_eigen_benchmark(job_dir: str | Path) -> CantileverEigenBenchmark:
    job_dir = Path(job_dir)
    grid_lines = _data_lines(job_dir / "grid.txt")
    material_lines = _data_lines(job_dir / "material.txt")
    section_lines = _data_lines(job_dir / "section.txt")
    x_coords = [float(parts[1]) for parts in grid_lines]
    length_m = max(x_coords) - min(x_coords)
    youngs_modulus_pa = float(material_lines[0][1])
    density_kg_per_m3 = float(material_lines[0][4])
    area_m2 = float(section_lines[0][1])
    second_moment_m4 = float(section_lines[0][3])
    return CantileverEigenBenchmark(
        length_m=length_m,
        youngs_modulus_pa=youngs_modulus_pa,
        density_kg_per_m3=density_kg_per_m3,
        area_m2=area_m2,
        second_moment_m4=second_moment_m4,
    )


def cantilever_eigen_frequencies_hz(benchmark: CantileverEigenBenchmark, *, num_modes: int = 4) -> list[float]:
    scale = sqrt(
        benchmark.youngs_modulus_pa * benchmark.second_moment_m4
        / (benchmark.density_kg_per_m3 * benchmark.area_m2)
    )
    length_sq = benchmark.length_m**2
    frequencies: list[float] = []
    for beta in CANTILEVER_BETA_ROOTS[:num_modes]:
        omega_n = (beta**2 / length_sq) * scale
        frequencies.append(omega_n / (2.0 * pi))
    return frequencies

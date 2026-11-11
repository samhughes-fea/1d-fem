from __future__ import annotations

from dataclasses import dataclass
from math import pi
from pathlib import Path

# Project-wide shear correction factor: 5/6 used throughout the FEM codebase.
SHEAR_CORRECTION_FACTOR = 5.0 / 6.0


@dataclass(frozen=True)
class TimoshenkoColumnBenchmark:
    length_m: float
    youngs_modulus_pa: float
    shear_modulus_pa: float
    area_m2: float
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


def load_timoshenko_column_benchmark(job_dir: str | Path) -> TimoshenkoColumnBenchmark:
    job_dir = Path(job_dir)
    grid_lines = _data_lines(job_dir / "grid.txt")
    material_lines = _data_lines(job_dir / "material.txt")
    section_lines = _data_lines(job_dir / "section.txt")
    point_load_lines = _data_lines(job_dir / "point_load.txt")
    x_coords = [float(parts[1]) for parts in grid_lines]
    length_m = max(x_coords) - min(x_coords)
    youngs_modulus_pa = float(material_lines[0][1])
    shear_modulus_pa = float(material_lines[0][2])
    area_m2 = float(section_lines[0][1])
    second_moment_m4 = float(section_lines[0][4])
    reference_axial_load_n = abs(float(point_load_lines[0][3]))
    return TimoshenkoColumnBenchmark(
        length_m=length_m,
        youngs_modulus_pa=youngs_modulus_pa,
        shear_modulus_pa=shear_modulus_pa,
        area_m2=area_m2,
        second_moment_m4=second_moment_m4,
        effective_length_factor=2.0,
        reference_axial_load_n=reference_axial_load_n,
    )


def euler_critical_load_n(benchmark: TimoshenkoColumnBenchmark) -> float:
    """Euler–Bernoulli critical load for the column (ignoring shear deformation)."""
    return (pi ** 2 * benchmark.youngs_modulus_pa * benchmark.second_moment_m4) / (
        (benchmark.effective_length_factor * benchmark.length_m) ** 2
    )


def engesser_critical_load_n(benchmark: TimoshenkoColumnBenchmark) -> float:
    """Timoshenko/Engesser critical load including shear deformation.

    Pcr_T = Pcr_Euler / (1 + Pcr_Euler / (κ * A * G))

    The denominator accounts for the shear stiffness of the column: as P approaches
    the shear stiffness κAG the column becomes shear-critical rather than bending-critical.
    """
    p_euler = euler_critical_load_n(benchmark)
    shear_stiffness = SHEAR_CORRECTION_FACTOR * benchmark.area_m2 * benchmark.shear_modulus_pa
    return p_euler / (1.0 + p_euler / shear_stiffness)


def engesser_buckling_load_factor(benchmark: TimoshenkoColumnBenchmark) -> float:
    return engesser_critical_load_n(benchmark) / benchmark.reference_axial_load_n

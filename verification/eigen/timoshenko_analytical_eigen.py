from __future__ import annotations

from dataclasses import dataclass
from math import pi, sqrt
from pathlib import Path

from verification.eigen.analytical_eigen import CANTILEVER_BETA_ROOTS

# Project-wide shear correction factor: 5/6 is used throughout the FEM codebase
# for the benchmark sections (rectangular and compact I-profile equivalents).
SHEAR_CORRECTION_FACTOR = 5.0 / 6.0


@dataclass(frozen=True)
class TimoshenkoCantileverBenchmark:
    length_m: float
    youngs_modulus_pa: float
    shear_modulus_pa: float
    density_kg_per_m3: float
    area_m2: float
    second_moment_m4: float


def _data_lines(path: str | Path) -> list[list[str]]:
    lines: list[list[str]] = []
    for raw_line in Path(path).read_text(encoding="utf-8").splitlines():
        stripped = raw_line.strip()
        if not stripped or stripped.startswith("#") or stripped.startswith("["):
            continue
        lines.append(stripped.split())
    return lines


def load_timoshenko_cantilever_benchmark(job_dir: str | Path) -> TimoshenkoCantileverBenchmark:
    job_dir = Path(job_dir)
    grid_lines = _data_lines(job_dir / "grid.txt")
    material_lines = _data_lines(job_dir / "material.txt")
    section_lines = _data_lines(job_dir / "section.txt")
    x_coords = [float(parts[1]) for parts in grid_lines]
    length_m = max(x_coords) - min(x_coords)
    youngs_modulus_pa = float(material_lines[0][1])
    shear_modulus_pa = float(material_lines[0][2])
    density_kg_per_m3 = float(material_lines[0][4])
    area_m2 = float(section_lines[0][1])
    second_moment_m4 = float(section_lines[0][3])
    return TimoshenkoCantileverBenchmark(
        length_m=length_m,
        youngs_modulus_pa=youngs_modulus_pa,
        shear_modulus_pa=shear_modulus_pa,
        density_kg_per_m3=density_kg_per_m3,
        area_m2=area_m2,
        second_moment_m4=second_moment_m4,
    )


def timoshenko_eigen_frequencies_hz(
    benchmark: TimoshenkoCantileverBenchmark, *, num_modes: int = 4
) -> list[float]:
    """Cantilever natural frequencies with Timoshenko shear-deformation and rotary inertia corrections.

    Correction (Traill-Nash & Collar 1953 / Weaver, Timoshenko & Young):
        ω_n_T = ω_n_EB / sqrt(1 + Φ_n)
    where
        Φ_n = β_n² * (I / (A * L²)) * (1 + E / (κ * G))
    accounts for both rotary inertia (I/(A*L²) term) and shear deformation (E/(κG) factor).
    """
    L2 = benchmark.length_m ** 2
    eb_scale = sqrt(
        benchmark.youngs_modulus_pa * benchmark.second_moment_m4
        / (benchmark.density_kg_per_m3 * benchmark.area_m2)
    )
    i_over_al2 = benchmark.second_moment_m4 / (benchmark.area_m2 * L2)
    shear_factor = 1.0 + benchmark.youngs_modulus_pa / (
        SHEAR_CORRECTION_FACTOR * benchmark.shear_modulus_pa
    )
    frequencies: list[float] = []
    for beta in CANTILEVER_BETA_ROOTS[:num_modes]:
        omega_n_eb = (beta ** 2 / L2) * eb_scale
        phi = beta ** 2 * i_over_al2 * shear_factor
        omega_n_t = omega_n_eb / sqrt(1.0 + phi)
        frequencies.append(omega_n_t / (2.0 * pi))
    return frequencies

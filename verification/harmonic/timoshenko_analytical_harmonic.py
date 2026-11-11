from __future__ import annotations

from dataclasses import dataclass
from math import pi, sqrt
from pathlib import Path

from verification.eigen.analytical_eigen import CANTILEVER_BETA_ROOTS

# First cantilever beta root — mode 1 dominates the SDOF equivalence.
_BETA_1 = CANTILEVER_BETA_ROOTS[0]

# Project-wide shear correction factor: 5/6 used throughout the FEM codebase.
SHEAR_CORRECTION_FACTOR = 5.0 / 6.0


@dataclass(frozen=True)
class TimoshenkoHarmonicBenchmark:
    length_m: float
    youngs_modulus_pa: float
    shear_modulus_pa: float
    density_kg_per_m3: float
    area_m2: float
    second_moment_m4: float
    force_amplitude_n: float
    frequency_min_hz: float
    frequency_max_hz: float
    num_frequency_points: int
    damping_ratio: float


def _token_lines(path: str | Path) -> list[list[str]]:
    """Parse grid/material/section files: skip # comments and [header] lines, split tokens."""
    lines: list[list[str]] = []
    for raw_line in Path(path).read_text(encoding="utf-8").splitlines():
        stripped = raw_line.strip()
        if not stripped or stripped.startswith("#") or stripped.startswith("["):
            continue
        lines.append(stripped.split())
    return lines


def _text_lines(path: str | Path) -> list[str]:
    """Parse simulation_settings/point_load files: skip # comments, return stripped lines."""
    return [
        line.strip()
        for line in Path(path).read_text(encoding="utf-8").splitlines()
        if line.strip() and not line.strip().startswith("#")
    ]


def load_timoshenko_harmonic_benchmark(job_dir: str | Path) -> TimoshenkoHarmonicBenchmark:
    job_dir = Path(job_dir)
    grid_lines = _token_lines(job_dir / "grid.txt")
    material_lines = _token_lines(job_dir / "material.txt")
    section_lines = _token_lines(job_dir / "section.txt")
    sim_lines = _text_lines(job_dir / "simulation_settings.txt")
    point_load_lines = _text_lines(job_dir / "point_load.txt")
    x_coords = [float(parts[1]) for parts in grid_lines]
    length_m = max(x_coords) - min(x_coords)
    youngs_modulus_pa = float(material_lines[0][1])
    shear_modulus_pa = float(material_lines[0][2])
    density_kg_per_m3 = float(material_lines[0][4])
    area_m2 = float(section_lines[0][1])
    second_moment_m4 = float(section_lines[0][3])
    cfg: dict[str, str] = {}
    for line in sim_lines:
        if "=" in line:
            key, value = line.split("=", 1)
            cfg[key.strip()] = value.strip()
    load_tokens = point_load_lines[-1].split()
    force_amplitude_n = abs(float(load_tokens[4]))
    return TimoshenkoHarmonicBenchmark(
        length_m=length_m,
        youngs_modulus_pa=youngs_modulus_pa,
        shear_modulus_pa=shear_modulus_pa,
        density_kg_per_m3=density_kg_per_m3,
        area_m2=area_m2,
        second_moment_m4=second_moment_m4,
        force_amplitude_n=force_amplitude_n,
        frequency_min_hz=float(cfg["frequency_min_hz"]),
        frequency_max_hz=float(cfg["frequency_max_hz"]),
        num_frequency_points=int(cfg["num_frequency_points"]),
        damping_ratio=float(cfg["modal_damping_ratio"]),
    )


def _timoshenko_sdof_params(
    benchmark: TimoshenkoHarmonicBenchmark,
) -> tuple[float, float, float]:
    """Return (stiffness_N_per_m, mass_kg, damping_N_s_per_m) for the Timoshenko SDOF.

    Timoshenko tip stiffness for a cantilever:
        1/k_T = L³/(3EI) + L/(κAG)

    Mode-1 natural frequency with shear + rotary inertia correction:
        ω_n_T = ω_n_EB / sqrt(1 + β₁² * (I/(AL²)) * (1 + E/(κG)))

    Effective modal mass consistent with the SDOF model:
        m_eff = k_T / ω_n_T²
    """
    L = benchmark.length_m
    EI = benchmark.youngs_modulus_pa * benchmark.second_moment_m4
    kAG = SHEAR_CORRECTION_FACTOR * benchmark.area_m2 * benchmark.shear_modulus_pa
    rho_A = benchmark.density_kg_per_m3 * benchmark.area_m2

    stiffness = 1.0 / (L ** 3 / (3.0 * EI) + L / kAG)

    omega_n_eb = (_BETA_1 ** 2 / L ** 2) * sqrt(EI / rho_A)
    phi = _BETA_1 ** 2 * (benchmark.second_moment_m4 / (benchmark.area_m2 * L ** 2)) * (
        1.0 + benchmark.youngs_modulus_pa / (SHEAR_CORRECTION_FACTOR * benchmark.shear_modulus_pa)
    )
    omega_n_t = omega_n_eb / sqrt(1.0 + phi)

    mass = stiffness / omega_n_t ** 2
    damping = 2.0 * benchmark.damping_ratio * mass * omega_n_t
    return stiffness, mass, damping


def timoshenko_sdof_natural_frequency_hz(benchmark: TimoshenkoHarmonicBenchmark) -> float:
    stiffness, mass, _ = _timoshenko_sdof_params(benchmark)
    return sqrt(stiffness / mass) / (2.0 * pi)


def _frequency_grid_hz(benchmark: TimoshenkoHarmonicBenchmark) -> list[float]:
    if benchmark.num_frequency_points == 1:
        return [benchmark.frequency_min_hz]
    step = (benchmark.frequency_max_hz - benchmark.frequency_min_hz) / (
        benchmark.num_frequency_points - 1
    )
    return [benchmark.frequency_min_hz + i * step for i in range(benchmark.num_frequency_points)]


def timoshenko_sdof_frequency_response(
    benchmark: TimoshenkoHarmonicBenchmark,
) -> list[tuple[float, float, float]]:
    """SDOF frequency response using Timoshenko-corrected stiffness, natural frequency, and mass.

    Returns list of (frequency_hz, response_real, response_imag).
    """
    stiffness, mass, damping = _timoshenko_sdof_params(benchmark)
    rows: list[tuple[float, float, float]] = []
    for frequency_hz in _frequency_grid_hz(benchmark):
        omega = 2.0 * pi * frequency_hz
        real_impedance = stiffness - mass * omega ** 2
        imag_impedance = damping * omega
        denom = real_impedance ** 2 + imag_impedance ** 2
        response_real = benchmark.force_amplitude_n * real_impedance / denom
        response_imag = -benchmark.force_amplitude_n * imag_impedance / denom
        rows.append((frequency_hz, response_real, response_imag))
    return rows

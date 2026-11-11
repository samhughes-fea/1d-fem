from __future__ import annotations

from dataclasses import dataclass
from math import pi
from pathlib import Path


@dataclass(frozen=True)
class HarmonicSdofBenchmark:
    force_amplitude_n: float
    frequency_min_hz: float
    frequency_max_hz: float
    num_frequency_points: int
    natural_frequency_hz: float
    damping_ratio: float
    static_gain_m_per_n: float


def _data_lines(path: str | Path) -> list[str]:
    return [
        line.strip()
        for line in Path(path).read_text(encoding="utf-8").splitlines()
        if line.strip() and not line.strip().startswith("#")
    ]


def load_harmonic_sdof_benchmark(job_dir: str | Path) -> HarmonicSdofBenchmark:
    job_dir = Path(job_dir)
    sim_lines = _data_lines(job_dir / "simulation_settings.txt")
    point_load_lines = _data_lines(job_dir / "point_load.txt")
    cfg: dict[str, str] = {}
    for line in sim_lines:
        if "=" in line:
            key, value = line.split("=", 1)
            cfg[key.strip()] = value.strip()
    load_tokens = point_load_lines[-1].split()
    force_amplitude_n = abs(float(load_tokens[4]))
    return HarmonicSdofBenchmark(
        force_amplitude_n=force_amplitude_n,
        frequency_min_hz=float(cfg["frequency_min_hz"]),
        frequency_max_hz=float(cfg["frequency_max_hz"]),
        num_frequency_points=int(cfg["num_frequency_points"]),
        natural_frequency_hz=10.0,
        damping_ratio=float(cfg["modal_damping_ratio"]),
        static_gain_m_per_n=2.5e-6,
    )


def harmonic_sdof_frequency_grid_hz(benchmark: HarmonicSdofBenchmark) -> list[float]:
    if benchmark.num_frequency_points == 1:
        return [benchmark.frequency_min_hz]
    step = (benchmark.frequency_max_hz - benchmark.frequency_min_hz) / (benchmark.num_frequency_points - 1)
    return [benchmark.frequency_min_hz + i * step for i in range(benchmark.num_frequency_points)]


def harmonic_sdof_frequency_response(benchmark: HarmonicSdofBenchmark) -> list[tuple[float, float, float]]:
    rows: list[tuple[float, float, float]] = []
    omega_n = 2.0 * pi * benchmark.natural_frequency_hz
    stiffness = benchmark.force_amplitude_n / benchmark.static_gain_m_per_n
    mass = stiffness / (omega_n**2)
    damping = 2.0 * benchmark.damping_ratio * mass * omega_n
    for frequency_hz in harmonic_sdof_frequency_grid_hz(benchmark):
        omega = 2.0 * pi * frequency_hz
        real = stiffness - mass * omega**2
        imag = damping * omega
        denom = real**2 + imag**2
        response_real = benchmark.force_amplitude_n * real / denom
        response_imag = -benchmark.force_amplitude_n * imag / denom
        rows.append((frequency_hz, response_real, response_imag))
    return rows

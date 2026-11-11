"""SDOF analytical reference for the cantilever transient benchmark.

Models the cantilever beam as an undamped SDOF oscillator using first-mode
effective stiffness (k = 3EI/L³) and effective mass (m_eff = 33/140·ρAL).
The forcing is sinusoidal at the tip (F_y·sin(Ω·t)).

Intended for use in the analytical_only_pass path of run_transient_validation
when Abaqus reference artifacts are not yet available.
"""
from __future__ import annotations

import math
from pathlib import Path


def load_transient_sdof_params(job_dir: str | Path) -> dict[str, float]:
    """Derive SDOF parameters from the cantilever transient benchmark job directory."""
    job_dir = Path(job_dir)
    E, rho = _parse_material(job_dir / "material.txt")
    A, I_z = _parse_section(job_dir / "section.txt")
    L = _parse_beam_length(job_dir / "grid.txt")
    freq_hz = _parse_force_frequency_hz(job_dir / "simulation_settings.txt")
    F0 = _parse_force_amplitude(job_dir / "point_load.txt")
    EI = E * I_z
    k = 3.0 * EI / L**3
    m_eff = (33.0 / 140.0) * rho * A * L
    omega_n = math.sqrt(k / m_eff)
    omega_force = 2.0 * math.pi * freq_hz
    return {
        "L": L,
        "EI": EI,
        "k": k,
        "m_eff": m_eff,
        "omega_n": omega_n,
        "omega_force": omega_force,
        "F0": F0,
    }


def forced_sdof_steady_state_amplitude(params: dict[str, float]) -> float:
    """Steady-state tip displacement amplitude for undamped SDOF sinusoidal forcing.

    Uses the dynamic amplification factor D = 1 / |1 - r²| where r = Ω/ω_n.
    Returns the absolute static deflection multiplied by D.
    """
    r = params["omega_force"] / params["omega_n"]
    denom = abs(1.0 - r**2)
    if denom < 1e-12:
        return float("inf")
    static_deflection = params["F0"] / params["k"]
    return abs(static_deflection / denom)


def _parse_material(path: Path) -> tuple[float, float]:
    for line in path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line or line.startswith("["):
            continue
        parts = line.split()
        if len(parts) >= 5:
            return float(parts[1]), float(parts[4])
    raise ValueError(f"Cannot parse material from {path}")


def _parse_section(path: Path) -> tuple[float, float]:
    for line in path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line or line.startswith("["):
            continue
        parts = line.split()
        if len(parts) >= 5:
            return float(parts[1]), float(parts[4])
    raise ValueError(f"Cannot parse section from {path}")


def _parse_beam_length(path: Path) -> float:
    max_x = 0.0
    for line in path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line or line.startswith("["):
            continue
        parts = line.split()
        if len(parts) >= 2:
            try:
                max_x = max(max_x, float(parts[1]))
            except ValueError:
                pass
    return max_x


def _parse_force_frequency_hz(path: Path) -> float:
    for line in path.read_text(encoding="utf-8").splitlines():
        if "force_analytic_frequency_hz" in line and "=" in line:
            return float(line.split("=")[1].strip())
    raise ValueError(f"Cannot find force_analytic_frequency_hz in {path}")


def _parse_force_amplitude(path: Path) -> float:
    for line in path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line or line.startswith("["):
            continue
        parts = line.split()
        if len(parts) >= 9:
            return abs(float(parts[4]))
    raise ValueError(f"Cannot parse force amplitude from {path}")

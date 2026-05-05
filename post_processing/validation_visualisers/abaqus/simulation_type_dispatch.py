from __future__ import annotations

from typing import Any, Dict

from pre_processing.parsing.simulation_settings_parser import parse_simulation_settings


def parse_validation_simulation_type(job_dir: str) -> str:
    settings = parse_simulation_settings(f"{job_dir}/simulation_settings.txt")
    return str(settings.get("type", "static")).strip().lower()


def build_validation_dispatch_payload(job_dir: str) -> Dict[str, Any]:
    settings = parse_simulation_settings(f"{job_dir}/simulation_settings.txt")
    simulation_type = str(settings.get("type", "static")).strip().lower()
    return {
        "simulation_type": simulation_type,
        "simulation_settings": settings,
        "simulation_settings_path": f"{job_dir}/simulation_settings.txt",
    }

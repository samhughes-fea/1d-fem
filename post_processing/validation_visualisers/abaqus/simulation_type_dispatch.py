from __future__ import annotations

from typing import Any, Dict

from pre_processing.parsing.simulation_settings_parser import parse_simulation_settings


def _artifact_contract_for_simulation_type(simulation_type: str) -> dict:
    sim = str(simulation_type).strip().lower()
    if sim == "transient":
        return {
            "contract_name": "transient_reference",
            "expected_files": ["U_global.csv", "rotation_source.txt", "transient_reference_contract.txt"],
        }
    if sim == "harmonic":
        return {
            "contract_name": "harmonic_reference",
            "expected_files": ["U_global.csv", "rotation_source.txt", "frequency_response.csv"],
        }
    if sim == "eigen":
        return {
            "contract_name": "eigen_reference",
            "expected_files": ["eigen_frequencies.csv", "mode_shapes.csv"],
        }
    if sim == "buckling":
        return {
            "contract_name": "linear_buckling_reference",
            "expected_files": ["buckling_load_factors.csv", "buckling_mode_shapes.csv"],
        }
    return {
        "contract_name": "static_reference",
        "expected_files": ["U_global.csv", "rotation_source.txt"],
    }


def parse_validation_simulation_type(job_dir: str) -> str:
    settings = parse_simulation_settings(f"{job_dir}/simulation_settings.txt")
    return str(settings.get("type", "static")).strip().lower()


def build_validation_dispatch_payload(job_dir: str) -> Dict[str, Any]:
    settings = parse_simulation_settings(f"{job_dir}/simulation_settings.txt")
    simulation_type = str(settings.get("type", "static")).strip().lower()
    contract = _artifact_contract_for_simulation_type(simulation_type)
    return {
        "simulation_type": simulation_type,
        "simulation_settings": settings,
        "simulation_settings_path": f"{job_dir}/simulation_settings.txt",
        "artifact_contract": contract,
    }

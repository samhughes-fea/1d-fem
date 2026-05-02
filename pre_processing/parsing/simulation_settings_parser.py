# pre_processing\parsing\simulation_settings_parser.py

import os
import logging
import re

VALID_SIMULATION_TYPES = {"static", "static_nonlinear", "modal", "dynamic"}

def _get_defaults():
    """Return default values for all configuration sections."""
    return {
        "type": "static",
        "solver": {
            "type": "cg",
            "tolerance": 1e-6,
            "max_iterations": 1000,
            "restart": 20,
            "ilu_drop_tol": 1e-6,
            "ilu_fill_factor": 1.0,
            "disable_scaling": False,
        },
        "condensation": {
            "base_tol": 1e-12,
        },
        "parallel": {
            "num_processes": "auto",
            "enable_parallel_instantiation": True,
            "enable_parallel_computation": True,
        },
        "modal": {
            "num_modes": 10,
        },
        "dynamic": {
            "time_step": 0.001,
            "end_time": 1.0,
            "scheme": "newmark",
        },
        "newton": {
            "tolerance": 1e-8,
            "max_iterations": 50,
            "tolerance_delta_u": 1e-10,
            "relative_tolerance": None,
            "relative_reference": "first_residual",
        },
        "nonlinear": {
            "num_increments": 1,
            "load_factors": None,
            "line_search": False,
            "line_search_max_backtracks": 6,
            "line_search_shrink": 0.5,
        },
    }

def _validate_parallel_config(parallel_config):
    """Validate parallel configuration parameters."""
    if "num_processes" in parallel_config:
        num_proc = parallel_config["num_processes"]
        if num_proc != "auto" and not isinstance(num_proc, (int, str)):
            raise ValueError("num_processes must be 'auto' or an integer")
        if isinstance(num_proc, str) and num_proc != "auto":
            try:
                num_proc = int(num_proc)
                if num_proc < 1:
                    raise ValueError("num_processes must be >= 1")
            except ValueError:
                raise ValueError(f"Invalid num_processes value: {num_proc}")
        elif isinstance(num_proc, int) and num_proc < 1:
            raise ValueError("num_processes must be >= 1")

    for key in ["enable_parallel_instantiation", "enable_parallel_computation"]:
        if key in parallel_config:
            val = parallel_config[key]
            if isinstance(val, str):
                parallel_config[key] = val.lower() in ("true", "1", "yes", "on")
            elif not isinstance(val, bool):
                raise ValueError(f"{key} must be a boolean or boolean-like string")

def _parse_key_value(line):
    """Parse a key = value line, handling case insensitivity."""
    if "=" not in line:
        return None, None
    parts = line.split("=", 1)
    key = parts[0].strip().lower()
    value = parts[1].strip()
    return key, value

def _convert_value(value, expected_type):
    """Convert string value to expected type."""
    if expected_type == bool:
        if isinstance(value, bool):
            return value
        return value.lower() in ("true", "1", "yes", "on")
    elif expected_type == int:
        return int(float(value))  # Handle "1e6" style numbers
    elif expected_type == float:
        return float(value)
    else:
        return value

def parse_simulation_settings(file_path):
    """
    Parses simulation settings from a structured simulation_settings.txt file.

    Supports multiple sections:
    - [Simulation][Type]: simulation type (static, modal, dynamic)
    - [Solver]: solver configuration (type, tolerance, max_iterations, etc.)
    - [Condensation]: condensation configuration (base_tol)
    - [Parallel]: parallel processing configuration

    Parameters
    ----------
    file_path : str
        Path to the simulation settings file.

    Returns
    -------
    dict
        Dictionary with keys: 'type', 'solver', 'condensation', 'parallel'
        All sections have defaults if not specified.

    Raises
    ------
    FileNotFoundError
        If the file does not exist.

    ValueError
        If the simulation type is missing or unrecognized, or if config values are invalid.

    Example
    -------
    File: simulation_settings.txt
        [Simulation]
        [Type]
        Static
        
        [Solver]
        type = cg
        tolerance = 1e-6
        max_iterations = 1000
        
        [Parallel]
        num_processes = auto
        enable_parallel_computation = true
    """
    if not os.path.exists(file_path):
        logging.error(f"Simulation settings file not found: {file_path}")
        raise FileNotFoundError(f"{file_path} not found")

    defaults = _get_defaults()
    simulation_settings = {
        "type": defaults["type"],
        "solver": defaults["solver"].copy(),
        "condensation": defaults["condensation"].copy(),
        "parallel": defaults["parallel"].copy(),
        "modal": defaults["modal"].copy(),
        "dynamic": defaults["dynamic"].copy(),
        "newton": defaults["newton"].copy(),
        "nonlinear": defaults["nonlinear"].copy(),
    }

    current_section = None
    type_found = False

    with open(file_path, 'r') as f:
        for line_number, raw_line in enumerate(f, 1):
            # Remove inline comments and strip whitespace
            line = raw_line.split("#")[0].strip()
            if not line:
                continue

            # Detect section headers (case-insensitive)
            section_match = re.match(r"\[(.+?)\]", line, re.IGNORECASE)
            if section_match:
                section_name = section_match.group(1).strip().lower()
                
                if section_name == "simulation":
                    current_section = "simulation"
                elif section_name == "type":
                    if current_section == "simulation":
                        current_section = "type"
                    else:
                        current_section = None
                elif section_name == "solver":
                    current_section = "solver"
                elif section_name == "condensation":
                    current_section = "condensation"
                elif section_name == "parallel":
                    current_section = "parallel"
                elif section_name == "modal":
                    current_section = "modal"
                elif section_name == "dynamic":
                    current_section = "dynamic"
                elif section_name == "newton":
                    current_section = "newton"
                elif section_name == "nonlinear":
                    current_section = "nonlinear"
                else:
                    current_section = None
                continue

            # Parse content based on current section
            if current_section == "type":
                sim_type = line.lower()
                if sim_type not in VALID_SIMULATION_TYPES:
                    logging.error(f"[Simulation] Line {line_number}: Invalid simulation type '{line}'. "
                                  f"Expected one of {list(VALID_SIMULATION_TYPES)}.")
                    raise ValueError(f"Invalid simulation type: '{line}'")
                simulation_settings["type"] = sim_type
                type_found = True
                current_section = None  # Reset after reading type
            
            elif current_section == "solver":
                key, value = _parse_key_value(line)
                if key and value:
                    if key == "type":
                        simulation_settings["solver"]["type"] = value.lower()
                    elif key == "tolerance":
                        simulation_settings["solver"]["tolerance"] = _convert_value(value, float)
                    elif key == "max_iterations":
                        simulation_settings["solver"]["max_iterations"] = _convert_value(value, int)
                    elif key == "restart":
                        simulation_settings["solver"]["restart"] = _convert_value(value, int)
                    elif key == "ilu_drop_tol":
                        simulation_settings["solver"]["ilu_drop_tol"] = _convert_value(value, float)
                    elif key == "ilu_fill_factor":
                        simulation_settings["solver"]["ilu_fill_factor"] = _convert_value(value, float)
                    elif key == "disable_scaling":
                        simulation_settings["solver"]["disable_scaling"] = _convert_value(value, bool)
            
            elif current_section == "condensation":
                key, value = _parse_key_value(line)
                if key and value:
                    if key == "base_tol":
                        simulation_settings["condensation"]["base_tol"] = _convert_value(value, float)
            
            elif current_section == "parallel":
                key, value = _parse_key_value(line)
                if key and value:
                    if key == "num_processes":
                        simulation_settings["parallel"]["num_processes"] = value.lower() if value.lower() == "auto" else _convert_value(value, int)
                    elif key == "enable_parallel_instantiation":
                        simulation_settings["parallel"]["enable_parallel_instantiation"] = _convert_value(value, bool)
                    elif key == "enable_parallel_computation":
                        simulation_settings["parallel"]["enable_parallel_computation"] = _convert_value(value, bool)

            elif current_section == "modal":
                key, value = _parse_key_value(line)
                if key and value:
                    if key == "num_modes":
                        simulation_settings["modal"]["num_modes"] = _convert_value(value, int)

            elif current_section == "dynamic":
                key, value = _parse_key_value(line)
                if key and value:
                    if key == "time_step":
                        simulation_settings["dynamic"]["time_step"] = _convert_value(value, float)
                    elif key == "end_time":
                        simulation_settings["dynamic"]["end_time"] = _convert_value(value, float)
                    elif key == "scheme":
                        simulation_settings["dynamic"]["scheme"] = value.lower()

            elif current_section == "newton":
                key, value = _parse_key_value(line)
                if key and value:
                    if key == "tolerance":
                        simulation_settings["newton"]["tolerance"] = _convert_value(value, float)
                    elif key == "max_iterations":
                        simulation_settings["newton"]["max_iterations"] = _convert_value(value, int)
                    elif key == "tolerance_delta_u":
                        simulation_settings["newton"]["tolerance_delta_u"] = _convert_value(value, float)
                    elif key == "relative_tolerance":
                        simulation_settings["newton"]["relative_tolerance"] = _convert_value(value, float)
                    elif key == "relative_reference":
                        simulation_settings["newton"]["relative_reference"] = value.strip().lower()

            elif current_section == "nonlinear":
                key, value = _parse_key_value(line)
                if key and value:
                    if key == "num_increments":
                        simulation_settings["nonlinear"]["num_increments"] = _convert_value(value, int)
                    elif key == "load_factors":
                        parts = [p.strip() for p in value.split(",") if p.strip()]
                        simulation_settings["nonlinear"]["load_factors"] = [
                            float(p) for p in parts
                        ]
                    elif key == "line_search":
                        simulation_settings["nonlinear"]["line_search"] = _convert_value(value, bool)
                    elif key == "line_search_max_backtracks":
                        simulation_settings["nonlinear"]["line_search_max_backtracks"] = _convert_value(
                            value, int
                        )
                    elif key == "line_search_shrink":
                        simulation_settings["nonlinear"]["line_search_shrink"] = _convert_value(value, float)

    # Validate simulation type was found (backward compatibility: allow missing if defaults are acceptable)
    if not type_found and simulation_settings["type"] == defaults["type"]:
        logging.warning("Simulation type not specified, using default: 'static'")
    
    # Validate parallel config
    try:
        _validate_parallel_config(simulation_settings["parallel"])
    except ValueError as e:
        logging.error(f"Invalid parallel configuration: {e}")
        raise

    return simulation_settings

# ------------------------------------------------
# Standalone execution for direct testing
# ------------------------------------------------
if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    
    test_file = os.path.join("jobs", "base", "simulation_settings.txt")

    if not os.path.exists(test_file):
        logging.error(f"Test file '{test_file}' not found. Please ensure it exists before running.")
    else:
        try:
            simulation_settings = parse_simulation_settings(test_file)
            print("\n------------- Parsed Simulation Settings -------------\n")
            print(simulation_settings)
        except Exception as e:
            logging.error(f"❌ Error parsing simulation settings file: {e}")

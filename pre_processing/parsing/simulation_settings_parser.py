# pre_processing\parsing\simulation_settings_parser.py

import os
import logging
import re

from pre_processing.parsing.simulation_settings_resolution import (
    VALID_TYPE_LINE_INPUTS,
    finalize_simulation_settings,
)

VALID_SIMULATION_TYPES = VALID_TYPE_LINE_INPUTS

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
            # vibration | buckling — buckling uses linear static prestress then (K_L + λ K_σ) eigenproblem
            "analysis": "vibration",
            # linear_static: run reference linear solve with job loads for internal forces (recommended)
            "buckling_prestress": "linear_static",
            # Scale reference loads before prestress solve (default 1.0)
            "buckling_load_factor": 1.0,
            # Twin nonlinear mesh for nonlinear_static when element.txt lists linear beams only
            "buckling_nonlinear_prestress_twins": False,
        },
        "dynamic": {
            "time_step": 0.001,
            "end_time": 1.0,
            "scheme": "newmark",
        },
        # Optional: strain/stress recovery from FormulationResultSet for modal/dynamic snapshots.
        "post_processing": {
            "run_secondary_tertiary_modal": False,
            "run_secondary_tertiary_dynamic": False,
            "run_secondary_tertiary_harmonic": False,
            "harmonic_frequency_index": 0,
            "harmonic_secondary_tertiary_all_frequencies": False,
            "harmonic_secondary_tertiary_frequency_indices": None,
            "harmonic_secondary_tertiary_displacement_component": "real",
            "modal_mode_index": 0,
            "modal_amplitude": 1.0,
            # Buckling snapshot: eigenmode shape vs prestress displacement field.
            "buckling_displacement": "mode",
            "buckling_mode_index": 0,
            # Row index into U(t) from Newmark; negative counts from end (-1 = last step).
            "dynamic_time_index": -1,
            # Optional comma-separated list of row indices (overrides single index when set).
            "dynamic_time_indices": None,
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
            # Optional: CorotationalBeamElement3D — finite_difference | elastic_material
            "corotational_tangent_mode": None,
            # Optional: GeometricallyExactShearDeformableBeam3D — use parent TL strain hook only
            "gesdb_tl_fallback": False,
            # Optional: GeometricallyExactShearDeformableBeam3D — tl_locked (default) | native
            "gesdb_kernel": "tl_locked",
        },
        # Taxonomy §1–§5 sections (optional ``enabled`` flags; see finalize_simulation_settings).
        "static": {
            "enabled": False,
        },
        "eigen": {
            "enabled": False,
            "num_modes": 10,
            "fixed_node_id": None,
            "dense_threshold": 512,
        },
        "transient": {
            "enabled": False,
            "time_step": 0.001,
            "end_time": 1.0,
            "scheme": "newmark",
            "load_scale": 1.0,
            "load_ramp": False,
            "fixed_node_id": None,
            "force_time_series_file": None,
            "force_analytic": None,
            "force_analytic_amplitude": None,
            "force_analytic_frequency_hz": None,
            "force_analytic_phase_rad": 0.0,
            "force_analytic_t_start": None,
            "force_analytic_t_end": None,
            "rayleigh_alpha": None,
            "rayleigh_beta": None,
        },
        "harmonic": {
            "enabled": False,
            # §4 harmonic sweep — defaults applied at runtime if omitted (see harmonic_simulation.effective_harmonic_config)
            "frequency_min_hz": None,
            "frequency_max_hz": None,
            "num_frequency_points": None,
            "modal_damping_ratio": None,
            "rayleigh_alpha": None,
            "rayleigh_beta": None,
            "load_phase_rad": None,
            "parallel_frequency_sweep": None,
            "mp_damping_reference": None,
            "use_modal_superposition": None,
            "modal_superposition_num_modes": None,
            "prescribed_motion_phase_rad": None,
            "harmonic_linear_solver": None,
            "damping_zeta_table_file": None,
            "harmonic_modal_basis_dir": None,
            "harmonic_modal_basis_job_name": None,
            "fixed_node_id": None,
        },
        "buckling": {
            "enabled": False,
            "num_modes": 10,
            "buckling_prestress": "linear_static",
            "buckling_load_factor": 1.0,
            "buckling_nonlinear_prestress_twins": False,
            # When true, run_job dispatches to NonlinearBucklingSimulationRunner (MVP stub; not linear eigen).
            "nonlinear_buckling": False,
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
    - [Simulation][Type]: canonical types static, eigen, transient, harmonic, buckling (legacy aliases modal, dynamic, static_nonlinear)
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
        Dictionary with keys including ``type``, ``solver``, ``condensation``, ``parallel``,
        ``modal``, ``dynamic``, ``static``, ``eigen``, ``transient``, ``harmonic``, ``buckling``,
        ``newton``, ``nonlinear``, ``post_processing``.
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
        "post_processing": defaults["post_processing"].copy(),
        "static": defaults["static"].copy(),
        "eigen": defaults["eigen"].copy(),
        "transient": defaults["transient"].copy(),
        "harmonic": defaults["harmonic"].copy(),
        "buckling": defaults["buckling"].copy(),
        "_modal_section_in_input": False,
        "_dynamic_section_in_input": False,
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
                    simulation_settings["_modal_section_in_input"] = True
                elif section_name == "dynamic":
                    current_section = "dynamic"
                    simulation_settings["_dynamic_section_in_input"] = True
                elif section_name == "newton":
                    current_section = "newton"
                elif section_name == "nonlinear":
                    current_section = "nonlinear"
                elif section_name.replace(" ", "").replace("_", "") == "postprocessing":
                    current_section = "post_processing"
                elif section_name == "static":
                    current_section = "static_taxonomy"
                elif section_name == "eigen":
                    current_section = "eigen"
                elif section_name == "transient":
                    current_section = "transient"
                elif section_name == "harmonic":
                    current_section = "harmonic"
                elif section_name == "buckling":
                    current_section = "buckling_taxonomy"
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
                    elif key == "dense_threshold":
                        simulation_settings["modal"]["dense_threshold"] = _convert_value(value, int)
                    elif key == "analysis":
                        av = value.strip().lower()
                        if av not in ("vibration", "buckling"):
                            raise ValueError(f"modal.analysis must be 'vibration' or 'buckling', got {value!r}")
                        simulation_settings["modal"]["analysis"] = av
                    elif key == "buckling_prestress":
                        bp = value.strip().lower()
                        if bp not in ("linear_static", "nonlinear_static", "none"):
                            raise ValueError(
                                "modal.buckling_prestress must be 'linear_static', "
                                f"'nonlinear_static', or 'none', got {value!r}"
                            )
                        simulation_settings["modal"]["buckling_prestress"] = bp
                    elif key == "buckling_load_factor":
                        simulation_settings["modal"]["buckling_load_factor"] = _convert_value(value, float)
                    elif key == "buckling_nonlinear_prestress_twins":
                        simulation_settings["modal"]["buckling_nonlinear_prestress_twins"] = _convert_value(
                            value, bool
                        )

            elif current_section == "dynamic":
                key, value = _parse_key_value(line)
                if key and value:
                    if key == "time_step":
                        simulation_settings["dynamic"]["time_step"] = _convert_value(value, float)
                    elif key == "end_time":
                        simulation_settings["dynamic"]["end_time"] = _convert_value(value, float)
                    elif key == "scheme":
                        simulation_settings["dynamic"]["scheme"] = value.lower()
                    elif key == "load_scale":
                        simulation_settings["dynamic"]["load_scale"] = _convert_value(value, float)
                    elif key == "load_ramp":
                        simulation_settings["dynamic"]["load_ramp"] = _convert_value(value, bool)
                    elif key == "fixed_node_id":
                        simulation_settings["dynamic"]["fixed_node_id"] = _convert_value(value, int)
                    elif key == "force_time_series_file":
                        simulation_settings["dynamic"]["force_time_series_file"] = value.strip()
                    elif key == "force_analytic":
                        simulation_settings["dynamic"]["force_analytic"] = value.strip().lower()
                    elif key == "force_analytic_amplitude":
                        simulation_settings["dynamic"]["force_analytic_amplitude"] = _convert_value(
                            value, float
                        )
                    elif key == "force_analytic_frequency_hz":
                        simulation_settings["dynamic"]["force_analytic_frequency_hz"] = _convert_value(
                            value, float
                        )
                    elif key == "force_analytic_phase_rad":
                        simulation_settings["dynamic"]["force_analytic_phase_rad"] = _convert_value(
                            value, float
                        )
                    elif key == "force_analytic_t_start":
                        simulation_settings["dynamic"]["force_analytic_t_start"] = _convert_value(
                            value, float
                        )
                    elif key == "force_analytic_t_end":
                        simulation_settings["dynamic"]["force_analytic_t_end"] = _convert_value(
                            value, float
                        )
                    elif key in ("rayleigh_alpha", "rayleigh_m"):
                        simulation_settings["dynamic"]["rayleigh_alpha"] = _convert_value(value, float)
                    elif key in ("rayleigh_beta", "rayleigh_k"):
                        simulation_settings["dynamic"]["rayleigh_beta"] = _convert_value(value, float)

            elif current_section == "post_processing":
                key, value = _parse_key_value(line)
                if key and value:
                    pp = simulation_settings["post_processing"]
                    if key in (
                        "run_secondary_tertiary_modal",
                        "run_secondary_tertiary_eigen",
                        "run_secondary_tertiary_buckling",
                    ):
                        pp["run_secondary_tertiary_modal"] = _convert_value(value, bool)
                    elif key == "run_secondary_tertiary_dynamic":
                        pp["run_secondary_tertiary_dynamic"] = _convert_value(value, bool)
                    elif key == "run_secondary_tertiary_harmonic":
                        pp["run_secondary_tertiary_harmonic"] = _convert_value(value, bool)
                    elif key == "harmonic_frequency_index":
                        pp["harmonic_frequency_index"] = _convert_value(value, int)
                    elif key == "harmonic_secondary_tertiary_all_frequencies":
                        pp["harmonic_secondary_tertiary_all_frequencies"] = _convert_value(value, bool)
                    elif key == "harmonic_secondary_tertiary_frequency_indices":
                        pp["harmonic_secondary_tertiary_frequency_indices"] = [
                            int(x.strip()) for x in value.split(",") if x.strip()
                        ]
                    elif key == "harmonic_secondary_tertiary_displacement_component":
                        hc = value.strip().lower()
                        if hc not in ("real", "imag", "both"):
                            raise ValueError(
                                "post_processing.harmonic_secondary_tertiary_displacement_component "
                                "must be 'real', 'imag', or 'both', "
                                f"got {value!r}"
                            )
                        pp["harmonic_secondary_tertiary_displacement_component"] = hc
                    elif key == "modal_mode_index":
                        pp["modal_mode_index"] = _convert_value(value, int)
                    elif key == "modal_amplitude":
                        pp["modal_amplitude"] = _convert_value(value, float)
                    elif key == "buckling_displacement":
                        bd = value.strip().lower()
                        if bd not in ("mode", "prestress"):
                            raise ValueError(
                                "post_processing.buckling_displacement must be 'mode' or 'prestress', "
                                f"got {value!r}"
                            )
                        pp["buckling_displacement"] = bd
                    elif key == "buckling_mode_index":
                        pp["buckling_mode_index"] = _convert_value(value, int)
                    elif key == "dynamic_time_index":
                        pp["dynamic_time_index"] = _convert_value(value, int)
                    elif key == "dynamic_time_indices":
                        parts = [p.strip() for p in value.split(",") if p.strip()]
                        pp["dynamic_time_indices"] = [_convert_value(p, int) for p in parts]

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
                    elif key == "corotational_tangent_mode":
                        simulation_settings["nonlinear"]["corotational_tangent_mode"] = value.strip().lower()
                    elif key == "gesdb_tl_fallback":
                        simulation_settings["nonlinear"]["gesdb_tl_fallback"] = _convert_value(value, bool)
                    elif key == "gesdb_kernel":
                        gk = value.strip().lower().replace("-", "_")
                        if gk not in (
                            "tl_locked",
                            "tl",
                            "locked",
                            "chord_tl",
                            "native",
                            "engineering",
                            "gesdb_native",
                        ):
                            raise ValueError(
                                "nonlinear.gesdb_kernel must be 'tl_locked' or 'native', "
                                f"got {value!r}"
                            )
                        simulation_settings["nonlinear"]["gesdb_kernel"] = value.strip()

            elif current_section == "static_taxonomy":
                key, value = _parse_key_value(line)
                if key and value:
                    if key == "enabled":
                        simulation_settings["static"]["enabled"] = _convert_value(value, bool)

            elif current_section == "eigen":
                key, value = _parse_key_value(line)
                if key and value:
                    if key == "enabled":
                        simulation_settings["eigen"]["enabled"] = _convert_value(value, bool)
                    elif key == "num_modes":
                        simulation_settings["eigen"]["num_modes"] = _convert_value(value, int)
                    elif key == "fixed_node_id":
                        simulation_settings["eigen"]["fixed_node_id"] = _convert_value(value, int)
                    elif key == "dense_threshold":
                        simulation_settings["eigen"]["dense_threshold"] = _convert_value(value, int)

            elif current_section == "transient":
                key, value = _parse_key_value(line)
                if key and value:
                    tr = simulation_settings["transient"]
                    if key == "enabled":
                        tr["enabled"] = _convert_value(value, bool)
                    elif key == "time_step":
                        tr["time_step"] = _convert_value(value, float)
                    elif key == "end_time":
                        tr["end_time"] = _convert_value(value, float)
                    elif key == "scheme":
                        tr["scheme"] = value.lower()
                    elif key == "load_scale":
                        tr["load_scale"] = _convert_value(value, float)
                    elif key == "load_ramp":
                        tr["load_ramp"] = _convert_value(value, bool)
                    elif key == "fixed_node_id":
                        tr["fixed_node_id"] = _convert_value(value, int)
                    elif key == "force_time_series_file":
                        tr["force_time_series_file"] = value.strip()
                    elif key == "force_analytic":
                        tr["force_analytic"] = value.strip().lower()
                    elif key == "force_analytic_amplitude":
                        tr["force_analytic_amplitude"] = _convert_value(value, float)
                    elif key == "force_analytic_frequency_hz":
                        tr["force_analytic_frequency_hz"] = _convert_value(value, float)
                    elif key == "force_analytic_phase_rad":
                        tr["force_analytic_phase_rad"] = _convert_value(value, float)
                    elif key == "force_analytic_t_start":
                        tr["force_analytic_t_start"] = _convert_value(value, float)
                    elif key == "force_analytic_t_end":
                        tr["force_analytic_t_end"] = _convert_value(value, float)
                    elif key in ("rayleigh_alpha", "rayleigh_m"):
                        tr["rayleigh_alpha"] = _convert_value(value, float)
                    elif key in ("rayleigh_beta", "rayleigh_k"):
                        tr["rayleigh_beta"] = _convert_value(value, float)

            elif current_section == "harmonic":
                key, value = _parse_key_value(line)
                if key and value:
                    hm = simulation_settings["harmonic"]
                    if key == "enabled":
                        hm["enabled"] = _convert_value(value, bool)
                    elif key in ("frequency_min_hz", "frequency_min"):
                        hm["frequency_min_hz"] = _convert_value(value, float)
                    elif key in ("frequency_max_hz", "frequency_max"):
                        hm["frequency_max_hz"] = _convert_value(value, float)
                    elif key in ("num_frequency_points", "num_points"):
                        hm["num_frequency_points"] = _convert_value(value, int)
                    elif key in ("modal_damping_ratio", "damping_ratio"):
                        hm["modal_damping_ratio"] = _convert_value(value, float)
                    elif key in ("rayleigh_alpha", "rayleigh_m"):
                        hm["rayleigh_alpha"] = _convert_value(value, float)
                    elif key in ("rayleigh_beta", "rayleigh_k"):
                        hm["rayleigh_beta"] = _convert_value(value, float)
                    elif key in ("load_phase_rad", "load_phase"):
                        hm["load_phase_rad"] = _convert_value(value, float)
                    elif key == "parallel_frequency_sweep":
                        hm["parallel_frequency_sweep"] = _convert_value(value, bool)
                    elif key in ("mp_damping_reference", "mass_proportional_damping_reference"):
                        hm["mp_damping_reference"] = value.strip().lower()
                    elif key == "use_modal_superposition":
                        hm["use_modal_superposition"] = _convert_value(value, bool)
                    elif key in ("modal_superposition_num_modes", "modal_superposition_modes"):
                        hm["modal_superposition_num_modes"] = _convert_value(value, int)
                    elif key in ("prescribed_motion_phase_rad", "prescribed_motion_phase"):
                        hm["prescribed_motion_phase_rad"] = _convert_value(value, float)
                    elif key == "harmonic_linear_solver":
                        hm["harmonic_linear_solver"] = value.strip().lower()
                    elif key in ("damping_zeta_table_file", "zeta_table_file"):
                        hm["damping_zeta_table_file"] = value.strip()
                    elif key in ("harmonic_modal_basis_dir", "modal_basis_dir"):
                        hm["harmonic_modal_basis_dir"] = value.strip()
                    elif key in ("harmonic_modal_basis_job_name", "modal_basis_job_name"):
                        hm["harmonic_modal_basis_job_name"] = value.strip()
                    elif key == "fixed_node_id":
                        hm["fixed_node_id"] = _convert_value(value, int)

            elif current_section == "buckling_taxonomy":
                key, value = _parse_key_value(line)
                if key and value:
                    bk = simulation_settings["buckling"]
                    if key == "enabled":
                        bk["enabled"] = _convert_value(value, bool)
                    elif key == "num_modes":
                        bk["num_modes"] = _convert_value(value, int)
                    elif key == "buckling_prestress":
                        bp = value.strip().lower()
                        if bp not in ("linear_static", "nonlinear_static", "none"):
                            raise ValueError(
                                "buckling.buckling_prestress must be 'linear_static', "
                                f"'nonlinear_static', or 'none', got {value!r}"
                            )
                        bk["buckling_prestress"] = bp
                    elif key == "buckling_load_factor":
                        bk["buckling_load_factor"] = _convert_value(value, float)
                    elif key == "buckling_nonlinear_prestress_twins":
                        bk["buckling_nonlinear_prestress_twins"] = _convert_value(value, bool)
                    elif key in ("nonlinear_buckling", "use_nonlinear_buckling"):
                        bk["nonlinear_buckling"] = _convert_value(value, bool)

    # Validate parallel config
    try:
        _validate_parallel_config(simulation_settings["parallel"])
    except ValueError as e:
        logging.error(f"Invalid parallel configuration: {e}")
        raise

    finalize_simulation_settings(simulation_settings, type_line_explicit=type_found)

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

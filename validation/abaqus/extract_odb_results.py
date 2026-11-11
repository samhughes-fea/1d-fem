# post_processing/validation_visualisers/abaqus/extract_odb_results.py
"""
Extract displacement (and optionally section forces) from an Abaqus ODB to CSV.
Designed to be run from inside Abaqus Python (abaqus cae noGUI=extract_odb_results.py)
or called from a generated run script.

Usage (from Abaqus Python):
  abaqus cae noGUI=extract_odb_results.py -- odb_path=<path> out_dir=<path>

Or import and call extract_odb_to_csv(odb_path, out_dir).
"""
from __future__ import annotations

import os
import sys
from collections import defaultdict


def _debug_log(out_dir: str, message: str) -> None:
    print(message)
    try:
        with open(os.path.join(out_dir, "run_log.txt"), "a") as f:
            f.write(str(message) + "\n")
    except Exception:
        pass


def _write_empty_harmonic_export(out_dir: str) -> None:
    with open(os.path.join(out_dir, "frequency_response.csv"), "w") as f:
        f.write("frequency_hz,response_real,response_imag\n")


def _write_empty_eigen_exports(out_dir: str) -> None:
    with open(os.path.join(out_dir, "eigen_frequencies.csv"), "w") as f:
        f.write("mode_index,frequency_hz\n")
    with open(os.path.join(out_dir, "mode_shapes.csv"), "w") as f:
        f.write("mode_index,node_label,dof,value\n")


def _write_empty_buckling_exports(out_dir: str) -> None:
    with open(os.path.join(out_dir, "buckling_load_factors.csv"), "w") as f:
        f.write("mode_index,load_factor\n")
    with open(os.path.join(out_dir, "buckling_mode_shapes.csv"), "w") as f:
        f.write("mode_index,node_label,dof,value\n")


def _write_empty_transient_export(out_dir: str) -> None:
    with open(os.path.join(out_dir, "transient_reference.csv"), "w") as f:
        f.write("time,tip_displacement\n")


def _write_empty_nl_buckling_export(out_dir: str) -> None:
    with open(os.path.join(out_dir, "continuation_history.csv"), "w") as f:
        f.write("load_factor,tip_displacement,converged\n")


def _extract_buckling_load_factors(step, out_dir: str) -> bool:
    _debug_log(out_dir, "DEBUG buckling load factors: start")
    frames = _spectral_frames_with_nontrivial_mode_data(step)
    _debug_log(out_dir, f"DEBUG buckling load factors: frames={len(frames)}")
    if len(frames) == 0:
        return False
    rows: list[tuple[int, float]] = []
    for i, frame in enumerate(frames, start=1):
        load_factor = _spectral_frequency_like_value(frame, prefer_description=True)
        if load_factor is None:
            continue
        rows.append((i, float(load_factor)))
    if len(rows) == 0:
        return False
    csv_path = os.path.join(out_dir, "buckling_load_factors.csv")
    with open(csv_path, "w") as f:
        f.write("mode_index,load_factor\n")
        for mode_index, load_factor in rows:
            f.write(f"{mode_index},{load_factor}\n")
    return True


def _buckling_numeric_token_from_description(description: str) -> float | None:
    cleaned = description.replace(",", " ").replace("=", " ").split()
    for token in reversed(cleaned):
        try:
            return float(token)
        except (TypeError, ValueError):
            continue
    return None


def _spectral_frames_with_mode_data(step) -> list:
    frames = list(getattr(step, "frames", []))
    return [frame for frame in frames if _spectral_frame_has_mode_data(frame)]


def _spectral_frames_with_nontrivial_mode_data(step) -> list:
    frames = _spectral_frames_with_mode_data(step)
    nontrivial = [frame for frame in frames if _spectral_frame_has_nontrivial_mode_shape(frame)]
    return nontrivial if len(nontrivial) > 0 else frames


def _spectral_frame_has_mode_data(frame) -> bool:
    field_outputs = getattr(frame, "fieldOutputs", None)
    u_field = _field_output_get(field_outputs, "U")
    if u_field is not None:
        values = list(getattr(u_field, "values", []))
        if len(values) > 0:
            return True
    return _spectral_frequency_like_value(frame, prefer_description=True) is not None


def _spectral_frame_has_nontrivial_mode_shape(frame) -> bool:
    u_field = _field_output_get(getattr(frame, "fieldOutputs", None), "U")
    if u_field is None:
        return False
    for value in list(getattr(u_field, "values", [])):
        raw_data = getattr(value, "data", None)
        data = list(raw_data) if raw_data is not None else []
        for component in data[:3]:
            if abs(float(component)) > 0.0:
                return True
    return False


def _field_output_get(field_outputs, key: str):
    if field_outputs is None:
        return None
    getter = getattr(field_outputs, "get", None)
    if callable(getter):
        try:
            return getter(key)
        except Exception:
            pass
    try:
        return field_outputs[key]
    except Exception:
        return None


def _field_outputs_has_key(field_outputs, key: str) -> bool:
    if field_outputs is None:
        return False
    value = _field_output_get(field_outputs, key)
    return value is not None


def _spectral_frequency_like_value(frame, *, prefer_description: bool) -> float | None:
    description_value = None
    description = getattr(frame, "description", None)
    if isinstance(description, str):
        description_value = _buckling_numeric_token_from_description(description)
    if prefer_description:
        ordered_candidates = [
            description_value,
            getattr(frame, "eigenvalue", None),
            getattr(frame, "frequency", None),
            getattr(frame, "frameValue", None),
        ]
    else:
        ordered_candidates = [
            getattr(frame, "frequency", None),
            getattr(frame, "frameValue", None),
            getattr(frame, "eigenvalue", None),
            description_value,
        ]
    for candidate in ordered_candidates:
        if candidate is None:
            continue
        try:
            return float(candidate)
        except (TypeError, ValueError):
            continue
    return None


def _extract_buckling_mode_shapes(step, out_dir: str) -> bool:
    _debug_log(out_dir, "DEBUG buckling mode shapes: start")
    frames = _spectral_frames_with_nontrivial_mode_data(step)
    _debug_log(out_dir, f"DEBUG buckling mode shapes: frames={len(frames)}")
    if len(frames) == 0:
        return False
    rows: list[tuple[int, int, str, float]] = []
    for mode_index, frame in enumerate(frames, start=1):
        u_field = _field_output_get(getattr(frame, "fieldOutputs", None), "U")
        if u_field is None:
            continue
        for value in list(getattr(u_field, "values", [])):
            node_label = int(getattr(value, "nodeLabel", 0))
            raw_data = getattr(value, "data", None)
            data = list(raw_data) if raw_data is not None else []
            for dof_index, component in enumerate(data[:3], start=1):
                rows.append((mode_index, node_label, f"U{dof_index}", float(component)))
    if len(rows) == 0:
        return False
    csv_path = os.path.join(out_dir, "buckling_mode_shapes.csv")
    with open(csv_path, "w") as f:
        f.write("mode_index,node_label,dof,value\n")
        for mode_index, node_label, dof, value in rows:
            f.write(f"{mode_index},{node_label},{dof},{value}\n")
    return True


def _maybe_extract_buckling_load_factors(step, out_dir: str, *, simulation_type: str) -> bool:
    if str(simulation_type).strip().lower() != "buckling":
        return False
    return _extract_buckling_load_factors(step, out_dir)


def _maybe_extract_buckling_mode_shapes(step, out_dir: str, *, simulation_type: str) -> bool:
    if str(simulation_type).strip().lower() != "buckling":
        return False
    return _extract_buckling_mode_shapes(step, out_dir)


def _write_missing_spectral_artifact_placeholders(
    out_dir: str,
    *,
    simulation_type: str,
    eigen_frequencies_written: bool = False,
    eigen_mode_shapes_written: bool = False,
    buckling_load_factors_written: bool = False,
    buckling_mode_shapes_written: bool = False,
) -> None:
    sim = str(simulation_type).strip().lower()
    if sim == "eigen":
        if not eigen_frequencies_written:
            with open(os.path.join(out_dir, "eigen_frequencies.csv"), "w") as f:
                f.write("mode_index,frequency_hz\n")
        if not eigen_mode_shapes_written:
            with open(os.path.join(out_dir, "mode_shapes.csv"), "w") as f:
                f.write("mode_index,node_label,dof,value\n")
    elif sim == "buckling":
        if not buckling_load_factors_written:
            with open(os.path.join(out_dir, "buckling_load_factors.csv"), "w") as f:
                f.write("mode_index,load_factor\n")
        if not buckling_mode_shapes_written:
            with open(os.path.join(out_dir, "buckling_mode_shapes.csv"), "w") as f:
                f.write("mode_index,node_label,dof,value\n")


def _extract_eigen_frequencies(step, out_dir: str) -> bool:
    _debug_log(out_dir, "DEBUG eigen frequencies: start")
    frames = _spectral_frames_with_nontrivial_mode_data(step)
    _debug_log(out_dir, f"DEBUG eigen frequencies: frames={len(frames)}")
    if len(frames) == 0:
        return False
    rows: list[tuple[int, float]] = []
    for i, frame in enumerate(frames, start=1):
        frequency = _spectral_frequency_like_value(frame, prefer_description=False)
        if frequency is None:
            continue
        rows.append((i, float(frequency)))
    if len(rows) == 0:
        return False
    csv_path = os.path.join(out_dir, "eigen_frequencies.csv")
    with open(csv_path, "w") as f:
        f.write("mode_index,frequency_hz\n")
        for mode_index, frequency_hz in rows:
            f.write(f"{mode_index},{frequency_hz}\n")
    return True


def _extract_eigen_mode_shapes(step, out_dir: str) -> bool:
    _debug_log(out_dir, "DEBUG eigen mode shapes: start")
    frames = _spectral_frames_with_nontrivial_mode_data(step)
    _debug_log(out_dir, f"DEBUG eigen mode shapes: frames={len(frames)}")
    if len(frames) == 0:
        return False
    rows: list[tuple[int, int, str, float]] = []
    for mode_index, frame in enumerate(frames, start=1):
        u_field = _field_output_get(getattr(frame, "fieldOutputs", None), "U")
        if u_field is None:
            continue
        for value in list(getattr(u_field, "values", [])):
            node_label = int(getattr(value, "nodeLabel", 0))
            raw_data = getattr(value, "data", None)
            data = list(raw_data) if raw_data is not None else []
            for dof_index, component in enumerate(data[:3], start=1):
                rows.append((mode_index, node_label, f"U{dof_index}", float(component)))
    if len(rows) == 0:
        return False
    csv_path = os.path.join(out_dir, "mode_shapes.csv")
    with open(csv_path, "w") as f:
        f.write("mode_index,node_label,dof,value\n")
        for mode_index, node_label, dof, value in rows:
            f.write(f"{mode_index},{node_label},{dof},{value}\n")
    return True


def _maybe_extract_eigen_frequencies(step, out_dir: str, *, simulation_type: str) -> bool:
    if str(simulation_type).strip().lower() != "eigen":
        return False
    return _extract_eigen_frequencies(step, out_dir)


def _maybe_extract_eigen_mode_shapes(step, out_dir: str, *, simulation_type: str) -> bool:
    if str(simulation_type).strip().lower() != "eigen":
        return False
    return _extract_eigen_mode_shapes(step, out_dir)


def _tip_node_label_from_u_field(u_values) -> int | None:
    node_labels = sorted({v.nodeLabel for v in u_values})
    return max(node_labels) if node_labels else None


def _tip_u2_from_field_values(u_values, tip_node: int | None) -> float:
    if tip_node is None:
        return 0.0
    vals = [v for v in u_values if v.nodeLabel == tip_node]
    return float(vals[0].data[1]) if len(vals) > 0 and len(vals[0].data) >= 2 else 0.0


def _write_missing_reference_exports(out_dir: str, *, simulation_type: str) -> None:
    sim = str(simulation_type).strip().lower()
    if sim == "harmonic":
        _write_empty_harmonic_export(out_dir)
    elif sim == "eigen":
        _write_empty_eigen_exports(out_dir)
    elif sim == "buckling":
        _write_empty_buckling_exports(out_dir)
    elif sim == "transient":
        _write_empty_transient_export(out_dir)
    elif sim == "nonlinear_buckling":
        _write_empty_nl_buckling_export(out_dir)


def _extract_transient_tip_displacement_history(step, out_dir: str) -> bool:
    _debug_log(out_dir, "DEBUG transient tip history: start")
    frames = list(getattr(step, "frames", []))
    _debug_log(out_dir, f"DEBUG transient tip history: frames={len(frames)}")
    if len(frames) == 0:
        return False
    tip_node: int | None = None
    rows: list[tuple[float, float]] = []
    for frame in frames:
        u_field = _field_output_get(getattr(frame, "fieldOutputs", None), "U")
        if u_field is None:
            continue
        u_values = list(getattr(u_field, "values", []))
        if tip_node is None:
            tip_node = _tip_node_label_from_u_field(u_values)
        time = float(getattr(frame, "frameValue", 0.0))
        tip_u2 = _tip_u2_from_field_values(u_values, tip_node)
        rows.append((time, tip_u2))
    if len(rows) == 0:
        return False
    csv_path = os.path.join(out_dir, "transient_reference.csv")
    with open(csv_path, "w") as f:
        f.write("time,tip_displacement\n")
        for time, tip_u2 in rows:
            f.write(f"{time:.12e},{tip_u2:.12e}\n")
    return True


def _maybe_extract_transient_tip_displacement_history(step, out_dir: str, *, simulation_type: str) -> bool:
    if str(simulation_type).strip().lower() != "transient":
        return False
    return _extract_transient_tip_displacement_history(step, out_dir)


def _extract_nl_buckling_continuation_history(step, out_dir: str) -> bool:
    """Extract arc-length continuation path from a Riks ODB step.

    In Abaqus Riks, frame.frameValue is the load proportionality factor (LPF).
    Every frame in the ODB corresponds to a converged increment.
    """
    _debug_log(out_dir, "DEBUG nl buckling continuation: start")
    frames = list(getattr(step, "frames", []))
    _debug_log(out_dir, f"DEBUG nl buckling continuation: frames={len(frames)}")
    if len(frames) == 0:
        return False
    tip_node: int | None = None
    rows: list[tuple[float, float]] = []
    for frame in frames:
        u_field = _field_output_get(getattr(frame, "fieldOutputs", None), "U")
        if u_field is None:
            continue
        u_values = list(getattr(u_field, "values", []))
        if tip_node is None:
            tip_node = _tip_node_label_from_u_field(u_values)
        load_factor = float(getattr(frame, "frameValue", 0.0))
        tip_u2 = _tip_u2_from_field_values(u_values, tip_node)
        rows.append((load_factor, tip_u2))
    if len(rows) == 0:
        return False
    csv_path = os.path.join(out_dir, "continuation_history.csv")
    with open(csv_path, "w") as f:
        f.write("load_factor,tip_displacement,converged\n")
        for load_factor, tip_u2 in rows:
            f.write(f"{load_factor:.12e},{tip_u2:.12e},True\n")
    return True


def _maybe_extract_nl_buckling_continuation_history(step, out_dir: str, *, simulation_type: str) -> bool:
    if str(simulation_type).strip().lower() != "nonlinear_buckling":
        return False
    return _extract_nl_buckling_continuation_history(step, out_dir)


def _extract_harmonic_frequency_response(step, out_dir: str) -> bool:
    _debug_log(out_dir, "DEBUG harmonic export: start")
    rows = _harmonic_rows_from_step(step, out_dir=out_dir)
    _debug_log(out_dir, f"DEBUG harmonic export: rows_built={len(rows)}")
    if len(rows) == 0:
        return False
    csv_path = os.path.join(out_dir, "frequency_response.csv")
    with open(csv_path, "w") as f:
        f.write("frequency_hz,response_real,response_imag\n")
        for frequency, response_real, response_imag in rows:
            f.write(f"{float(frequency)},{response_real},{response_imag}\n")
    return True


def _harmonic_frame_frequency(frame, *, frame_index: int) -> float:
    frequency = getattr(frame, "frequency", None)
    if frequency is None:
        frequency = getattr(frame, "frameValue", None)
    if frequency is None:
        description = getattr(frame, "description", None)
        if isinstance(description, str):
            frequency = _buckling_numeric_token_from_description(description)
    if frequency is None:
        frequency = float(frame_index)
    return float(frequency)


def _harmonic_response_component_index_from_value(value) -> int | None:
    raw_data = getattr(value, "data", None)
    data = list(raw_data) if raw_data is not None else []
    if len(data) >= 2:
        return 1
    if len(data) >= 1:
        return 0
    return None


def _harmonic_response_row_from_frame(frame, *, frame_index: int, out_dir: str | None = None) -> tuple[float, float, float]:
    frequency = _harmonic_frame_frequency(frame, frame_index=frame_index)
    response_real, response_imag = _harmonic_response_components_from_frame(frame, out_dir=out_dir)
    return frequency, response_real, response_imag


def _harmonic_rows_from_step(step, *, out_dir: str | None = None) -> list[tuple[float, float, float]]:
    rows: list[tuple[float, float, float]] = []
    for i, frame in enumerate(list(getattr(step, "frames", []))):
        rows.append(_harmonic_response_row_from_frame(frame, frame_index=i, out_dir=out_dir))
    return rows


def _harmonic_has_meaningful_rows(rows: list[tuple[float, float, float]]) -> bool:
    for _, response_real, response_imag in rows:
        if abs(float(response_real)) > 0.0:
            return True
        if abs(float(response_imag)) > 0.0:
            return True
    return False


def _harmonic_response_components_from_frame(frame, *, out_dir: str | None = None) -> tuple[float, float]:
    response_real = 0.0
    response_imag = 0.0
    u_field = _field_output_get(getattr(frame, "fieldOutputs", None), "U")
    if out_dir is not None:
        _debug_log(out_dir, f"DEBUG harmonic components: u_field_is_none={u_field is None}")
    if u_field is None:
        return response_real, response_imag
    tip_value = _harmonic_tip_value_from_u_field(u_field)
    if out_dir is not None:
        _debug_log(out_dir, f"DEBUG harmonic components: tip_value_is_none={tip_value is None}")
    if tip_value is None:
        return response_real, response_imag
    raw_data = getattr(tip_value, "data", None)
    data = list(raw_data) if raw_data is not None else []
    component_index = _harmonic_response_component_index_from_value(tip_value)
    if out_dir is not None:
        _debug_log(out_dir, f"DEBUG harmonic components: data_len={len(data)} component_index={component_index}")
    if component_index is None:
        return response_real, response_imag
    if len(data) > component_index:
        response_real = float(data[component_index])
    response_imag = _harmonic_imag_component_from_value(tip_value, component_index=component_index)
    return response_real, response_imag


def _harmonic_tip_value_from_u_field(u_field):
    values = list(getattr(u_field, "values", []))
    if len(values) == 0:
        return None
    return max(values, key=lambda v: getattr(v, "nodeLabel", 0))


def _harmonic_imag_component_from_value(value, *, component_index: int) -> float:
    conj = getattr(value, "conjugateData", None)
    if conj is None:
        return 0.0
    conj_data = list(conj)
    if len(conj_data) > component_index:
        return float(conj_data[component_index])
    return 0.0


def _maybe_extract_harmonic_frequency_response(step, out_dir: str, *, simulation_type: str) -> bool:
    if str(simulation_type).strip().lower() != "harmonic":
        return False
    return _extract_harmonic_frequency_response(step, out_dir)


def extract_odb_to_csv(
    odb_path: str,
    out_dir: str,
    *,
    export_tip_history: bool = False,
    simulation_type: str = "static",
) -> None:
    """
    Read ODB at odb_path; write U_global.csv (and optionally section_forces.csv) to out_dir.
    Requires Abaqus odbAccess (run from Abaqus Python).
    """
    try:
        import odbAccess
    except ImportError:
        try:
            import odb
            odbAccess = odb
        except ImportError:
            raise RuntimeError("Neither odbAccess nor odb available; run this script from Abaqus Python.")

    if not os.path.exists(odb_path):
        raise FileNotFoundError(f"ODB not found: {odb_path}")

    odb = odbAccess.openOdb(path=odb_path)
    os.makedirs(out_dir, exist_ok=True)

    # Last step, last frame
    step_names = list(odb.steps.keys())
    if not step_names:
        odb.close()
        raise ValueError("No steps in ODB")
    step = odb.steps[step_names[-1]]
    frame = step.frames[-1]
    harmonic_export_written = _maybe_extract_harmonic_frequency_response(
        step,
        out_dir,
        simulation_type=simulation_type,
    )
    eigen_export_written = _maybe_extract_eigen_frequencies(
        step,
        out_dir,
        simulation_type=simulation_type,
    )
    eigen_mode_shapes_written = _maybe_extract_eigen_mode_shapes(
        step,
        out_dir,
        simulation_type=simulation_type,
    )
    buckling_load_factors_written = _maybe_extract_buckling_load_factors(
        step,
        out_dir,
        simulation_type=simulation_type,
    )
    buckling_mode_shapes_written = _maybe_extract_buckling_mode_shapes(
        step,
        out_dir,
        simulation_type=simulation_type,
    )
    transient_written = _maybe_extract_transient_tip_displacement_history(
        step,
        out_dir,
        simulation_type=simulation_type,
    )
    nl_buckling_written = _maybe_extract_nl_buckling_continuation_history(
        step,
        out_dir,
        simulation_type=simulation_type,
    )
    if not any(
        [
            harmonic_export_written,
            eigen_export_written,
            eigen_mode_shapes_written,
            buckling_load_factors_written,
            buckling_mode_shapes_written,
            transient_written,
            nl_buckling_written,
        ]
    ):
        _write_missing_reference_exports(out_dir, simulation_type=simulation_type)
    else:
        _write_missing_spectral_artifact_placeholders(
            out_dir,
            simulation_type=simulation_type,
            eigen_frequencies_written=eigen_export_written,
            eigen_mode_shapes_written=eigen_mode_shapes_written,
            buckling_load_factors_written=buckling_load_factors_written,
            buckling_mode_shapes_written=buckling_mode_shapes_written,
        )
    odb_var_list = list(frame.fieldOutputs.keys())
    has_ur = _field_outputs_has_key(frame.fieldOutputs, "UR")
    print("ODB field outputs:", odb_var_list)
    print("ODB has UR:", has_ur)

    # Displacements U (and UR if available) -> U_global.csv (same convention as FEM: 6 DOFs per node)
    if _field_outputs_has_key(frame.fieldOutputs, "U"):
        u_field = _field_output_get(frame.fieldOutputs, "U")
        node_list = list(u_field.values)
        # U is typically 3 components (U1,U2,U3); UR is separate in Abaqus
        ur_by_node = {}
        if _field_outputs_has_key(frame.fieldOutputs, "UR"):
            for v in _field_output_get(frame.fieldOutputs, "UR").values:
                nlab = v.nodeLabel
                ur_by_node[nlab] = (v.data[0], v.data[1], v.data[2]) if len(v.data) >= 3 else (0.0, 0.0, 0.0)
        csv_path = os.path.join(out_dir, "U_global.csv")
        with open(csv_path, "w") as f:
            f.write("Global DOF,Value\n")
            for v in node_list:
                node_label = v.nodeLabel
                u1, u2, u3 = v.data[0], v.data[1], v.data[2]
                if len(v.data) >= 6:
                    ur1, ur2, ur3 = v.data[3], v.data[4], v.data[5]
                else:
                    ur1, ur2, ur3 = ur_by_node.get(node_label, (0.0, 0.0, 0.0))
                if _field_outputs_has_key(frame.fieldOutputs, "UR"):
                    ur1, ur2, ur3 = -ur1, -ur2, -ur3
                gdof_base = (node_label - 1) * 6
                for d, val in enumerate([u1, u2, u3, ur1, ur2, ur3]):
                    f.write(f"{gdof_base + d},{val}\n")
        rotation_src = "ODB" if _field_outputs_has_key(frame.fieldOutputs, "UR") else "none"
        rotation_src_path = os.path.join(out_dir, "rotation_source.txt")
        with open(rotation_src_path, "w") as rf:
            rf.write(rotation_src + "\n")
        if has_ur:
            print("U_global.csv: U and UR read from ODB and written (rotation sign flipped to match FEM).")
        else:
            print("U_global.csv: U read from ODB; rotation not in ODB, written as zero.")
        print(f"Wrote {csv_path}")

        if export_tip_history:
            tip_hist_path = os.path.join(out_dir, "tip_load_history.csv")
            tip_node = _tip_node_label_from_u_field(node_list)
            frames = list(step.frames)
            with open(tip_hist_path, "w") as tf:
                tf.write("frame_index,load_factor,tip_displacement\n")
                denom = max(1, len(frames) - 1)
                for i, fr in enumerate(frames):
                    if not _field_outputs_has_key(fr.fieldOutputs, "U") or tip_node is None:
                        tf.write(f"{i},{float(i)/denom:.12e},0.0\n")
                        continue
                    tip_u = _tip_u2_from_field_values(_field_output_get(fr.fieldOutputs, "U").values, tip_node)
                    tf.write(f"{i},{float(i)/denom:.12e},{tip_u:.12e}\n")
            print(f"Wrote {tip_hist_path}")

    # Section forces: write section_forces.csv (x, N, Vy, Vz, T, My, Mz) and nodal_section_forces.csv (FEM format)
    if _field_outputs_has_key(frame.fieldOutputs, "SF"):
        sf_field = _field_output_get(frame.fieldOutputs, "SF")
        sm_field = _field_output_get(frame.fieldOutputs, "SM")
        try:
            inst = odb.rootAssembly.instances[list(odb.rootAssembly.instances.keys())[0]]
            node_coords = {n.label: n.coordinates for n in inst.nodes}
            node_labels_sorted = sorted(node_coords.keys())
            elem_conn = {}
            for c in inst.elementConnectivity:
                elem_conn[c[0]] = (c[1], c[2])
        except Exception:
            node_coords = {}
            node_labels_sorted = []
            elem_conn = {}
        rows = []
        elem_sf = {}  # elabel -> (n, vy, vz, t, my, mz)
        sf_has_six = False
        used_sm = False
        for val in sf_field.values:
            elabel = val.elementLabel
            data = list(val.data)
            if len(data) >= 6:
                n, vy, vz, t, my, mz = data[0], data[1], data[2], data[3], data[4], data[5]
                sf_has_six = True
            elif len(data) >= 3 and sm_field is not None:
                n, vy, vz = data[0], data[1], data[2]
                sm_vals = [v for v in sm_field.values if v.elementLabel == elabel]
                if len(sm_vals) > 0 and len(sm_vals[0].data) >= 3:
                    t, my, mz = sm_vals[0].data[0], sm_vals[0].data[1], sm_vals[0].data[2]
                    used_sm = True
                else:
                    t, my, mz = 0.0, 0.0, 0.0
            else:
                n = data[0] if len(data) > 0 else 0.0
                vy = data[1] if len(data) > 1 else 0.0
                vz = data[2] if len(data) > 2 else 0.0
                t, my, mz = 0.0, 0.0, 0.0
            elem_sf[elabel] = (n, vy, vz, t, my, mz)
            try:
                conn = elem_conn.get(elabel)
                if conn is not None and len(node_coords) > 0:
                    n1, n2 = conn[0], conn[1]
                    x_center = (node_coords[n1][0] + node_coords[n2][0]) / 2.0
                else:
                    x_center = float(elabel)
            except Exception:
                x_center = float(elabel)
            rows.append((x_center, n, vy, vz, t, my, mz))
        if len(rows) > 0:
            csv_path = os.path.join(out_dir, "section_forces.csv")
            with open(csv_path, "w") as f:
                f.write("x,N,Vy,Vz,T,My,Mz\n")
                for r in rows:
                    f.write(f"{r[0]},{r[1]},{r[2]},{r[3]},{r[4]},{r[5]},{r[6]}\n")
            if used_sm:
                print("Section forces: SF + SM written to section_forces.csv, nodal_section_forces.csv.")
            elif sf_has_six:
                print("Section forces: SF (6 components) written to section_forces.csv, nodal_section_forces.csv.")
            else:
                print("Section forces: SF only (T,My,Mz from SF or zero); written to section_forces.csv, nodal_section_forces.csv.")
            print(f"Wrote {csv_path}")
        # nodal_section_forces.csv: same format as FEM (one row per node, N,Vy,Vz,T,My,Mz; node order 1..n)
        if len(node_labels_sorted) > 0 and len(elem_conn) > 0 and len(elem_sf) > 0:
            nodal_sum = defaultdict(lambda: [0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
            nodal_count = defaultdict(int)
            for elabel, (n, vy, vz, t, my, mz) in elem_sf.items():
                conn = elem_conn.get(elabel)
                if conn is None:
                    continue
                for nid in conn:
                    s = nodal_sum[nid]
                    s[0] += n
                    s[1] += vy
                    s[2] += vz
                    s[3] += t
                    s[4] += my
                    s[5] += mz
                    nodal_count[nid] += 1
            nodal_path = os.path.join(out_dir, "nodal_section_forces.csv")
            with open(nodal_path, "w") as f:
                f.write("# column_order=resultant\n")
                f.write("N,Vy,Vz,T,My,Mz\n")
                for nid in node_labels_sorted:
                    c = nodal_count[nid]
                    s = nodal_sum[nid]
                    if c > 0:
                        f.write("%.12e,%.12e,%.12e,%.12e,%.12e,%.12e\n" % (
                            s[0] / c, s[1] / c, s[2] / c, s[3] / c, s[4] / c, s[5] / c))
                    else:
                        f.write("0.000000000000e+00,0.000000000000e+00,0.000000000000e+00,0.000000000000e+00,0.000000000000e+00,0.000000000000e+00\n")
            print(f"Wrote {nodal_path}")

    odb.close()


def main() -> None:
    # Parse -- key=value args after --
    args = {}
    for arg in sys.argv[1:]:
        if "=" in arg:
            k, v = arg.split("=", 1)
            args[k.strip()] = v.strip()
    odb_path = args.get("odb_path")
    out_dir = args.get("out_dir")
    if not odb_path or not out_dir:
        print("Usage: abaqus cae noGUI=extract_odb_results.py -- odb_path=<path> out_dir=<path>", file=sys.stderr)
        sys.exit(1)
    export_tip_history = str(args.get("export_tip_history", "false")).strip().lower() in ("1", "true", "yes", "on")
    simulation_type = str(args.get("simulation_type", "static")).strip().lower()
    extract_odb_to_csv(odb_path, out_dir, export_tip_history=export_tip_history, simulation_type=simulation_type)


if __name__ == "__main__":
    main()

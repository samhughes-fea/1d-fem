# post_processing/validation_visualisers/abaqus/job_to_abaqus_script.py
"""
Generate an Abaqus CAE Python script from a job directory.
The script creates model, wire part, beam section, material, BCs, loads,
runs the job, and exports ODB results to CSV.
Supports transient validation script generation contract metadata for
external-reference workflows.
Run with project Python (abqpy): python run_<job>.py; abqpy's saveAs() launches Abaqus.
"""
from __future__ import annotations

import argparse
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

# Project root for importing parsers
SCRIPT_DIR = Path(__file__).resolve().parent
VALIDATION_DIR = SCRIPT_DIR.parent
PROJECT_ROOT = VALIDATION_DIR.parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from pre_processing.parsing.grid_parser import GridParser
from pre_processing.parsing.element_parser import ElementParser
from pre_processing.parsing.material_parser import MaterialParser
from pre_processing.parsing.section_parser import SectionParser
from pre_processing.parsing.prescribed_displacement_parser import parse_prescribed_displacement
from pre_processing.element_library.beam_warping import mesh_uses_warping_dof
from pre_processing.parsing.point_load_parser import parse_point_load

from post_processing.validation_visualisers.abaqus.config import (
    JOBS_DIR,
    ABAQUS_GENERATED_DIR,
    ABAQUS_RESULTS_DIR,
    ELEMENT_TYPE_MAP,
    SUPPORTED_ELEMENT_TYPES,
)
from post_processing.validation_visualisers.abaqus.simulation_type_dispatch import build_validation_dispatch_payload
from pre_processing.element_library.linear.beam.zero_order_shear_deformation_theory.euler_bernoulli.utilities.interpolate_loads import (
    LoadInterpolationOperator,
)
from pre_processing.element_library.linear.beam.zero_order_shear_deformation_theory.euler_bernoulli.utilities.shape_functions import (
    ShapeFunctionOperator,
)


def _equivalent_nodal_from_distributed_load(
    coords: list,
    connectivity: list,
    distributed_loads: list,
) -> list:
    """
    Compute equivalent nodal forces (and moments) from distributed load q(x) so that
    Abaqus matches the FEM. The job file distributed_load.txt stores q(x) in N/m
    (and moment per length) at sample points; the FEM integrates N^T q over each
    element. We replicate that integration here and output (node_label, Fx,Fy,Fz,
    Mx,My,Mz) per node (node_label 1-based).
    """
    if not distributed_loads or not connectivity or not coords:
        return []
    coords = np.asarray(coords, dtype=np.float64)
    dl_array = np.array(
        [tuple(float(r[k]) for k in range(9)) for r in distributed_loads if len(r) >= 9],
        dtype=np.float64,
    )
    if dl_array.size == 0:
        return []
    try:
        interpolator = LoadInterpolationOperator(
            distributed_loads_array=dl_array,
            boundary_mode="error",
            interpolation_order="cubic",
            n_gauss_points=3,
        )
    except Exception:
        return []
    xi_gauss, weights = np.polynomial.legendre.leggauss(3)
    nodal = defaultdict(lambda: np.zeros(6, dtype=np.float64))
    for elem in connectivity:
        if len(elem) != 2:
            continue
        na, nb = int(elem[0]), int(elem[1])
        xa, xb = coords[na][0], coords[nb][0]
        x_start = min(xa, xb)
        L = max(xa, xb) - x_start
        if L <= 0:
            continue
        shape_op = ShapeFunctionOperator(element_length=L)
        x_gauss = (xi_gauss + 1) * (L / 2) + x_start
        try:
            q_gauss = interpolator.interpolate(x_gauss)
        except ValueError:
            continue
        N = np.stack(
            [
                shape_op.natural_coordinate_form(xi)[0][0]
                for xi in xi_gauss
            ]
        )
        Fe = np.einsum("gij,gj,g->i", N, q_gauss, weights) * (L / 2)
        nodal[na][:] += Fe[:6]
        nodal[nb][:] += Fe[6:12]
    return [
        (int(n) + 1, float(v[0]), float(v[1]), float(v[2]), float(v[3]), float(v[4]), float(v[5]))
        for n, v in sorted(nodal.items())
    ]


def _parse_job(job_dir: Path) -> dict:
    """Parse job directory and return a single dict with all data for script generation."""
    job_dir = Path(job_dir)
    job_str = str(job_dir)
    grid = GridParser(str(job_dir / "grid.txt"), job_str).parse()
    element = ElementParser(str(job_dir / "element.txt"), job_str).parse()
    material = MaterialParser(str(job_dir / "material.txt"), job_str).parse()
    section = SectionParser(str(job_dir / "section.txt"), job_str).parse()

    gd = grid["grid_dictionary"]
    ed = element["element_dictionary"]
    md = material["material_dictionary"]
    sd = section["section_dictionary"]

    coords = gd["coordinates"].tolist()
    node_ids = gd["ids"].tolist()
    connectivity = ed["connectivity"].tolist()
    types_arr = ed["types"]
    # Single element type for whole model (all elements same in our jobs)
    elem_type = str(types_arr[0]) if len(types_arr) > 0 else "EulerBernoulliBeamElement3D"
    if elem_type not in SUPPORTED_ELEMENT_TYPES:
        raise ValueError(
            f"Unsupported element type for Abaqus validation: {elem_type}. "
            f"Supported: {sorted(SUPPORTED_ELEMENT_TYPES)}"
        )
    abaqus_elem = ELEMENT_TYPE_MAP[elem_type]

    # Material: use first element's E, G, nu, rho
    E = float(md["E"][0])
    G = float(md["G"][0])
    nu = float(md["nu"][0])
    rho = float(md["rho"][0]) if "rho" in md else 0.0

    # Section: A, I_y -> I11, I_z -> I22, J_t -> J (first element)
    A = float(sd["A"][0])
    I_y = float(sd["I_y"][0])
    I_z = float(sd["I_z"][0])
    J_t = float(sd["J_t"][0])

    # Prescribed displacement (7 DOF/node when warping mesh is on)
    pres_path = job_dir / "prescribed_displacement.txt"
    if pres_path.is_file():
        _dpn = 7 if mesh_uses_warping_dof(ed) else 6
        pres = parse_prescribed_displacement(str(pres_path), dof_per_node=_dpn)
    else:
        pres = {"node_id": [], "dof": [], "value": [], "dof_index": []}

    # Point loads
    point_path = job_dir / "point_load.txt"
    if point_path.is_file():
        try:
            pl_array = parse_point_load(str(point_path))
            point_loads = pl_array.tolist() if pl_array is not None and pl_array.size > 0 else []
        except Exception:
            point_loads = []
    else:
        point_loads = []

    # Distributed load (optional; first phase we may approximate as equivalent nodal or skip)
    dist_path = job_dir / "distributed_load.txt"
    distributed_loads = []
    if dist_path.is_file():
        try:
            from pre_processing.parsing.distributed_load_parser import parse_distributed_load
            dl_array = parse_distributed_load(str(dist_path))
            if dl_array is not None and dl_array.size > 0:
                distributed_loads = dl_array.tolist()
        except Exception:
            pass

    dispatch = build_validation_dispatch_payload(str(job_dir))

    return {
        "job_name": job_dir.name,
        "coords": coords,
        "node_ids": node_ids,
        "connectivity": connectivity,
        "abaqus_elem": abaqus_elem,
        "E": E,
        "G": G,
        "nu": nu,
        "rho": rho,
        "A": A,
        "I_y": I_y,
        "I_z": I_z,
        "J_t": J_t,
        "prescribed": pres,
        "point_loads": point_loads,
        "distributed_loads": distributed_loads,
        "simulation_settings_path": str(job_dir / "simulation_settings.txt"),
        "simulation_type": dispatch["simulation_type"],
        "simulation_settings": dispatch["simulation_settings"],
        "artifact_contract": dispatch["artifact_contract"],
    }


def _generate_script_content(data: dict, out_csv_dir: str) -> str:
    """Return the full Abaqus Python script as a string."""
    job_name = data["job_name"]
    coords = data["coords"]
    abaqus_elem = data["abaqus_elem"]
    E = data["E"]
    nu = data["nu"]
    rho = data["rho"]
    A = data["A"]
    I_y = data["I_y"]
    I_z = data["I_z"]
    J_t = data["J_t"]
    prescribed = data["prescribed"]
    point_loads = data["point_loads"]
    distributed_loads = data.get("distributed_loads") or []
    simulation_type = str(data.get("simulation_type", "static")).strip().lower()
    simulation_settings = data.get("simulation_settings") or {}
    artifact_contract = data.get("artifact_contract") or {}
    artifact_contract_name = str(artifact_contract.get("contract_name", "static_reference"))
    expected_files = list(artifact_contract.get("expected_files") or [])
    nonlinear_cfg = simulation_settings.get("nonlinear") or {}
    nonlinear_num_increments = int(nonlinear_cfg.get("num_increments", 1))

    # Build prescribed list as literal for generated script (Abaqus has no access to 'prescribed' dict)
    node_ids = prescribed.get("node_id", [])
    dof_indices = prescribed.get("dof_index", [])
    values = prescribed.get("value", [])
    prescribed_tuples = [(int(n), int(d), float(v)) for n, d, v in zip(node_ids, dof_indices, values)]

    # Distributed loads: use UDL (LineLoad) when Fy is constant; else equivalent nodal forces (triangular/parabolic).
    # For non-UDL, distributed_load.txt stores q(x) in N/m; we integrate N^T q over each element to match the FEM.
    distributed_equivalent_nodal = []
    connectivity = data.get("connectivity") or []
    if distributed_loads and len(coords) > 0:
        fy_vals = [row[4] for row in distributed_loads if len(row) >= 5]
        if fy_vals:
            tol = 1e-6 * (max(abs(v) for v in fy_vals) + 1e-10)
            use_udl = (max(fy_vals) - min(fy_vals)) <= tol
            if not use_udl:
                distributed_equivalent_nodal = _equivalent_nodal_from_distributed_load(
                    coords, connectivity, distributed_loads
                )

    # Abaqus script runs inside Abaqus CAE; no project imports available
    # Embed extraction logic and output path
    out_csv_dir_escaped = out_csv_dir.replace("\\", "\\\\")

    step_block = (
        'model.ImplicitDynamicsStep(name="Step-1", previous="Initial", timePeriod=0.05, maxNumInc=1000, initialInc=0.01, minInc=1e-06, nlgeom=OFF)'
        if simulation_type == "transient"
        else ('model.StaticStep(name="Step-1", previous="Initial", description="Nonlinear static")' if simulation_type == "static" and nonlinear_num_increments > 1 else 'model.StaticStep(name="Step-1", previous="Initial", description="Static")')
    )
    script = (
        _build_script_preamble(
        job_name=job_name,
        out_csv_dir=out_csv_dir,
        coords=coords,
        abaqus_elem=abaqus_elem,
        E=E,
        nu=nu,
        rho=rho,
        A=A,
        I_y=I_y,
        I_z=I_z,
        J_t=J_t,
        prescribed_tuples=prescribed_tuples,
        point_loads=point_loads,
        distributed_loads=distributed_loads,
        distributed_equivalent_nodal=distributed_equivalent_nodal,
        out_csv_dir_escaped=out_csv_dir_escaped,
        artifact_contract_name=artifact_contract_name,
        expected_files=expected_files,
        )
        + _build_step_and_model_block(step_block)
        + _build_loads_and_job_block(job_name)
        + _build_results_export_prologue_block()
        + _build_results_export_displacement_block()
        + f'''
# Request U, UR (deflection/rotation) and SF, SM (section forces and section moments).
# UR is nodal rotation; SF = N,Vy,Vz; SM = T,My,Mz. If model.FieldOutputRequest is missing (e.g. Abaqus 2021), try step-based API.
FIELD_OUTPUT_USED = "none"
try:
    model.FieldOutputRequest(name="F-Output-1", createStepName="Step-1", variables=("U", "UR", "SF", "SM"))
    FIELD_OUTPUT_USED = "U, UR, SF, SM"
    _log("Field output requested: U, UR, SF, SM")
except Exception as e_fld:
    _log("Field output request (U, UR, SF, SM) failed: " + str(e_fld))
    step = model.steps["Step-1"]
    try:
        step.FieldOutputRequest(name="F-Output-1", variables=("U", "UR", "SF", "SM"))
        FIELD_OUTPUT_USED = "U, UR, SF, SM"
        _log("Field output requested: U, UR, SF, SM (via step)")
    except Exception:
        try:
            step.FieldOutputRequest(name="F-Output-1", variables=("U", "UR"))
            FIELD_OUTPUT_USED = "U, UR"
            _log("Field output requested: U, UR (SF/SM omitted, via step)")
        except Exception as e_ur:
            _log("Field output request (U, UR) failed: " + str(e_ur))
            try:
                step.FieldOutputRequest(name="F-Output-1", variables=("SF", "SM"))
                FIELD_OUTPUT_USED = "SF, SM only"
                _log("Field output requested: SF, SM only (via step)")
            except Exception:
                try:
                    step.FieldOutputRequest(name="F-Output-1", variables=("SF",))
                    FIELD_OUTPUT_USED = "SF only"
                    _log("Field output requested: SF only (via step)")
                except Exception as e2:
                    _log("Field output request failed (no U, UR, SF, or SM): " + str(e2))
    _log("(Abaqus may still write U/UR by default; see ODB field outputs after the job.)")

# --- BCs: encastre at first node (typical cantilever) ---
try:
    n_set = assembly.instances["Beam-1"].nodes.sequenceFromLabels(labels=(1,))
    region = assembly.Set(nodes=n_set, name="Fixed")
    model.DisplacementBC(name="Encastre", createStepName="Initial", region=region, u1=SET, u2=SET, u3=SET, ur1=SET, ur2=SET, ur3=SET)
    model.boundaryConditions["Encastre"].setValuesInStep(stepName="Step-1", u1=0, u2=0, u3=0, ur1=0, ur2=0, ur3=0)
except Exception as e:
    print("BC warning: " + str(e))

# --- Distributed load: equivalent nodal (triangular/parabolic) or LineLoad (UDL) ---
# Job file inputs represented identically: apply only non-zero (Fx,Fy,Fz) and (Mx,My,Mz); Abaqus forbids zero-magnitude ConcentratedForce
if DISTRIBUTED_EQUIVALENT_NODAL:
    try:
        _tol = 1e-12
        for i, (node_label, fx, fy, fz, mx, my, mz) in enumerate(DISTRIBUTED_EQUIVALENT_NODAL):
            n_set = assembly.instances["Beam-1"].nodes.sequenceFromLabels(labels=(node_label,))
            region = assembly.Set(nodes=n_set, name="DistLoadNode_" + str(i))
            if abs(fx) > _tol or abs(fy) > _tol or abs(fz) > _tol:
                model.ConcentratedForce(name="DistCF_" + str(i), createStepName="Step-1", region=region, cf1=fx, cf2=fy, cf3=fz)
            if abs(mx) > _tol or abs(my) > _tol or abs(mz) > _tol:
                model.Moment(name="DistMoment_" + str(i), createStepName="Step-1", region=region, cm1=mx, cm2=my, cm3=mz)
    except Exception as e:
        print("Distributed equivalent nodal load warning: " + str(e))
elif DISTRIBUTED_LOADS:
    try:
        fy_vals = [row[4] for row in DISTRIBUTED_LOADS if len(row) >= 5]
        if fy_vals:
            comp2_mag = sum(fy_vals) / len(fy_vals)
            beam_edges = assembly.instances["Beam-1"].edges
            region = assembly.Set(edges=beam_edges, name="BeamEdges")
            model.LineLoad(name="LineLoad-1", createStepName="Step-1", region=region, comp1=0.0, comp2=comp2_mag, comp3=0.0, distributionType=UNIFORM)
    except Exception as e:
        print("Distributed load warning: " + str(e))

# --- Point loads (apply at nearest node to given x,y,z) ---
for i, row in enumerate(POINT_LOADS):
    if len(row) < 9:
        continue
    x, y, z, Fx, Fy, Fz, Mx, My, Mz = row[0:9]
    idx = min(range(len(COORDS)), key=lambda i: (COORDS[i][0]-x)**2 + (COORDS[i][1]-y)**2 + (COORDS[i][2]-z)**2)
    node_label = idx + 1
    try:
        n_set = assembly.instances["Beam-1"].nodes.sequenceFromLabels(labels=(node_label,))
        region = assembly.Set(nodes=n_set, name="LoadNode_" + str(i))
        model.ConcentratedForce(name="CF_" + str(i), createStepName="Step-1", region=region, cf1=Fx, cf2=Fy, cf3=Fz)
        if Mx != 0 or My != 0 or Mz != 0:
            model.Moment(name="Moment_" + str(i), createStepName="Step-1", region=region, cm1=Mx, cm2=My, cm3=Mz)
    except Exception as e:
        print("Load warning: " + str(e))

# --- Job ---
jobName = "{job_name.replace('-', '_')}_abaqus"
mdb.Job(name=jobName, model=modelName, description="", type=ANALYSIS)
mdb.jobs[jobName].submit(consistencyChecking=OFF)
mdb.jobs[jobName].waitForCompletion()

# --- Export ODB to CSV (displacements) ---
try:
    odb_path = mdb.jobs[jobName].path + ".odb"
except AttributeError:
    odb_path = os.path.join(os.getcwd(), jobName + ".odb")
if os.path.exists(odb_path):
    import odbAccess
    odb = odbAccess.openOdb(path=odb_path)
    step = odb.steps["Step-1"]
    frame = step.frames[-1]
    odb_var_list = list(frame.fieldOutputs.keys())
    _log("ODB field outputs: " + str(odb_var_list))
    _log("ODB has UR: " + str("UR" in frame.fieldOutputs))
    field = frame.fieldOutputs["U"]
    node_list = field.values
    ur_by_node = {{}}
    if "UR" in frame.fieldOutputs:
        for _v in frame.fieldOutputs["UR"].values:
            nlab = _v.nodeLabel
            ur_by_node[nlab] = (_v.data[0], _v.data[1], _v.data[2]) if len(_v.data) >= 3 else (0.0, 0.0, 0.0)
    if not os.path.exists(OUT_CSV_DIR):
        os.makedirs(OUT_CSV_DIR)
    csv_path = os.path.join(OUT_CSV_DIR, "U_global.csv")
    with open(csv_path, "w") as f:
        f.write("Global DOF,Value\\n")
        for v in node_list:
            nodeLabel = v.nodeLabel
            u1, u2, u3 = v.data[0], v.data[1], v.data[2]
            if len(v.data) >= 6:
                ur1, ur2, ur3 = v.data[3], v.data[4], v.data[5]
            else:
                ur1, ur2, ur3 = ur_by_node.get(nodeLabel, (0.0, 0.0, 0.0))
            if "UR" in frame.fieldOutputs:
                ur1, ur2, ur3 = -ur1, -ur2, -ur3
            gdof_base = (nodeLabel - 1) * 6
            for d, val in enumerate([u1, u2, u3, ur1, ur2, ur3]):
                f.write(str(gdof_base + d) + "," + str(val) + "\\n")
    if "{simulation_type}" == "transient":
        with open(os.path.join(OUT_CSV_DIR, "transient_reference_contract.txt"), "w") as _tf:
            _tf.write("simulation_type=transient\\n")
            _tf.write("step_name=Step-1\\n")
            _tf.write("expected_reference=time_history_csv\\n")
    if "{simulation_type}" == "static" and {nonlinear_num_increments} > 1:
        try:
            from post_processing.validation_visualisers.abaqus.extract_odb_results import extract_odb_to_csv as _extract_tip
            _extract_tip(odb_path, OUT_CSV_DIR, export_tip_history=True)
        except Exception as _e_tip:
            _log("tip_load_history export failed: " + str(_e_tip))
    rotation_src = "ODB" if "UR" in frame.fieldOutputs else "none"
    with open(os.path.join(OUT_CSV_DIR, "rotation_source.txt"), "w") as _rf:
        _rf.write(rotation_src + "\\n")
    with open(os.path.join(OUT_CSV_DIR, "artifact_contract.txt"), "w") as _cf:
        _cf.write("contract_name=" + ARTIFACT_CONTRACT_NAME + "\\n")
        _cf.write("expected_files=" + ",".join(EXPECTED_ARTIFACT_FILES) + "\\n")
    # --- Export section forces (x, N, Vy, Vz, T, My, Mz) for comparison ---
    if "SF" in frame.fieldOutputs:
        sf_field = frame.fieldOutputs["SF"]
        sm_field = frame.fieldOutputs["SM"] if "SM" in frame.fieldOutputs else None
        sf_path = os.path.join(OUT_CSV_DIR, "section_forces.csv")
        with open(sf_path, "w") as sf_file:
            sf_file.write("x,N,Vy,Vz,T,My,Mz\\n")
            for val in sf_field.values:
                elabel = val.elementLabel
                data = list(val.data)
                if len(data) >= 6:
                    n, vy, vz, t, my, mz = data[0], data[1], data[2], data[3], data[4], data[5]
                elif len(data) >= 3 and sm_field is not None:
                    n, vy, vz = data[0], data[1], data[2]
                    sm_vals = [v for v in sm_field.values if v.elementLabel == elabel]
                    t, my, mz = (sm_vals[0].data[0], sm_vals[0].data[1], sm_vals[0].data[2]) if sm_vals and len(sm_vals[0].data) >= 3 else (0.0, 0.0, 0.0)
                else:
                    n, vy, vz, t, my, mz = data[0] if len(data) > 0 else 0, data[1] if len(data) > 1 else 0, data[2] if len(data) > 2 else 0, 0.0, 0.0, 0.0
                i1, i2 = elabel - 1, min(elabel, len(COORDS) - 1)
                x_center = (float(COORDS[i1][0]) + float(COORDS[i2][0])) / 2.0
                sf_file.write(str(x_center) + "," + str(n) + "," + str(vy) + "," + str(vz) + "," + str(t) + "," + str(my) + "," + str(mz) + "\\n")
    odb.close()
    # Copy .odb, .inp, .sta, .msg into result dir so each job dir is self-contained
    try:
        base = odb_path.replace(".odb", "")
        for ext in (".odb", ".inp", ".sta", ".msg"):
            src = base + ext
            if os.path.exists(src):
                shutil.copy2(src, os.path.join(OUT_CSV_DIR, os.path.basename(src)))
    except Exception as e_copy:
        _log("Copy inp/odb/sta/msg to result dir failed: " + str(e_copy))

    # --- End ---
'''
    )
    return script


def _build_results_export_prologue_block() -> str:
    return '''
# --- Export ODB to CSV (displacements) ---
try:
    odb_path = mdb.jobs[jobName].path + ".odb"
except AttributeError:
    odb_path = os.path.join(os.getcwd(), jobName + ".odb")
if os.path.exists(odb_path):
    import odbAccess
    odb = odbAccess.openOdb(path=odb_path)
    step = odb.steps["Step-1"]
    frame = step.frames[-1]
    odb_var_list = list(frame.fieldOutputs.keys())
    _log("ODB field outputs: " + str(odb_var_list))
    _log("ODB has UR: " + str("UR" in frame.fieldOutputs))
    field = frame.fieldOutputs["U"]
    node_list = field.values
    ur_by_node = {}
    if "UR" in frame.fieldOutputs:
        for _v in frame.fieldOutputs["UR"].values:
            nlab = _v.nodeLabel
            ur_by_node[nlab] = (_v.data[0], _v.data[1], _v.data[2]) if len(_v.data) >= 3 else (0.0, 0.0, 0.0)
'''


def _build_results_export_displacement_block() -> str:
    return '''
    if not os.path.exists(OUT_CSV_DIR):
        os.makedirs(OUT_CSV_DIR)
    csv_path = os.path.join(OUT_CSV_DIR, "U_global.csv")
    with open(csv_path, "w") as f:
        f.write("Global DOF,Value\\n")
        for v in node_list:
            nodeLabel = v.nodeLabel
            u1, u2, u3 = v.data[0], v.data[1], v.data[2]
            if len(v.data) >= 6:
                ur1, ur2, ur3 = v.data[3], v.data[4], v.data[5]
            else:
                ur1, ur2, ur3 = ur_by_node.get(nodeLabel, (0.0, 0.0, 0.0))
            if "UR" in frame.fieldOutputs:
                ur1, ur2, ur3 = -ur1, -ur2, -ur3
            gdof_base = (nodeLabel - 1) * 6
            for d, val in enumerate([u1, u2, u3, ur1, ur2, ur3]):
                f.write(str(gdof_base + d) + "," + str(val) + "\\n")
'''


def _build_script_preamble(
    *,
    job_name: str,
    out_csv_dir: str,
    coords: list,
    abaqus_elem: str,
    E: float,
    nu: float,
    rho: float,
    A: float,
    I_y: float,
    I_z: float,
    J_t: float,
    prescribed_tuples: list,
    point_loads: list,
    distributed_loads: list,
    distributed_equivalent_nodal: list,
    out_csv_dir_escaped: str,
    artifact_contract_name: str,
    expected_files: list,
) -> str:
    return f'''# -*- coding: utf-8 -*-
from abaqus import *
from abaqusConstants import *
from driverUtils import *
import regionToolset
import mesh
import os
import sys
import shutil

executeOnCaeStartup()

# Detect if running inside Abaqus (real kernel) vs project Python (abqpy stubs)
try:
    import odbAccess
    IN_ABAQUS = True
except Exception:
    IN_ABAQUS = False

# --- Data (from job dir) ---
COORDS = {coords}
ABAQUS_ELEM = "{abaqus_elem}"
E_MODULUS = {E}
NU = {nu}
RHO = {rho}
A_AREA = {A}
I11 = {I_y}
I22 = {I_z}
J_TORSION = {J_t}
PRESCRIBED_NODES_DOF_VALUES = {prescribed_tuples}
POINT_LOADS = {point_loads}
DISTRIBUTED_LOADS = {distributed_loads}
DISTRIBUTED_EQUIVALENT_NODAL = {distributed_equivalent_nodal}
OUT_CSV_DIR = r"{out_csv_dir_escaped}"
ARTIFACT_CONTRACT_NAME = "{artifact_contract_name}"
EXPECTED_ARTIFACT_FILES = {expected_files}

# Default encastre at node 0 if no prescribed DOFs given
if not PRESCRIBED_NODES_DOF_VALUES and COORDS:
    for dof_idx in range(6):
        PRESCRIBED_NODES_DOF_VALUES.append((0, dof_idx, 0.0))

# When run by project Python (abqpy), only trigger Abaqus launch; model/job/export run only inside Abaqus
if not IN_ABAQUS:
    mdb.saveAs(os.path.join(OUT_CSV_DIR, "model.cae"))
    sys.exit(0)

'''


def _build_step_and_model_block(step_block: str) -> str:
    return f'''
# --- Model ---
modelName = "Model-1"
try:
    del mdb.models[modelName]
except KeyError:
    pass
mdb.Model(name=modelName)
model = mdb.models[modelName]

# --- Part: 3D wire from polyline ---
partName = "Beam"
part = model.Part(name=partName, dimensionality=THREE_D, type=DEFORMABLE_BODY)
points = tuple(tuple(float(c) for c in pt) for pt in COORDS)
part.WirePolyLine(mergeType=IMPRINT, meshable=ON, points=points)

# --- Material ---
matName = "Steel"
model.Material(name=matName)
model.materials[matName].Elastic(table=((E_MODULUS, NU),))
if RHO > 0:
    model.materials[matName].Density(table=((RHO,),))

# --- Beam section: use Rectangular profile (dimensions from A, I22) so Abaqus accepts it ---
# Match cross-section A and I22 (bending); I11 and J will differ for non-rectangular sections.
import math
a_rect = A_AREA * math.sqrt(A_AREA / (12.0 * I22)) if I22 > 0 else math.sqrt(A_AREA)
b_rect = math.sqrt(12.0 * I22 / A_AREA) if A_AREA > 0 and I22 > 0 else math.sqrt(A_AREA)
profileName = "BeamRectProfile"
model.RectangularProfile(name=profileName, a=a_rect, b=b_rect)
sectionName = "BeamSection"
model.BeamSection(
    name=sectionName,
    integration=DURING_ANALYSIS,
    poissonRatio=NU,
    profile=profileName,
    material=matName,
    temperatureVar=LINEAR,
)
region = part.Set(edges=part.edges, name="BeamSet")
part.SectionAssignment(region=region, sectionName=sectionName, offset=0.0, offsetType=MIDDLE_SURFACE)
try:
    part.beamOrientations.assign(region=region, method=N1_COSINES, n1=(0.0, 1.0, 0.0))
except AttributeError:
    pass  # Some Abaqus versions use different API; default orientation may apply

# --- Element type and mesh (seed, setElementType, then generateMesh) ---
part.seedPart(size=1e10, minSizeFactor=0.1, deviationFactor=0.1)
part.setElementType(regions=(part.edges,), elemTypes=(mesh.ElemType(elemCode=ABAQUS_ELEM, elemLibrary=STANDARD),))
part.generateMesh()

# --- Assembly ---
assembly = model.rootAssembly
assembly.DatumCsysByDefault(CARTESIAN)
assembly.Instance(name="Beam-1", part=part, dependent=ON)

# --- Step ---
{step_block}
if not os.path.exists(OUT_CSV_DIR):
    os.makedirs(OUT_CSV_DIR)
def _log(msg):
    print(msg)
    try:
        with open(os.path.join(OUT_CSV_DIR, "run_log.txt"), "a") as _f:
            _f.write(str(msg) + "\\n")
    except Exception:
        pass

'''


def _build_loads_and_job_block(job_name: str) -> str:
    return f'''
# --- BCs: encastre at first node (typical cantilever) ---
try:
    n_set = assembly.instances["Beam-1"].nodes.sequenceFromLabels(labels=(1,))
    region = assembly.Set(nodes=n_set, name="Fixed")
    model.DisplacementBC(name="Encastre", createStepName="Initial", region=region, u1=SET, u2=SET, u3=SET, ur1=SET, ur2=SET, ur3=SET)
    model.boundaryConditions["Encastre"].setValuesInStep(stepName="Step-1", u1=0, u2=0, u3=0, ur1=0, ur2=0, ur3=0)
except Exception as e:
    print("BC warning: " + str(e))

# --- Distributed load: equivalent nodal (triangular/parabolic) or LineLoad (UDL) ---
# Job file inputs represented identically: apply only non-zero (Fx,Fy,Fz) and (Mx,My,Mz); Abaqus forbids zero-magnitude ConcentratedForce
if DISTRIBUTED_EQUIVALENT_NODAL:
    try:
        _tol = 1e-12
        for i, (node_label, fx, fy, fz, mx, my, mz) in enumerate(DISTRIBUTED_EQUIVALENT_NODAL):
            n_set = assembly.instances["Beam-1"].nodes.sequenceFromLabels(labels=(node_label,))
            region = assembly.Set(nodes=n_set, name="DistLoadNode_" + str(i))
            if abs(fx) > _tol or abs(fy) > _tol or abs(fz) > _tol:
                model.ConcentratedForce(name="DistCF_" + str(i), createStepName="Step-1", region=region, cf1=fx, cf2=fy, cf3=fz)
            if abs(mx) > _tol or abs(my) > _tol or abs(mz) > _tol:
                model.Moment(name="DistMoment_" + str(i), createStepName="Step-1", region=region, cm1=mx, cm2=my, cm3=mz)
    except Exception as e:
        print("Distributed equivalent nodal load warning: " + str(e))
elif DISTRIBUTED_LOADS:
    try:
        fy_vals = [row[4] for row in DISTRIBUTED_LOADS if len(row) >= 5]
        if fy_vals:
            comp2_mag = sum(fy_vals) / len(fy_vals)
            beam_edges = assembly.instances["Beam-1"].edges
            region = assembly.Set(edges=beam_edges, name="BeamEdges")
            model.LineLoad(name="LineLoad-1", createStepName="Step-1", region=region, comp1=0.0, comp2=comp2_mag, comp3=0.0, distributionType=UNIFORM)
    except Exception as e:
        print("Distributed load warning: " + str(e))

# --- Point loads (apply at nearest node to given x,y,z) ---
for i, row in enumerate(POINT_LOADS):
    if len(row) < 9:
        continue
    x, y, z, Fx, Fy, Fz, Mx, My, Mz = row[0:9]
    idx = min(range(len(COORDS)), key=lambda i: (COORDS[i][0]-x)**2 + (COORDS[i][1]-y)**2 + (COORDS[i][2]-z)**2)
    node_label = idx + 1
    try:
        n_set = assembly.instances["Beam-1"].nodes.sequenceFromLabels(labels=(node_label,))
        region = assembly.Set(nodes=n_set, name="LoadNode_" + str(i))
        model.ConcentratedForce(name="CF_" + str(i), createStepName="Step-1", region=region, cf1=Fx, cf2=Fy, cf3=Fz)
        if Mx != 0 or My != 0 or Mz != 0:
            model.Moment(name="Moment_" + str(i), createStepName="Step-1", region=region, cm1=Mx, cm2=My, cm3=Mz)
    except Exception as e:
        print("Load warning: " + str(e))

# --- Job ---
jobName = "{job_name.replace('-', '_')}_abaqus"
mdb.Job(name=jobName, model=modelName, description="", type=ANALYSIS)
mdb.jobs[jobName].submit(consistencyChecking=OFF)
mdb.jobs[jobName].waitForCompletion()

'''


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate Abaqus CAE script from job dir")
    parser.add_argument("--job-dir", type=str, required=True, help="Path to job dir (e.g. jobs/job_0000_n8)")
    parser.add_argument("--output", type=str, default=None, help="Output script path (default: abaqus/generated/run_<job_name>.py)")
    args = parser.parse_args()

    job_dir = Path(args.job_dir)
    if not job_dir.is_absolute():
        job_dir = PROJECT_ROOT / job_dir
    if not job_dir.is_dir():
        print(f"Error: job directory not found: {job_dir}", file=sys.stderr)
        sys.exit(1)

    data = _parse_job(job_dir)
    job_name = data["job_name"]

    ABAQUS_GENERATED_DIR.mkdir(parents=True, exist_ok=True)
    out_script = Path(args.output) if args.output else ABAQUS_GENERATED_DIR / f"run_{job_name}.py"
    out_script = out_script.resolve()
    out_script.parent.mkdir(parents=True, exist_ok=True)

    out_csv_dir = str(ABAQUS_RESULTS_DIR / job_name)
    script_content = _generate_script_content(data, out_csv_dir)

    with open(out_script, "w", encoding="utf-8") as f:
        f.write(script_content)

    print(f"Generated: {out_script}")
    print(f"Abaqus results will be written to: {out_csv_dir}")


if __name__ == "__main__":
    main()

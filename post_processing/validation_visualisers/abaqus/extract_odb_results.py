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


def extract_odb_to_csv(odb_path: str, out_dir: str) -> None:
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
    odb_var_list = list(frame.fieldOutputs.keys())
    has_ur = "UR" in frame.fieldOutputs
    print("ODB field outputs:", odb_var_list)
    print("ODB has UR:", has_ur)

    # Displacements U (and UR if available) -> U_global.csv (same convention as FEM: 6 DOFs per node)
    if "U" in frame.fieldOutputs:
        u_field = frame.fieldOutputs["U"]
        node_list = list(u_field.values)
        # U is typically 3 components (U1,U2,U3); UR is separate in Abaqus
        ur_by_node = {}
        if "UR" in frame.fieldOutputs:
            for v in frame.fieldOutputs["UR"].values:
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
                if "UR" in frame.fieldOutputs:
                    ur1, ur2, ur3 = -ur1, -ur2, -ur3
                gdof_base = (node_label - 1) * 6
                for d, val in enumerate([u1, u2, u3, ur1, ur2, ur3]):
                    f.write(f"{gdof_base + d},{val}\n")
        rotation_src = "ODB" if "UR" in frame.fieldOutputs else "none"
        rotation_src_path = os.path.join(out_dir, "rotation_source.txt")
        with open(rotation_src_path, "w") as rf:
            rf.write(rotation_src + "\n")
        if has_ur:
            print("U_global.csv: U and UR read from ODB and written (rotation sign flipped to match FEM).")
        else:
            print("U_global.csv: U read from ODB; rotation not in ODB, written as zero.")
        print(f"Wrote {csv_path}")

    # Section forces: write section_forces.csv (x, N, Vy, Vz, T, My, Mz) and nodal_section_forces.csv (FEM format)
    if "SF" in frame.fieldOutputs:
        sf_field = frame.fieldOutputs["SF"]
        sm_field = frame.fieldOutputs.get("SM")
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
            elif len(data) >= 3 and sm_field:
                n, vy, vz = data[0], data[1], data[2]
                sm_vals = [v for v in sm_field.values if v.elementLabel == elabel]
                if sm_vals and len(sm_vals[0].data) >= 3:
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
                if conn and node_coords:
                    n1, n2 = conn[0], conn[1]
                    x_center = (node_coords[n1][0] + node_coords[n2][0]) / 2.0
                else:
                    x_center = float(elabel)
            except Exception:
                x_center = float(elabel)
            rows.append((x_center, n, vy, vz, t, my, mz))
        if rows:
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
        if node_labels_sorted and elem_conn and elem_sf:
            nodal_sum = defaultdict(lambda: [0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
            nodal_count = defaultdict(int)
            for elabel, (n, vy, vz, t, my, mz) in elem_sf.items():
                conn = elem_conn.get(elabel)
                if not conn:
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
    extract_odb_to_csv(odb_path, out_dir)


if __name__ == "__main__":
    main()

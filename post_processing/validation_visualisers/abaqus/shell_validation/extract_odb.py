# post_processing/validation_visualisers/abaqus/shell_validation/extract_odb.py
"""Write the station displacements of a shell validation ODB to CSV (run with Abaqus Python).

Usage::

    abaqus python extract_odb.py <job>.odb <out>.csv

Writes one row per node and frame: ``kind, frame, value, label, u1, u2, u3``. ``kind`` is
``buckle`` (``value`` is the eigenvalue of the mode) or ``static`` (``value`` is the
step time, which is the load factor). Only the ``STATIONS`` node set is written.

Kept to Python 2.7 syntax so that it also runs under the Abaqus releases that ship it.
"""

import re
import sys


def _node_set(odb, name):
    assembly = odb.rootAssembly
    if name in assembly.nodeSets.keys():
        return assembly.nodeSets[name]
    for inst in assembly.instances.values():
        if name in inst.nodeSets.keys():
            return inst.nodeSets[name]
    raise KeyError("node set %s not found in the ODB" % name)


def extract(odb_path, out_csv, set_name="STATIONS"):
    from odbAccess import openOdb

    odb = openOdb(path=odb_path, readOnly=True)
    region = _node_set(odb, set_name)
    step = odb.steps[list(odb.steps.keys())[-1]]
    buckle = "BUCKLE" in str(step.procedure).upper()
    out = open(out_csv, "w")
    out.write("kind,frame,value,label,u1,u2,u3\n")
    n = 0
    for i, frame in enumerate(step.frames):
        if buckle:
            m = re.search(r"EigenValue\s*=\s*([-+0-9.eE]+)", frame.description)
            if not m:
                continue          # the base state
            value = float(m.group(1))
        else:
            value = float(frame.frameValue)
        field = frame.fieldOutputs["U"].getSubset(region=region)
        for v in field.values:
            d = v.data
            out.write("%s,%d,%.12g,%d,%.12g,%.12g,%.12g\n"
                      % ("buckle" if buckle else "static", i, value, v.nodeLabel, d[0], d[1], d[2]))
        n += 1
    out.close()
    odb.close()
    print("wrote %d frames to %s" % (n, out_csv))


if __name__ == "__main__":
    args = [a for a in sys.argv[1:] if not a.startswith("-")]
    if len(args) < 2:
        print(__doc__)
        sys.exit(2)
    extract(args[0], args[1])

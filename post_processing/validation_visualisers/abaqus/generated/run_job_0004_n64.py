# -*- coding: utf-8 -*-
# Auto-generated Abaqus CAE script for job: job_0004_n64
# Run with project Python (abqpy): python this_file; abqpy saveAs() launches Abaqus.
# Outputs CSV to: C:\Users\s1834431\Programs\fem_model\post_processing\validation_visualisers\abaqus_results\job_0004_n64

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
COORDS = [[0.0, 0.0, 0.0], [0.03125, 0.0, 0.0], [0.0625, 0.0, 0.0], [0.09375, 0.0, 0.0], [0.125, 0.0, 0.0], [0.15625, 0.0, 0.0], [0.1875, 0.0, 0.0], [0.21875, 0.0, 0.0], [0.25, 0.0, 0.0], [0.28125, 0.0, 0.0], [0.3125, 0.0, 0.0], [0.34375, 0.0, 0.0], [0.375, 0.0, 0.0], [0.40625, 0.0, 0.0], [0.4375, 0.0, 0.0], [0.46875, 0.0, 0.0], [0.5, 0.0, 0.0], [0.53125, 0.0, 0.0], [0.5625, 0.0, 0.0], [0.59375, 0.0, 0.0], [0.625, 0.0, 0.0], [0.65625, 0.0, 0.0], [0.6875, 0.0, 0.0], [0.71875, 0.0, 0.0], [0.75, 0.0, 0.0], [0.78125, 0.0, 0.0], [0.8125, 0.0, 0.0], [0.84375, 0.0, 0.0], [0.875, 0.0, 0.0], [0.90625, 0.0, 0.0], [0.9375, 0.0, 0.0], [0.96875, 0.0, 0.0], [1.0, 0.0, 0.0], [1.03125, 0.0, 0.0], [1.0625, 0.0, 0.0], [1.09375, 0.0, 0.0], [1.125, 0.0, 0.0], [1.15625, 0.0, 0.0], [1.1875, 0.0, 0.0], [1.21875, 0.0, 0.0], [1.25, 0.0, 0.0], [1.28125, 0.0, 0.0], [1.3125, 0.0, 0.0], [1.34375, 0.0, 0.0], [1.375, 0.0, 0.0], [1.40625, 0.0, 0.0], [1.4375, 0.0, 0.0], [1.46875, 0.0, 0.0], [1.5, 0.0, 0.0], [1.53125, 0.0, 0.0], [1.5625, 0.0, 0.0], [1.59375, 0.0, 0.0], [1.625, 0.0, 0.0], [1.65625, 0.0, 0.0], [1.6875, 0.0, 0.0], [1.71875, 0.0, 0.0], [1.75, 0.0, 0.0], [1.78125, 0.0, 0.0], [1.8125, 0.0, 0.0], [1.84375, 0.0, 0.0], [1.875, 0.0, 0.0], [1.90625, 0.0, 0.0], [1.9375, 0.0, 0.0], [1.96875, 0.0, 0.0], [2.0, 0.0, 0.0]]
ABAQUS_ELEM = "B33"
E_MODULUS = 210000000000.0
NU = 0.3
RHO = 7850.0
A_AREA = 0.00131
I11 = 3.234e-07
I22 = 2.08769e-06
J_TORSION = 2.60673e-08
PRESCRIBED_NODES_DOF_VALUES = [(0, 0, 0.0), (0, 1, 0.0), (0, 2, 0.0), (0, 3, 0.0), (0, 4, 0.0), (0, 5, 0.0)]
POINT_LOADS = []
DISTRIBUTED_LOADS = [[0.0, 0.0, 0.0, 0.0, -0.0, 0.0, 0.0, 0.0, 0.0], [0.03125, 0.0, 0.0, 0.0, -7.8125, 0.0, 0.0, 0.0, 0.0], [0.0625, 0.0, 0.0, 0.0, -15.625, 0.0, 0.0, 0.0, 0.0], [0.09375, 0.0, 0.0, 0.0, -23.4375, 0.0, 0.0, 0.0, 0.0], [0.125, 0.0, 0.0, 0.0, -31.25, 0.0, 0.0, 0.0, 0.0], [0.15625, 0.0, 0.0, 0.0, -39.0625, 0.0, 0.0, 0.0, 0.0], [0.1875, 0.0, 0.0, 0.0, -46.875, 0.0, 0.0, 0.0, 0.0], [0.21875, 0.0, 0.0, 0.0, -54.6875, 0.0, 0.0, 0.0, 0.0], [0.25, 0.0, 0.0, 0.0, -62.5, 0.0, 0.0, 0.0, 0.0], [0.28125, 0.0, 0.0, 0.0, -70.3125, 0.0, 0.0, 0.0, 0.0], [0.3125, 0.0, 0.0, 0.0, -78.125, 0.0, 0.0, 0.0, 0.0], [0.34375, 0.0, 0.0, 0.0, -85.9375, 0.0, 0.0, 0.0, 0.0], [0.375, 0.0, 0.0, 0.0, -93.75, 0.0, 0.0, 0.0, 0.0], [0.40625, 0.0, 0.0, 0.0, -101.5625, 0.0, 0.0, 0.0, 0.0], [0.4375, 0.0, 0.0, 0.0, -109.375, 0.0, 0.0, 0.0, 0.0], [0.46875, 0.0, 0.0, 0.0, -117.1875, 0.0, 0.0, 0.0, 0.0], [0.5, 0.0, 0.0, 0.0, -125.0, 0.0, 0.0, 0.0, 0.0], [0.53125, 0.0, 0.0, 0.0, -132.8125, 0.0, 0.0, 0.0, 0.0], [0.5625, 0.0, 0.0, 0.0, -140.625, 0.0, 0.0, 0.0, 0.0], [0.59375, 0.0, 0.0, 0.0, -148.4375, 0.0, 0.0, 0.0, 0.0], [0.625, 0.0, 0.0, 0.0, -156.25, 0.0, 0.0, 0.0, 0.0], [0.65625, 0.0, 0.0, 0.0, -164.0625, 0.0, 0.0, 0.0, 0.0], [0.6875, 0.0, 0.0, 0.0, -171.875, 0.0, 0.0, 0.0, 0.0], [0.71875, 0.0, 0.0, 0.0, -179.6875, 0.0, 0.0, 0.0, 0.0], [0.75, 0.0, 0.0, 0.0, -187.5, 0.0, 0.0, 0.0, 0.0], [0.78125, 0.0, 0.0, 0.0, -195.3125, 0.0, 0.0, 0.0, 0.0], [0.8125, 0.0, 0.0, 0.0, -203.125, 0.0, 0.0, 0.0, 0.0], [0.84375, 0.0, 0.0, 0.0, -210.9375, 0.0, 0.0, 0.0, 0.0], [0.875, 0.0, 0.0, 0.0, -218.75, 0.0, 0.0, 0.0, 0.0], [0.90625, 0.0, 0.0, 0.0, -226.5625, 0.0, 0.0, 0.0, 0.0], [0.9375, 0.0, 0.0, 0.0, -234.375, 0.0, 0.0, 0.0, 0.0], [0.96875, 0.0, 0.0, 0.0, -242.1875, 0.0, 0.0, 0.0, 0.0], [1.0, 0.0, 0.0, 0.0, -250.0, 0.0, 0.0, 0.0, 0.0], [1.03125, 0.0, 0.0, 0.0, -257.8125, 0.0, 0.0, 0.0, 0.0], [1.0625, 0.0, 0.0, 0.0, -265.625, 0.0, 0.0, 0.0, 0.0], [1.09375, 0.0, 0.0, 0.0, -273.4375, 0.0, 0.0, 0.0, 0.0], [1.125, 0.0, 0.0, 0.0, -281.25, 0.0, 0.0, 0.0, 0.0], [1.15625, 0.0, 0.0, 0.0, -289.0625, 0.0, 0.0, 0.0, 0.0], [1.1875, 0.0, 0.0, 0.0, -296.875, 0.0, 0.0, 0.0, 0.0], [1.21875, 0.0, 0.0, 0.0, -304.6875, 0.0, 0.0, 0.0, 0.0], [1.25, 0.0, 0.0, 0.0, -312.5, 0.0, 0.0, 0.0, 0.0], [1.28125, 0.0, 0.0, 0.0, -320.3125, 0.0, 0.0, 0.0, 0.0], [1.3125, 0.0, 0.0, 0.0, -328.125, 0.0, 0.0, 0.0, 0.0], [1.34375, 0.0, 0.0, 0.0, -335.9375, 0.0, 0.0, 0.0, 0.0], [1.375, 0.0, 0.0, 0.0, -343.75, 0.0, 0.0, 0.0, 0.0], [1.40625, 0.0, 0.0, 0.0, -351.5625, 0.0, 0.0, 0.0, 0.0], [1.4375, 0.0, 0.0, 0.0, -359.375, 0.0, 0.0, 0.0, 0.0], [1.46875, 0.0, 0.0, 0.0, -367.1875, 0.0, 0.0, 0.0, 0.0], [1.5, 0.0, 0.0, 0.0, -375.0, 0.0, 0.0, 0.0, 0.0], [1.53125, 0.0, 0.0, 0.0, -382.8125, 0.0, 0.0, 0.0, 0.0], [1.5625, 0.0, 0.0, 0.0, -390.625, 0.0, 0.0, 0.0, 0.0], [1.59375, 0.0, 0.0, 0.0, -398.4375, 0.0, 0.0, 0.0, 0.0], [1.625, 0.0, 0.0, 0.0, -406.25, 0.0, 0.0, 0.0, 0.0], [1.65625, 0.0, 0.0, 0.0, -414.0625, 0.0, 0.0, 0.0, 0.0], [1.6875, 0.0, 0.0, 0.0, -421.875, 0.0, 0.0, 0.0, 0.0], [1.71875, 0.0, 0.0, 0.0, -429.6875, 0.0, 0.0, 0.0, 0.0], [1.75, 0.0, 0.0, 0.0, -437.5, 0.0, 0.0, 0.0, 0.0], [1.78125, 0.0, 0.0, 0.0, -445.3125, 0.0, 0.0, 0.0, 0.0], [1.8125, 0.0, 0.0, 0.0, -453.125, 0.0, 0.0, 0.0, 0.0], [1.84375, 0.0, 0.0, 0.0, -460.9375, 0.0, 0.0, 0.0, 0.0], [1.875, 0.0, 0.0, 0.0, -468.75, 0.0, 0.0, 0.0, 0.0], [1.90625, 0.0, 0.0, 0.0, -476.5625, 0.0, 0.0, 0.0, 0.0], [1.9375, 0.0, 0.0, 0.0, -484.375, 0.0, 0.0, 0.0, 0.0], [1.96875, 0.0, 0.0, 0.0, -492.1875, 0.0, 0.0, 0.0, 0.0], [2.0, 0.0, 0.0, 0.0, -500.0, 0.0, 0.0, 0.0, 0.0]]
DISTRIBUTED_EQUIVALENT_NODAL = [(1, 0.0, -0.036621093750000014, 0.0, 0.0, 0.0, 0.0), (2, 0.0, -0.24414062500000006, 0.0, 0.0, 0.0, 0.0), (3, 0.0, -0.4882812500000001, 0.0, 0.0, 0.0, 0.0), (4, 0.0, -0.7324218750000001, 0.0, 0.0, 0.0, 0.0), (5, 0.0, -0.9765625000000001, 0.0, 0.0, 0.0, 0.0), (6, 0.0, -1.220703125, 0.0, 0.0, 0.0, 0.0), (7, 0.0, -1.4648437500000002, 0.0, 0.0, 0.0, 0.0), (8, 0.0, -1.708984375, 0.0, 0.0, 0.0, 0.0), (9, 0.0, -1.9531250000000002, 0.0, 0.0, 0.0, 0.0), (10, 0.0, -2.197265625, 0.0, 0.0, 0.0, 0.0), (11, 0.0, -2.44140625, 0.0, 0.0, 0.0, 0.0), (12, 0.0, -2.685546875000001, 0.0, 0.0, 0.0, 0.0), (13, 0.0, -2.929687500000001, 0.0, 0.0, 0.0, 0.0), (14, 0.0, -3.1738281250000004, 0.0, 0.0, 0.0, 0.0), (15, 0.0, -3.417968750000001, 0.0, 0.0, 0.0, 0.0), (16, 0.0, -3.662109375000001, 0.0, 0.0, 0.0, 0.0), (17, 0.0, -3.906250000000001, 0.0, 0.0, 0.0, 0.0), (18, 0.0, -4.150390625000002, 0.0, 0.0, 0.0, 0.0), (19, 0.0, -4.394531250000001, 0.0, 0.0, 0.0, 0.0), (20, 0.0, -4.638671875, 0.0, 0.0, 0.0, 0.0), (21, 0.0, -4.882812500000002, 0.0, 0.0, 0.0, 0.0), (22, 0.0, -5.126953125000001, 0.0, 0.0, 0.0, 0.0), (23, 0.0, -5.371093750000002, 0.0, 0.0, 0.0, 0.0), (24, 0.0, -5.615234375000002, 0.0, 0.0, 0.0, 0.0), (25, 0.0, -5.859375, 0.0, 0.0, 0.0, 0.0), (26, 0.0, -6.103515625000002, 0.0, 0.0, 0.0, 0.0), (27, 0.0, -6.347656250000002, 0.0, 0.0, 0.0, 0.0), (28, 0.0, -6.591796875, 0.0, 0.0, 0.0, 0.0), (29, 0.0, -6.835937500000002, 0.0, 0.0, 0.0, 0.0), (30, 0.0, -7.080078125000002, 0.0, 0.0, 0.0, 0.0), (31, 0.0, -7.324218750000002, 0.0, 0.0, 0.0, 0.0), (32, 0.0, -7.568359375000002, 0.0, 0.0, 0.0, 0.0), (33, 0.0, -7.812500000000002, 0.0, 0.0, 0.0, 0.0), (34, 0.0, -8.056640625, 0.0, 0.0, 0.0, 0.0), (35, 0.0, -8.300781250000002, 0.0, 0.0, 0.0, 0.0), (36, 0.0, -8.544921875000002, 0.0, 0.0, 0.0, 0.0), (37, 0.0, -8.789062500000002, 0.0, 0.0, 0.0, 0.0), (38, 0.0, -9.033203125, 0.0, 0.0, 0.0, 0.0), (39, 0.0, -9.27734375, 0.0, 0.0, 0.0, 0.0), (40, 0.0, -9.521484375000002, 0.0, 0.0, 0.0, 0.0), (41, 0.0, -9.765625000000002, 0.0, 0.0, 0.0, 0.0), (42, 0.0, -10.009765625000002, 0.0, 0.0, 0.0, 0.0), (43, 0.0, -10.253906250000004, 0.0, 0.0, 0.0, 0.0), (44, 0.0, -10.498046875000004, 0.0, 0.0, 0.0, 0.0), (45, 0.0, -10.7421875, 0.0, 0.0, 0.0, 0.0), (46, 0.0, -10.986328125000002, 0.0, 0.0, 0.0, 0.0), (47, 0.0, -11.230468750000002, 0.0, 0.0, 0.0, 0.0), (48, 0.0, -11.474609375000002, 0.0, 0.0, 0.0, 0.0), (49, 0.0, -11.718750000000004, 0.0, 0.0, 0.0, 0.0), (50, 0.0, -11.962890625000004, 0.0, 0.0, 0.0, 0.0), (51, 0.0, -12.20703125, 0.0, 0.0, 0.0, 0.0), (52, 0.0, -12.451171875000002, 0.0, 0.0, 0.0, 0.0), (53, 0.0, -12.6953125, 0.0, 0.0, 0.0, 0.0), (54, 0.0, -12.939453125000002, 0.0, 0.0, 0.0, 0.0), (55, 0.0, -13.183593750000002, 0.0, 0.0, 0.0, 0.0), (56, 0.0, -13.427734375, 0.0, 0.0, 0.0, 0.0), (57, 0.0, -13.671875000000004, 0.0, 0.0, 0.0, 0.0), (58, 0.0, -13.916015625, 0.0, 0.0, 0.0, 0.0), (59, 0.0, -14.16015625, 0.0, 0.0, 0.0, 0.0), (60, 0.0, -14.404296875000002, 0.0, 0.0, 0.0, 0.0), (61, 0.0, -14.648437500000004, 0.0, 0.0, 0.0, 0.0), (62, 0.0, -14.892578125000002, 0.0, 0.0, 0.0, 0.0), (63, 0.0, -15.13671875, 0.0, 0.0, 0.0, 0.0), (64, 0.0, -15.380859375, 0.0, 0.0, 0.0, 0.0), (65, 0.0, -7.77587890625, 0.0, 0.0, 0.0, 0.0)]
OUT_CSV_DIR = r"C:\\Users\\s1834431\\Programs\\fem_model\\post_processing\\validation_visualisers\\abaqus_results\\job_0004_n64"

# Default encastre at node 0 if no prescribed DOFs given
if not PRESCRIBED_NODES_DOF_VALUES and COORDS:
    for dof_idx in range(6):
        PRESCRIBED_NODES_DOF_VALUES.append((0, dof_idx, 0.0))

# When run by project Python (abqpy), only trigger Abaqus launch; model/job/export run only inside Abaqus
if not IN_ABAQUS:
    mdb.saveAs(os.path.join(OUT_CSV_DIR, "model.cae"))
    sys.exit(0)

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
model.StaticStep(name="Step-1", previous="Initial", description="Static")
if not os.path.exists(OUT_CSV_DIR):
    os.makedirs(OUT_CSV_DIR)
def _log(msg):
    print(msg)
    try:
        with open(os.path.join(OUT_CSV_DIR, "run_log.txt"), "a") as _f:
            _f.write(str(msg) + "\n")
    except Exception:
        pass
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
jobName = "job_0004_n64_abaqus"
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
    ur_by_node = {}
    if "UR" in frame.fieldOutputs:
        for _v in frame.fieldOutputs["UR"].values:
            nlab = _v.nodeLabel
            ur_by_node[nlab] = (_v.data[0], _v.data[1], _v.data[2]) if len(_v.data) >= 3 else (0.0, 0.0, 0.0)
    if not os.path.exists(OUT_CSV_DIR):
        os.makedirs(OUT_CSV_DIR)
    csv_path = os.path.join(OUT_CSV_DIR, "U_global.csv")
    with open(csv_path, "w") as f:
        f.write("Global DOF,Value\n")
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
                f.write(str(gdof_base + d) + "," + str(val) + "\n")
    rotation_src = "ODB" if "UR" in frame.fieldOutputs else "none"
    with open(os.path.join(OUT_CSV_DIR, "rotation_source.txt"), "w") as _rf:
        _rf.write(rotation_src + "\n")
    # --- Export section forces (x, N, Vy, Vz, T, My, Mz) for comparison ---
    if "SF" in frame.fieldOutputs:
        sf_field = frame.fieldOutputs["SF"]
        sm_field = frame.fieldOutputs["SM"] if "SM" in frame.fieldOutputs else None
        sf_path = os.path.join(OUT_CSV_DIR, "section_forces.csv")
        with open(sf_path, "w") as sf_file:
            sf_file.write("x,N,Vy,Vz,T,My,Mz\n")
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
                sf_file.write(str(x_center) + "," + str(n) + "," + str(vy) + "," + str(vz) + "," + str(t) + "," + str(my) + "," + str(mz) + "\n")
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

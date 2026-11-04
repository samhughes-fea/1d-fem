# -*- coding: utf-8 -*-
# Auto-generated Abaqus CAE script for job: job_0011_n64
# Run with project Python (abqpy): python this_file; abqpy saveAs() launches Abaqus.
# Outputs CSV to: C:\Users\s1834431\Programs\fem_model\post_processing\validation_visualisers\abaqus_results\job_0011_n64

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
ABAQUS_ELEM = "B31"
E_MODULUS = 210000000000.0
NU = 0.3
RHO = 7850.0
A_AREA = 0.00131
I11 = 3.234e-07
I22 = 2.08769e-06
J_TORSION = 2.60673e-08
PRESCRIBED_NODES_DOF_VALUES = [(0, 0, 0.0), (0, 1, 0.0), (0, 2, 0.0), (0, 3, 0.0), (0, 4, 0.0), (0, 5, 0.0)]
POINT_LOADS = []
DISTRIBUTED_LOADS = [[0.0, 0.0, 0.0, 0.0, -0.0, 0.0, 0.0, 0.0, 0.0], [0.03125, 0.0, 0.0, 0.0, -0.12207, 0.0, 0.0, 0.0, 0.0], [0.0625, 0.0, 0.0, 0.0, -0.488281, 0.0, 0.0, 0.0, 0.0], [0.09375, 0.0, 0.0, 0.0, -1.098633, 0.0, 0.0, 0.0, 0.0], [0.125, 0.0, 0.0, 0.0, -1.953125, 0.0, 0.0, 0.0, 0.0], [0.15625, 0.0, 0.0, 0.0, -3.051758, 0.0, 0.0, 0.0, 0.0], [0.1875, 0.0, 0.0, 0.0, -4.394531, 0.0, 0.0, 0.0, 0.0], [0.21875, 0.0, 0.0, 0.0, -5.981445, 0.0, 0.0, 0.0, 0.0], [0.25, 0.0, 0.0, 0.0, -7.8125, 0.0, 0.0, 0.0, 0.0], [0.28125, 0.0, 0.0, 0.0, -9.887695, 0.0, 0.0, 0.0, 0.0], [0.3125, 0.0, 0.0, 0.0, -12.207031, 0.0, 0.0, 0.0, 0.0], [0.34375, 0.0, 0.0, 0.0, -14.770508, 0.0, 0.0, 0.0, 0.0], [0.375, 0.0, 0.0, 0.0, -17.578125, 0.0, 0.0, 0.0, 0.0], [0.40625, 0.0, 0.0, 0.0, -20.629883, 0.0, 0.0, 0.0, 0.0], [0.4375, 0.0, 0.0, 0.0, -23.925781, 0.0, 0.0, 0.0, 0.0], [0.46875, 0.0, 0.0, 0.0, -27.46582, 0.0, 0.0, 0.0, 0.0], [0.5, 0.0, 0.0, 0.0, -31.25, 0.0, 0.0, 0.0, 0.0], [0.53125, 0.0, 0.0, 0.0, -35.27832, 0.0, 0.0, 0.0, 0.0], [0.5625, 0.0, 0.0, 0.0, -39.550781, 0.0, 0.0, 0.0, 0.0], [0.59375, 0.0, 0.0, 0.0, -44.067383, 0.0, 0.0, 0.0, 0.0], [0.625, 0.0, 0.0, 0.0, -48.828125, 0.0, 0.0, 0.0, 0.0], [0.65625, 0.0, 0.0, 0.0, -53.833008, 0.0, 0.0, 0.0, 0.0], [0.6875, 0.0, 0.0, 0.0, -59.082031, 0.0, 0.0, 0.0, 0.0], [0.71875, 0.0, 0.0, 0.0, -64.575195, 0.0, 0.0, 0.0, 0.0], [0.75, 0.0, 0.0, 0.0, -70.3125, 0.0, 0.0, 0.0, 0.0], [0.78125, 0.0, 0.0, 0.0, -76.293945, 0.0, 0.0, 0.0, 0.0], [0.8125, 0.0, 0.0, 0.0, -82.519531, 0.0, 0.0, 0.0, 0.0], [0.84375, 0.0, 0.0, 0.0, -88.989258, 0.0, 0.0, 0.0, 0.0], [0.875, 0.0, 0.0, 0.0, -95.703125, 0.0, 0.0, 0.0, 0.0], [0.90625, 0.0, 0.0, 0.0, -102.661133, 0.0, 0.0, 0.0, 0.0], [0.9375, 0.0, 0.0, 0.0, -109.863281, 0.0, 0.0, 0.0, 0.0], [0.96875, 0.0, 0.0, 0.0, -117.30957, 0.0, 0.0, 0.0, 0.0], [1.0, 0.0, 0.0, 0.0, -125.0, 0.0, 0.0, 0.0, 0.0], [1.03125, 0.0, 0.0, 0.0, -132.93457, 0.0, 0.0, 0.0, 0.0], [1.0625, 0.0, 0.0, 0.0, -141.113281, 0.0, 0.0, 0.0, 0.0], [1.09375, 0.0, 0.0, 0.0, -149.536133, 0.0, 0.0, 0.0, 0.0], [1.125, 0.0, 0.0, 0.0, -158.203125, 0.0, 0.0, 0.0, 0.0], [1.15625, 0.0, 0.0, 0.0, -167.114258, 0.0, 0.0, 0.0, 0.0], [1.1875, 0.0, 0.0, 0.0, -176.269531, 0.0, 0.0, 0.0, 0.0], [1.21875, 0.0, 0.0, 0.0, -185.668945, 0.0, 0.0, 0.0, 0.0], [1.25, 0.0, 0.0, 0.0, -195.3125, 0.0, 0.0, 0.0, 0.0], [1.28125, 0.0, 0.0, 0.0, -205.200195, 0.0, 0.0, 0.0, 0.0], [1.3125, 0.0, 0.0, 0.0, -215.332031, 0.0, 0.0, 0.0, 0.0], [1.34375, 0.0, 0.0, 0.0, -225.708008, 0.0, 0.0, 0.0, 0.0], [1.375, 0.0, 0.0, 0.0, -236.328125, 0.0, 0.0, 0.0, 0.0], [1.40625, 0.0, 0.0, 0.0, -247.192383, 0.0, 0.0, 0.0, 0.0], [1.4375, 0.0, 0.0, 0.0, -258.300781, 0.0, 0.0, 0.0, 0.0], [1.46875, 0.0, 0.0, 0.0, -269.65332, 0.0, 0.0, 0.0, 0.0], [1.5, 0.0, 0.0, 0.0, -281.25, 0.0, 0.0, 0.0, 0.0], [1.53125, 0.0, 0.0, 0.0, -293.09082, 0.0, 0.0, 0.0, 0.0], [1.5625, 0.0, 0.0, 0.0, -305.175781, 0.0, 0.0, 0.0, 0.0], [1.59375, 0.0, 0.0, 0.0, -317.504883, 0.0, 0.0, 0.0, 0.0], [1.625, 0.0, 0.0, 0.0, -330.078125, 0.0, 0.0, 0.0, 0.0], [1.65625, 0.0, 0.0, 0.0, -342.895508, 0.0, 0.0, 0.0, 0.0], [1.6875, 0.0, 0.0, 0.0, -355.957031, 0.0, 0.0, 0.0, 0.0], [1.71875, 0.0, 0.0, 0.0, -369.262695, 0.0, 0.0, 0.0, 0.0], [1.75, 0.0, 0.0, 0.0, -382.8125, 0.0, 0.0, 0.0, 0.0], [1.78125, 0.0, 0.0, 0.0, -396.606445, 0.0, 0.0, 0.0, 0.0], [1.8125, 0.0, 0.0, 0.0, -410.644531, 0.0, 0.0, 0.0, 0.0], [1.84375, 0.0, 0.0, 0.0, -424.926758, 0.0, 0.0, 0.0, 0.0], [1.875, 0.0, 0.0, 0.0, -439.453125, 0.0, 0.0, 0.0, 0.0], [1.90625, 0.0, 0.0, 0.0, -454.223633, 0.0, 0.0, 0.0, 0.0], [1.9375, 0.0, 0.0, 0.0, -469.238281, 0.0, 0.0, 0.0, 0.0], [1.96875, 0.0, 0.0, 0.0, -484.49707, 0.0, 0.0, 0.0, 0.0], [2.0, 0.0, 0.0, 0.0, -500.0, 0.0, 0.0, 0.0, 0.0]]
DISTRIBUTED_EQUIVALENT_NODAL = [(1, 0.0, -0.00025431152214625496, 0.0, 0.0, 0.0, 0.0), (2, 0.0, -0.004323314583333334, 0.0, 0.0, 0.0, 0.0), (3, 0.0, -0.015767408789656462, 0.0, 0.0, 0.0, 0.0), (4, 0.0, -0.034840905414290836, 0.0, 0.0, 0.0, 0.0), (5, 0.0, -0.061543784240680216, 0.0, 0.0, 0.0, 0.0), (6, 0.0, -0.09587606168548835, 0.0, 0.0, 0.0, 0.0), (7, 0.0, -0.13783772120486643, 0.0, 0.0, 0.0, 0.0), (8, 0.0, -0.18742878365129598, 0.0, 0.0, 0.0, 0.0), (9, 0.0, -0.24464924934619972, 0.0, 0.0, 0.0, 0.0), (10, 0.0, -0.3094990961514053, 0.0, 0.0, 0.0, 0.0), (11, 0.0, -0.38197834620442944, 0.0, 0.0, 0.0, 0.0), (12, 0.0, -0.4620869991871275, 0.0, 0.0, 0.0, 0.0), (13, 0.0, -0.5498250342345612, 0.0, 0.0, 0.0, 0.0), (14, 0.0, -0.6451924679371279, 0.0, 0.0, 0.0, 0.0), (15, 0.0, -0.7481892837044273, 0.0, 0.0, 0.0, 0.0), (16, 0.0, -0.8588155024014138, 0.0, 0.0, 0.0, 0.0), (17, 0.0, -0.9770711243461685, 0.0, 0.0, 0.0, 0.0), (18, 0.0, -1.1029561274014141, 0.0, 0.0, 0.0, 0.0), (19, 0.0, -1.2364705337044273, 0.0, 0.0, 0.0, 0.0), (20, 0.0, -1.3776143429371284, 0.0, 0.0, 0.0, 0.0), (21, 0.0, -1.5263875342345616, 0.0, 0.0, 0.0, 0.0), (22, 0.0, -1.6827901241871284, 0.0, 0.0, 0.0, 0.0), (23, 0.0, -1.8468220962044275, 0.0, 0.0, 0.0, 0.0), (24, 0.0, -2.018483471151414, 0.0, 0.0, 0.0, 0.0), (25, 0.0, -2.197774249346169, 0.0, 0.0, 0.0, 0.0), (26, 0.0, -2.384694408651414, 0.0, 0.0, 0.0, 0.0), (27, 0.0, -2.5792439712044275, 0.0, 0.0, 0.0, 0.0), (28, 0.0, -2.781422936687129, 0.0, 0.0, 0.0, 0.0), (29, 0.0, -2.991231284234562, 0.0, 0.0, 0.0, 0.0), (30, 0.0, -3.208669030437129, 0.0, 0.0, 0.0, 0.0), (31, 0.0, -3.433736158704428, 0.0, 0.0, 0.0, 0.0), (32, 0.0, -3.6664326899014146, 0.0, 0.0, 0.0, 0.0), (33, 0.0, -3.9067586243461685, 0.0, 0.0, 0.0, 0.0), (34, 0.0, -4.154713939901415, 0.0, 0.0, 0.0, 0.0), (35, 0.0, -4.410298658704428, 0.0, 0.0, 0.0, 0.0), (36, 0.0, -4.673512780437129, 0.0, 0.0, 0.0, 0.0), (37, 0.0, -4.944356284234562, 0.0, 0.0, 0.0, 0.0), (38, 0.0, -5.22282918668713, 0.0, 0.0, 0.0, 0.0), (39, 0.0, -5.508931471204428, 0.0, 0.0, 0.0, 0.0), (40, 0.0, -5.802663158651415, 0.0, 0.0, 0.0, 0.0), (41, 0.0, -6.104024249346169, 0.0, 0.0, 0.0, 0.0), (42, 0.0, -6.413014721151415, 0.0, 0.0, 0.0, 0.0), (43, 0.0, -6.729634596204429, 0.0, 0.0, 0.0, 0.0), (44, 0.0, -7.053883874187129, 0.0, 0.0, 0.0, 0.0), (45, 0.0, -7.385762534234562, 0.0, 0.0, 0.0, 0.0), (46, 0.0, -7.725270592937128, 0.0, 0.0, 0.0, 0.0), (47, 0.0, -8.072408033704427, 0.0, 0.0, 0.0, 0.0), (48, 0.0, -8.427174877401415, 0.0, 0.0, 0.0, 0.0), (49, 0.0, -8.789571124346171, 0.0, 0.0, 0.0, 0.0), (50, 0.0, -9.159596752401415, 0.0, 0.0, 0.0, 0.0), (51, 0.0, -9.537251783704427, 0.0, 0.0, 0.0, 0.0), (52, 0.0, -9.92253621793713, 0.0, 0.0, 0.0, 0.0), (53, 0.0, -10.315450034234562, 0.0, 0.0, 0.0, 0.0), (54, 0.0, -10.715993249187129, 0.0, 0.0, 0.0, 0.0), (55, 0.0, -11.12416584620443, 0.0, 0.0, 0.0, 0.0), (56, 0.0, -11.539967846151406, 0.0, 0.0, 0.0, 0.0), (57, 0.0, -11.963399249346203, 0.0, 0.0, 0.0, 0.0), (58, 0.0, -12.3944600336513, 0.0, 0.0, 0.0, 0.0), (59, 0.0, -12.833150221204868, 0.0, 0.0, 0.0, 0.0), (60, 0.0, -13.279469811685491, 0.0, 0.0, 0.0, 0.0), (61, 0.0, -13.733418784240683, 0.0, 0.0, 0.0, 0.0), (62, 0.0, -14.194997155414292, 0.0, 0.0, 0.0, 0.0), (63, 0.0, -14.664204908789658, 0.0, 0.0, 0.0, 0.0), (64, 0.0, -15.141042064583337, 0.0, 0.0, 0.0, 0.0), (65, 0.0, -7.7395121240221485, 0.0, 0.0, 0.0, 0.0)]
OUT_CSV_DIR = r"C:\\Users\\s1834431\\Programs\\fem_model\\post_processing\\validation_visualisers\\abaqus_results\\job_0011_n64"

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
jobName = "job_0011_n64_abaqus"
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

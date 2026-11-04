# -*- coding: utf-8 -*-
# Auto-generated Abaqus CAE script for job: job_0008_n500
# Run with project Python (abqpy): python this_file; abqpy saveAs() launches Abaqus.
# Outputs CSV to: C:\Users\s1834431\Programs\fem_model\post_processing\validation_visualisers\abaqus_results\job_0008_n500

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
COORDS = [[0.0, 0.0, 0.0], [0.004, 0.0, 0.0], [0.008, 0.0, 0.0], [0.012, 0.0, 0.0], [0.016, 0.0, 0.0], [0.02, 0.0, 0.0], [0.024, 0.0, 0.0], [0.028, 0.0, 0.0], [0.032, 0.0, 0.0], [0.036, 0.0, 0.0], [0.04, 0.0, 0.0], [0.044, 0.0, 0.0], [0.048, 0.0, 0.0], [0.052, 0.0, 0.0], [0.056, 0.0, 0.0], [0.06, 0.0, 0.0], [0.064, 0.0, 0.0], [0.068, 0.0, 0.0], [0.072, 0.0, 0.0], [0.076, 0.0, 0.0], [0.08, 0.0, 0.0], [0.084, 0.0, 0.0], [0.088, 0.0, 0.0], [0.092, 0.0, 0.0], [0.096, 0.0, 0.0], [0.1, 0.0, 0.0], [0.104, 0.0, 0.0], [0.108, 0.0, 0.0], [0.112, 0.0, 0.0], [0.116, 0.0, 0.0], [0.12, 0.0, 0.0], [0.124, 0.0, 0.0], [0.128, 0.0, 0.0], [0.132, 0.0, 0.0], [0.136, 0.0, 0.0], [0.14, 0.0, 0.0], [0.144, 0.0, 0.0], [0.148, 0.0, 0.0], [0.152, 0.0, 0.0], [0.156, 0.0, 0.0], [0.16, 0.0, 0.0], [0.164, 0.0, 0.0], [0.168, 0.0, 0.0], [0.172, 0.0, 0.0], [0.176, 0.0, 0.0], [0.18, 0.0, 0.0], [0.184, 0.0, 0.0], [0.188, 0.0, 0.0], [0.192, 0.0, 0.0], [0.196, 0.0, 0.0], [0.2, 0.0, 0.0], [0.204, 0.0, 0.0], [0.208, 0.0, 0.0], [0.212, 0.0, 0.0], [0.216, 0.0, 0.0], [0.22, 0.0, 0.0], [0.224, 0.0, 0.0], [0.228, 0.0, 0.0], [0.232, 0.0, 0.0], [0.236, 0.0, 0.0], [0.24, 0.0, 0.0], [0.244, 0.0, 0.0], [0.248, 0.0, 0.0], [0.252, 0.0, 0.0], [0.256, 0.0, 0.0], [0.26, 0.0, 0.0], [0.264, 0.0, 0.0], [0.268, 0.0, 0.0], [0.272, 0.0, 0.0], [0.276, 0.0, 0.0], [0.28, 0.0, 0.0], [0.284, 0.0, 0.0], [0.288, 0.0, 0.0], [0.292, 0.0, 0.0], [0.296, 0.0, 0.0], [0.3, 0.0, 0.0], [0.304, 0.0, 0.0], [0.308, 0.0, 0.0], [0.312, 0.0, 0.0], [0.316, 0.0, 0.0], [0.32, 0.0, 0.0], [0.324, 0.0, 0.0], [0.328, 0.0, 0.0], [0.332, 0.0, 0.0], [0.336, 0.0, 0.0], [0.34, 0.0, 0.0], [0.344, 0.0, 0.0], [0.348, 0.0, 0.0], [0.352, 0.0, 0.0], [0.356, 0.0, 0.0], [0.36, 0.0, 0.0], [0.364, 0.0, 0.0], [0.368, 0.0, 0.0], [0.372, 0.0, 0.0], [0.376, 0.0, 0.0], [0.38, 0.0, 0.0], [0.384, 0.0, 0.0], [0.388, 0.0, 0.0], [0.392, 0.0, 0.0], [0.396, 0.0, 0.0], [0.4, 0.0, 0.0], [0.404, 0.0, 0.0], [0.408, 0.0, 0.0], [0.412, 0.0, 0.0], [0.416, 0.0, 0.0], [0.42, 0.0, 0.0], [0.424, 0.0, 0.0], [0.428, 0.0, 0.0], [0.432, 0.0, 0.0], [0.436, 0.0, 0.0], [0.44, 0.0, 0.0], [0.444, 0.0, 0.0], [0.448, 0.0, 0.0], [0.452, 0.0, 0.0], [0.456, 0.0, 0.0], [0.46, 0.0, 0.0], [0.464, 0.0, 0.0], [0.468, 0.0, 0.0], [0.472, 0.0, 0.0], [0.476, 0.0, 0.0], [0.48, 0.0, 0.0], [0.484, 0.0, 0.0], [0.488, 0.0, 0.0], [0.492, 0.0, 0.0], [0.496, 0.0, 0.0], [0.5, 0.0, 0.0], [0.504, 0.0, 0.0], [0.508, 0.0, 0.0], [0.512, 0.0, 0.0], [0.516, 0.0, 0.0], [0.52, 0.0, 0.0], [0.524, 0.0, 0.0], [0.528, 0.0, 0.0], [0.532, 0.0, 0.0], [0.536, 0.0, 0.0], [0.54, 0.0, 0.0], [0.544, 0.0, 0.0], [0.548, 0.0, 0.0], [0.552, 0.0, 0.0], [0.556, 0.0, 0.0], [0.56, 0.0, 0.0], [0.564, 0.0, 0.0], [0.568, 0.0, 0.0], [0.572, 0.0, 0.0], [0.576, 0.0, 0.0], [0.58, 0.0, 0.0], [0.584, 0.0, 0.0], [0.588, 0.0, 0.0], [0.592, 0.0, 0.0], [0.596, 0.0, 0.0], [0.6, 0.0, 0.0], [0.604, 0.0, 0.0], [0.608, 0.0, 0.0], [0.612, 0.0, 0.0], [0.616, 0.0, 0.0], [0.62, 0.0, 0.0], [0.624, 0.0, 0.0], [0.628, 0.0, 0.0], [0.632, 0.0, 0.0], [0.636, 0.0, 0.0], [0.64, 0.0, 0.0], [0.644, 0.0, 0.0], [0.648, 0.0, 0.0], [0.652, 0.0, 0.0], [0.656, 0.0, 0.0], [0.66, 0.0, 0.0], [0.664, 0.0, 0.0], [0.668, 0.0, 0.0], [0.672, 0.0, 0.0], [0.676, 0.0, 0.0], [0.68, 0.0, 0.0], [0.684, 0.0, 0.0], [0.688, 0.0, 0.0], [0.692, 0.0, 0.0], [0.696, 0.0, 0.0], [0.7, 0.0, 0.0], [0.704, 0.0, 0.0], [0.708, 0.0, 0.0], [0.712, 0.0, 0.0], [0.716, 0.0, 0.0], [0.72, 0.0, 0.0], [0.724, 0.0, 0.0], [0.728, 0.0, 0.0], [0.732, 0.0, 0.0], [0.736, 0.0, 0.0], [0.74, 0.0, 0.0], [0.744, 0.0, 0.0], [0.748, 0.0, 0.0], [0.752, 0.0, 0.0], [0.756, 0.0, 0.0], [0.76, 0.0, 0.0], [0.764, 0.0, 0.0], [0.768, 0.0, 0.0], [0.772, 0.0, 0.0], [0.776, 0.0, 0.0], [0.78, 0.0, 0.0], [0.784, 0.0, 0.0], [0.788, 0.0, 0.0], [0.792, 0.0, 0.0], [0.796, 0.0, 0.0], [0.8, 0.0, 0.0], [0.804, 0.0, 0.0], [0.808, 0.0, 0.0], [0.812, 0.0, 0.0], [0.816, 0.0, 0.0], [0.82, 0.0, 0.0], [0.824, 0.0, 0.0], [0.828, 0.0, 0.0], [0.832, 0.0, 0.0], [0.836, 0.0, 0.0], [0.84, 0.0, 0.0], [0.844, 0.0, 0.0], [0.848, 0.0, 0.0], [0.852, 0.0, 0.0], [0.856, 0.0, 0.0], [0.86, 0.0, 0.0], [0.864, 0.0, 0.0], [0.868, 0.0, 0.0], [0.872, 0.0, 0.0], [0.876, 0.0, 0.0], [0.88, 0.0, 0.0], [0.884, 0.0, 0.0], [0.888, 0.0, 0.0], [0.892, 0.0, 0.0], [0.896, 0.0, 0.0], [0.9, 0.0, 0.0], [0.904, 0.0, 0.0], [0.908, 0.0, 0.0], [0.912, 0.0, 0.0], [0.916, 0.0, 0.0], [0.92, 0.0, 0.0], [0.924, 0.0, 0.0], [0.928, 0.0, 0.0], [0.932, 0.0, 0.0], [0.936, 0.0, 0.0], [0.94, 0.0, 0.0], [0.944, 0.0, 0.0], [0.948, 0.0, 0.0], [0.952, 0.0, 0.0], [0.956, 0.0, 0.0], [0.96, 0.0, 0.0], [0.964, 0.0, 0.0], [0.968, 0.0, 0.0], [0.972, 0.0, 0.0], [0.976, 0.0, 0.0], [0.98, 0.0, 0.0], [0.984, 0.0, 0.0], [0.988, 0.0, 0.0], [0.992, 0.0, 0.0], [0.996, 0.0, 0.0], [1.0, 0.0, 0.0], [1.004, 0.0, 0.0], [1.008, 0.0, 0.0], [1.012, 0.0, 0.0], [1.016, 0.0, 0.0], [1.02, 0.0, 0.0], [1.024, 0.0, 0.0], [1.028, 0.0, 0.0], [1.032, 0.0, 0.0], [1.036, 0.0, 0.0], [1.04, 0.0, 0.0], [1.044, 0.0, 0.0], [1.048, 0.0, 0.0], [1.052, 0.0, 0.0], [1.056, 0.0, 0.0], [1.06, 0.0, 0.0], [1.064, 0.0, 0.0], [1.068, 0.0, 0.0], [1.072, 0.0, 0.0], [1.076, 0.0, 0.0], [1.08, 0.0, 0.0], [1.084, 0.0, 0.0], [1.088, 0.0, 0.0], [1.092, 0.0, 0.0], [1.096, 0.0, 0.0], [1.1, 0.0, 0.0], [1.104, 0.0, 0.0], [1.108, 0.0, 0.0], [1.112, 0.0, 0.0], [1.116, 0.0, 0.0], [1.12, 0.0, 0.0], [1.124, 0.0, 0.0], [1.128, 0.0, 0.0], [1.132, 0.0, 0.0], [1.136, 0.0, 0.0], [1.14, 0.0, 0.0], [1.144, 0.0, 0.0], [1.148, 0.0, 0.0], [1.152, 0.0, 0.0], [1.156, 0.0, 0.0], [1.16, 0.0, 0.0], [1.164, 0.0, 0.0], [1.168, 0.0, 0.0], [1.172, 0.0, 0.0], [1.176, 0.0, 0.0], [1.18, 0.0, 0.0], [1.184, 0.0, 0.0], [1.188, 0.0, 0.0], [1.192, 0.0, 0.0], [1.196, 0.0, 0.0], [1.2, 0.0, 0.0], [1.204, 0.0, 0.0], [1.208, 0.0, 0.0], [1.212, 0.0, 0.0], [1.216, 0.0, 0.0], [1.22, 0.0, 0.0], [1.224, 0.0, 0.0], [1.228, 0.0, 0.0], [1.232, 0.0, 0.0], [1.236, 0.0, 0.0], [1.24, 0.0, 0.0], [1.244, 0.0, 0.0], [1.248, 0.0, 0.0], [1.252, 0.0, 0.0], [1.256, 0.0, 0.0], [1.26, 0.0, 0.0], [1.264, 0.0, 0.0], [1.268, 0.0, 0.0], [1.272, 0.0, 0.0], [1.276, 0.0, 0.0], [1.28, 0.0, 0.0], [1.284, 0.0, 0.0], [1.288, 0.0, 0.0], [1.292, 0.0, 0.0], [1.296, 0.0, 0.0], [1.3, 0.0, 0.0], [1.304, 0.0, 0.0], [1.308, 0.0, 0.0], [1.312, 0.0, 0.0], [1.316, 0.0, 0.0], [1.32, 0.0, 0.0], [1.324, 0.0, 0.0], [1.328, 0.0, 0.0], [1.332, 0.0, 0.0], [1.336, 0.0, 0.0], [1.34, 0.0, 0.0], [1.344, 0.0, 0.0], [1.348, 0.0, 0.0], [1.352, 0.0, 0.0], [1.356, 0.0, 0.0], [1.36, 0.0, 0.0], [1.364, 0.0, 0.0], [1.368, 0.0, 0.0], [1.372, 0.0, 0.0], [1.376, 0.0, 0.0], [1.38, 0.0, 0.0], [1.384, 0.0, 0.0], [1.388, 0.0, 0.0], [1.392, 0.0, 0.0], [1.396, 0.0, 0.0], [1.4, 0.0, 0.0], [1.404, 0.0, 0.0], [1.408, 0.0, 0.0], [1.412, 0.0, 0.0], [1.416, 0.0, 0.0], [1.42, 0.0, 0.0], [1.424, 0.0, 0.0], [1.428, 0.0, 0.0], [1.432, 0.0, 0.0], [1.436, 0.0, 0.0], [1.44, 0.0, 0.0], [1.444, 0.0, 0.0], [1.448, 0.0, 0.0], [1.452, 0.0, 0.0], [1.456, 0.0, 0.0], [1.46, 0.0, 0.0], [1.464, 0.0, 0.0], [1.468, 0.0, 0.0], [1.472, 0.0, 0.0], [1.476, 0.0, 0.0], [1.48, 0.0, 0.0], [1.484, 0.0, 0.0], [1.488, 0.0, 0.0], [1.492, 0.0, 0.0], [1.496, 0.0, 0.0], [1.5, 0.0, 0.0], [1.504, 0.0, 0.0], [1.508, 0.0, 0.0], [1.512, 0.0, 0.0], [1.516, 0.0, 0.0], [1.52, 0.0, 0.0], [1.524, 0.0, 0.0], [1.528, 0.0, 0.0], [1.532, 0.0, 0.0], [1.536, 0.0, 0.0], [1.54, 0.0, 0.0], [1.544, 0.0, 0.0], [1.548, 0.0, 0.0], [1.552, 0.0, 0.0], [1.556, 0.0, 0.0], [1.56, 0.0, 0.0], [1.564, 0.0, 0.0], [1.568, 0.0, 0.0], [1.572, 0.0, 0.0], [1.576, 0.0, 0.0], [1.58, 0.0, 0.0], [1.584, 0.0, 0.0], [1.588, 0.0, 0.0], [1.592, 0.0, 0.0], [1.596, 0.0, 0.0], [1.6, 0.0, 0.0], [1.604, 0.0, 0.0], [1.608, 0.0, 0.0], [1.612, 0.0, 0.0], [1.616, 0.0, 0.0], [1.62, 0.0, 0.0], [1.624, 0.0, 0.0], [1.628, 0.0, 0.0], [1.632, 0.0, 0.0], [1.636, 0.0, 0.0], [1.64, 0.0, 0.0], [1.644, 0.0, 0.0], [1.648, 0.0, 0.0], [1.652, 0.0, 0.0], [1.656, 0.0, 0.0], [1.66, 0.0, 0.0], [1.664, 0.0, 0.0], [1.668, 0.0, 0.0], [1.672, 0.0, 0.0], [1.676, 0.0, 0.0], [1.68, 0.0, 0.0], [1.684, 0.0, 0.0], [1.688, 0.0, 0.0], [1.692, 0.0, 0.0], [1.696, 0.0, 0.0], [1.7, 0.0, 0.0], [1.704, 0.0, 0.0], [1.708, 0.0, 0.0], [1.712, 0.0, 0.0], [1.716, 0.0, 0.0], [1.72, 0.0, 0.0], [1.724, 0.0, 0.0], [1.728, 0.0, 0.0], [1.732, 0.0, 0.0], [1.736, 0.0, 0.0], [1.74, 0.0, 0.0], [1.744, 0.0, 0.0], [1.748, 0.0, 0.0], [1.752, 0.0, 0.0], [1.756, 0.0, 0.0], [1.76, 0.0, 0.0], [1.764, 0.0, 0.0], [1.768, 0.0, 0.0], [1.772, 0.0, 0.0], [1.776, 0.0, 0.0], [1.78, 0.0, 0.0], [1.784, 0.0, 0.0], [1.788, 0.0, 0.0], [1.792, 0.0, 0.0], [1.796, 0.0, 0.0], [1.8, 0.0, 0.0], [1.804, 0.0, 0.0], [1.808, 0.0, 0.0], [1.812, 0.0, 0.0], [1.816, 0.0, 0.0], [1.82, 0.0, 0.0], [1.824, 0.0, 0.0], [1.828, 0.0, 0.0], [1.832, 0.0, 0.0], [1.836, 0.0, 0.0], [1.84, 0.0, 0.0], [1.844, 0.0, 0.0], [1.848, 0.0, 0.0], [1.852, 0.0, 0.0], [1.856, 0.0, 0.0], [1.86, 0.0, 0.0], [1.864, 0.0, 0.0], [1.868, 0.0, 0.0], [1.872, 0.0, 0.0], [1.876, 0.0, 0.0], [1.88, 0.0, 0.0], [1.884, 0.0, 0.0], [1.888, 0.0, 0.0], [1.892, 0.0, 0.0], [1.896, 0.0, 0.0], [1.9, 0.0, 0.0], [1.904, 0.0, 0.0], [1.908, 0.0, 0.0], [1.912, 0.0, 0.0], [1.916, 0.0, 0.0], [1.92, 0.0, 0.0], [1.924, 0.0, 0.0], [1.928, 0.0, 0.0], [1.932, 0.0, 0.0], [1.936, 0.0, 0.0], [1.94, 0.0, 0.0], [1.944, 0.0, 0.0], [1.948, 0.0, 0.0], [1.952, 0.0, 0.0], [1.956, 0.0, 0.0], [1.96, 0.0, 0.0], [1.964, 0.0, 0.0], [1.968, 0.0, 0.0], [1.972, 0.0, 0.0], [1.976, 0.0, 0.0], [1.98, 0.0, 0.0], [1.984, 0.0, 0.0], [1.988, 0.0, 0.0], [1.992, 0.0, 0.0], [1.996, 0.0, 0.0], [2.0, 0.0, 0.0]]
ABAQUS_ELEM = "B31"
E_MODULUS = 210000000000.0
NU = 0.3
RHO = 7850.0
A_AREA = 0.00131
I11 = 3.234e-07
I22 = 2.08769e-06
J_TORSION = 2.60673e-08
PRESCRIBED_NODES_DOF_VALUES = [(0, 0, 0.0), (0, 1, 0.0), (0, 2, 0.0), (0, 3, 0.0), (0, 4, 0.0), (0, 5, 0.0)]
POINT_LOADS = [[0.5, 0.0, 0.0, 0.0, -500.0, 0.0, 0.0, 0.0, 0.0]]
DISTRIBUTED_LOADS = []
DISTRIBUTED_EQUIVALENT_NODAL = []
OUT_CSV_DIR = r"C:\\Users\\s1834431\\Programs\\fem_model\\post_processing\\validation_visualisers\\abaqus_results\\job_0008_n500"

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
jobName = "job_0008_n500_abaqus"
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

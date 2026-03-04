# -*- coding: utf-8 -*-
# Auto-generated Abaqus CAE script for job: job_0007_n100
# Run with project Python (abqpy): python this_file; abqpy saveAs() launches Abaqus.
# Outputs CSV to: C:\Users\s1834431\Programs\fem_model\post_processing\validation_visualisers\abaqus_results\job_0007_n100

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
COORDS = [[0.0, 0.0, 0.0], [0.02, 0.0, 0.0], [0.04, 0.0, 0.0], [0.06, 0.0, 0.0], [0.08, 0.0, 0.0], [0.1, 0.0, 0.0], [0.12, 0.0, 0.0], [0.14, 0.0, 0.0], [0.16, 0.0, 0.0], [0.18, 0.0, 0.0], [0.2, 0.0, 0.0], [0.22, 0.0, 0.0], [0.24, 0.0, 0.0], [0.26, 0.0, 0.0], [0.28, 0.0, 0.0], [0.3, 0.0, 0.0], [0.32, 0.0, 0.0], [0.34, 0.0, 0.0], [0.36, 0.0, 0.0], [0.38, 0.0, 0.0], [0.4, 0.0, 0.0], [0.42, 0.0, 0.0], [0.44, 0.0, 0.0], [0.46, 0.0, 0.0], [0.48, 0.0, 0.0], [0.5, 0.0, 0.0], [0.52, 0.0, 0.0], [0.54, 0.0, 0.0], [0.56, 0.0, 0.0], [0.58, 0.0, 0.0], [0.6, 0.0, 0.0], [0.62, 0.0, 0.0], [0.64, 0.0, 0.0], [0.66, 0.0, 0.0], [0.68, 0.0, 0.0], [0.7, 0.0, 0.0], [0.72, 0.0, 0.0], [0.74, 0.0, 0.0], [0.76, 0.0, 0.0], [0.78, 0.0, 0.0], [0.8, 0.0, 0.0], [0.82, 0.0, 0.0], [0.84, 0.0, 0.0], [0.86, 0.0, 0.0], [0.88, 0.0, 0.0], [0.9, 0.0, 0.0], [0.92, 0.0, 0.0], [0.94, 0.0, 0.0], [0.96, 0.0, 0.0], [0.98, 0.0, 0.0], [1.0, 0.0, 0.0], [1.02, 0.0, 0.0], [1.04, 0.0, 0.0], [1.06, 0.0, 0.0], [1.08, 0.0, 0.0], [1.1, 0.0, 0.0], [1.12, 0.0, 0.0], [1.14, 0.0, 0.0], [1.16, 0.0, 0.0], [1.18, 0.0, 0.0], [1.2, 0.0, 0.0], [1.22, 0.0, 0.0], [1.24, 0.0, 0.0], [1.26, 0.0, 0.0], [1.28, 0.0, 0.0], [1.3, 0.0, 0.0], [1.32, 0.0, 0.0], [1.34, 0.0, 0.0], [1.36, 0.0, 0.0], [1.38, 0.0, 0.0], [1.4, 0.0, 0.0], [1.42, 0.0, 0.0], [1.44, 0.0, 0.0], [1.46, 0.0, 0.0], [1.48, 0.0, 0.0], [1.5, 0.0, 0.0], [1.52, 0.0, 0.0], [1.54, 0.0, 0.0], [1.56, 0.0, 0.0], [1.58, 0.0, 0.0], [1.6, 0.0, 0.0], [1.62, 0.0, 0.0], [1.64, 0.0, 0.0], [1.66, 0.0, 0.0], [1.68, 0.0, 0.0], [1.7, 0.0, 0.0], [1.72, 0.0, 0.0], [1.74, 0.0, 0.0], [1.76, 0.0, 0.0], [1.78, 0.0, 0.0], [1.8, 0.0, 0.0], [1.82, 0.0, 0.0], [1.84, 0.0, 0.0], [1.86, 0.0, 0.0], [1.88, 0.0, 0.0], [1.9, 0.0, 0.0], [1.92, 0.0, 0.0], [1.94, 0.0, 0.0], [1.96, 0.0, 0.0], [1.98, 0.0, 0.0], [2.0, 0.0, 0.0]]
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
DISTRIBUTED_LOADS = [[0.0, 0.0, 0.0, 0.0, -0.0, 0.0, 0.0, 0.0, 0.0], [0.02, 0.0, 0.0, 0.0, -0.05, 0.0, 0.0, 0.0, 0.0], [0.04, 0.0, 0.0, 0.0, -0.2, 0.0, 0.0, 0.0, 0.0], [0.06, 0.0, 0.0, 0.0, -0.45, 0.0, 0.0, 0.0, 0.0], [0.08, 0.0, 0.0, 0.0, -0.8, 0.0, 0.0, 0.0, 0.0], [0.1, 0.0, 0.0, 0.0, -1.25, 0.0, 0.0, 0.0, 0.0], [0.12, 0.0, 0.0, 0.0, -1.8, 0.0, 0.0, 0.0, 0.0], [0.14, 0.0, 0.0, 0.0, -2.45, 0.0, 0.0, 0.0, 0.0], [0.16, 0.0, 0.0, 0.0, -3.2, 0.0, 0.0, 0.0, 0.0], [0.18, 0.0, 0.0, 0.0, -4.05, 0.0, 0.0, 0.0, 0.0], [0.2, 0.0, 0.0, 0.0, -5.0, 0.0, 0.0, 0.0, 0.0], [0.22, 0.0, 0.0, 0.0, -6.05, 0.0, 0.0, 0.0, 0.0], [0.24, 0.0, 0.0, 0.0, -7.2, 0.0, 0.0, 0.0, 0.0], [0.26, 0.0, 0.0, 0.0, -8.45, 0.0, 0.0, 0.0, 0.0], [0.28, 0.0, 0.0, 0.0, -9.8, 0.0, 0.0, 0.0, 0.0], [0.3, 0.0, 0.0, 0.0, -11.25, 0.0, 0.0, 0.0, 0.0], [0.32, 0.0, 0.0, 0.0, -12.8, 0.0, 0.0, 0.0, 0.0], [0.34, 0.0, 0.0, 0.0, -14.45, 0.0, 0.0, 0.0, 0.0], [0.36, 0.0, 0.0, 0.0, -16.2, 0.0, 0.0, 0.0, 0.0], [0.38, 0.0, 0.0, 0.0, -18.05, 0.0, 0.0, 0.0, 0.0], [0.4, 0.0, 0.0, 0.0, -20.0, 0.0, 0.0, 0.0, 0.0], [0.42, 0.0, 0.0, 0.0, -22.05, 0.0, 0.0, 0.0, 0.0], [0.44, 0.0, 0.0, 0.0, -24.2, 0.0, 0.0, 0.0, 0.0], [0.46, 0.0, 0.0, 0.0, -26.45, 0.0, 0.0, 0.0, 0.0], [0.48, 0.0, 0.0, 0.0, -28.8, 0.0, 0.0, 0.0, 0.0], [0.5, 0.0, 0.0, 0.0, -31.25, 0.0, 0.0, 0.0, 0.0], [0.52, 0.0, 0.0, 0.0, -33.8, 0.0, 0.0, 0.0, 0.0], [0.54, 0.0, 0.0, 0.0, -36.45, 0.0, 0.0, 0.0, 0.0], [0.56, 0.0, 0.0, 0.0, -39.2, 0.0, 0.0, 0.0, 0.0], [0.58, 0.0, 0.0, 0.0, -42.05, 0.0, 0.0, 0.0, 0.0], [0.6, 0.0, 0.0, 0.0, -45.0, 0.0, 0.0, 0.0, 0.0], [0.62, 0.0, 0.0, 0.0, -48.05, 0.0, 0.0, 0.0, 0.0], [0.64, 0.0, 0.0, 0.0, -51.2, 0.0, 0.0, 0.0, 0.0], [0.66, 0.0, 0.0, 0.0, -54.45, 0.0, 0.0, 0.0, 0.0], [0.68, 0.0, 0.0, 0.0, -57.8, 0.0, 0.0, 0.0, 0.0], [0.7, 0.0, 0.0, 0.0, -61.25, 0.0, 0.0, 0.0, 0.0], [0.72, 0.0, 0.0, 0.0, -64.8, 0.0, 0.0, 0.0, 0.0], [0.74, 0.0, 0.0, 0.0, -68.45, 0.0, 0.0, 0.0, 0.0], [0.76, 0.0, 0.0, 0.0, -72.2, 0.0, 0.0, 0.0, 0.0], [0.78, 0.0, 0.0, 0.0, -76.05, 0.0, 0.0, 0.0, 0.0], [0.8, 0.0, 0.0, 0.0, -80.0, 0.0, 0.0, 0.0, 0.0], [0.82, 0.0, 0.0, 0.0, -84.05, 0.0, 0.0, 0.0, 0.0], [0.84, 0.0, 0.0, 0.0, -88.2, 0.0, 0.0, 0.0, 0.0], [0.86, 0.0, 0.0, 0.0, -92.45, 0.0, 0.0, 0.0, 0.0], [0.88, 0.0, 0.0, 0.0, -96.8, 0.0, 0.0, 0.0, 0.0], [0.9, 0.0, 0.0, 0.0, -101.25, 0.0, 0.0, 0.0, 0.0], [0.92, 0.0, 0.0, 0.0, -105.8, 0.0, 0.0, 0.0, 0.0], [0.94, 0.0, 0.0, 0.0, -110.45, 0.0, 0.0, 0.0, 0.0], [0.96, 0.0, 0.0, 0.0, -115.2, 0.0, 0.0, 0.0, 0.0], [0.98, 0.0, 0.0, 0.0, -120.05, 0.0, 0.0, 0.0, 0.0], [1.0, 0.0, 0.0, 0.0, -125.0, 0.0, 0.0, 0.0, 0.0], [1.02, 0.0, 0.0, 0.0, -130.05, 0.0, 0.0, 0.0, 0.0], [1.04, 0.0, 0.0, 0.0, -135.2, 0.0, 0.0, 0.0, 0.0], [1.06, 0.0, 0.0, 0.0, -140.45, 0.0, 0.0, 0.0, 0.0], [1.08, 0.0, 0.0, 0.0, -145.8, 0.0, 0.0, 0.0, 0.0], [1.1, 0.0, 0.0, 0.0, -151.25, 0.0, 0.0, 0.0, 0.0], [1.12, 0.0, 0.0, 0.0, -156.8, 0.0, 0.0, 0.0, 0.0], [1.14, 0.0, 0.0, 0.0, -162.45, 0.0, 0.0, 0.0, 0.0], [1.16, 0.0, 0.0, 0.0, -168.2, 0.0, 0.0, 0.0, 0.0], [1.18, 0.0, 0.0, 0.0, -174.05, 0.0, 0.0, 0.0, 0.0], [1.2, 0.0, 0.0, 0.0, -180.0, 0.0, 0.0, 0.0, 0.0], [1.22, 0.0, 0.0, 0.0, -186.05, 0.0, 0.0, 0.0, 0.0], [1.24, 0.0, 0.0, 0.0, -192.2, 0.0, 0.0, 0.0, 0.0], [1.26, 0.0, 0.0, 0.0, -198.45, 0.0, 0.0, 0.0, 0.0], [1.28, 0.0, 0.0, 0.0, -204.8, 0.0, 0.0, 0.0, 0.0], [1.3, 0.0, 0.0, 0.0, -211.25, 0.0, 0.0, 0.0, 0.0], [1.32, 0.0, 0.0, 0.0, -217.8, 0.0, 0.0, 0.0, 0.0], [1.34, 0.0, 0.0, 0.0, -224.45, 0.0, 0.0, 0.0, 0.0], [1.36, 0.0, 0.0, 0.0, -231.2, 0.0, 0.0, 0.0, 0.0], [1.38, 0.0, 0.0, 0.0, -238.05, 0.0, 0.0, 0.0, 0.0], [1.4, 0.0, 0.0, 0.0, -245.0, 0.0, 0.0, 0.0, 0.0], [1.42, 0.0, 0.0, 0.0, -252.05, 0.0, 0.0, 0.0, 0.0], [1.44, 0.0, 0.0, 0.0, -259.2, 0.0, 0.0, 0.0, 0.0], [1.46, 0.0, 0.0, 0.0, -266.45, 0.0, 0.0, 0.0, 0.0], [1.48, 0.0, 0.0, 0.0, -273.8, 0.0, 0.0, 0.0, 0.0], [1.5, 0.0, 0.0, 0.0, -281.25, 0.0, 0.0, 0.0, 0.0], [1.52, 0.0, 0.0, 0.0, -288.8, 0.0, 0.0, 0.0, 0.0], [1.54, 0.0, 0.0, 0.0, -296.45, 0.0, 0.0, 0.0, 0.0], [1.56, 0.0, 0.0, 0.0, -304.2, 0.0, 0.0, 0.0, 0.0], [1.58, 0.0, 0.0, 0.0, -312.05, 0.0, 0.0, 0.0, 0.0], [1.6, 0.0, 0.0, 0.0, -320.0, 0.0, 0.0, 0.0, 0.0], [1.62, 0.0, 0.0, 0.0, -328.05, 0.0, 0.0, 0.0, 0.0], [1.64, 0.0, 0.0, 0.0, -336.2, 0.0, 0.0, 0.0, 0.0], [1.66, 0.0, 0.0, 0.0, -344.45, 0.0, 0.0, 0.0, 0.0], [1.68, 0.0, 0.0, 0.0, -352.8, 0.0, 0.0, 0.0, 0.0], [1.7, 0.0, 0.0, 0.0, -361.25, 0.0, 0.0, 0.0, 0.0], [1.72, 0.0, 0.0, 0.0, -369.8, 0.0, 0.0, 0.0, 0.0], [1.74, 0.0, 0.0, 0.0, -378.45, 0.0, 0.0, 0.0, 0.0], [1.76, 0.0, 0.0, 0.0, -387.2, 0.0, 0.0, 0.0, 0.0], [1.78, 0.0, 0.0, 0.0, -396.05, 0.0, 0.0, 0.0, 0.0], [1.8, 0.0, 0.0, 0.0, -405.0, 0.0, 0.0, 0.0, 0.0], [1.82, 0.0, 0.0, 0.0, -414.05, 0.0, 0.0, 0.0, 0.0], [1.84, 0.0, 0.0, 0.0, -423.2, 0.0, 0.0, 0.0, 0.0], [1.86, 0.0, 0.0, 0.0, -432.45, 0.0, 0.0, 0.0, 0.0], [1.88, 0.0, 0.0, 0.0, -441.8, 0.0, 0.0, 0.0, 0.0], [1.9, 0.0, 0.0, 0.0, -451.25, 0.0, 0.0, 0.0, 0.0], [1.92, 0.0, 0.0, 0.0, -460.8, 0.0, 0.0, 0.0, 0.0], [1.94, 0.0, 0.0, 0.0, -470.45, 0.0, 0.0, 0.0, 0.0], [1.96, 0.0, 0.0, 0.0, -480.2, 0.0, 0.0, 0.0, 0.0], [1.98, 0.0, 0.0, 0.0, -490.05, 0.0, 0.0, 0.0, 0.0], [2.0, 0.0, 0.0, 0.0, -500.0, 0.0, 0.0, 0.0, 0.0]]
DISTRIBUTED_EQUIVALENT_NODAL = [(1, 0.0, -6.666666666666668e-05, 0.0, 0.0, 0.0, 0.0), (2, 0.0, -0.0011333333333333334, 0.0, 0.0, 0.0, 0.0), (3, 0.0, -0.0041333333333333335, 0.0, 0.0, 0.0, 0.0), (4, 0.0, -0.009133333333333335, 0.0, 0.0, 0.0, 0.0), (5, 0.0, -0.016133333333333336, 0.0, 0.0, 0.0, 0.0), (6, 0.0, -0.025133333333333327, 0.0, 0.0, 0.0, 0.0), (7, 0.0, -0.03613333333333335, 0.0, 0.0, 0.0, 0.0), (8, 0.0, -0.04913333333333335, 0.0, 0.0, 0.0, 0.0), (9, 0.0, -0.0641333333333333, 0.0, 0.0, 0.0, 0.0), (10, 0.0, -0.08113333333333336, 0.0, 0.0, 0.0, 0.0), (11, 0.0, -0.10013333333333335, 0.0, 0.0, 0.0, 0.0), (12, 0.0, -0.12113333333333326, 0.0, 0.0, 0.0, 0.0), (13, 0.0, -0.14413333333333334, 0.0, 0.0, 0.0, 0.0), (14, 0.0, -0.16913333333333347, 0.0, 0.0, 0.0, 0.0), (15, 0.0, -0.19613333333333322, 0.0, 0.0, 0.0, 0.0), (16, 0.0, -0.22513333333333324, 0.0, 0.0, 0.0, 0.0), (17, 0.0, -0.25613333333333355, 0.0, 0.0, 0.0, 0.0), (18, 0.0, -0.2891333333333331, 0.0, 0.0, 0.0, 0.0), (19, 0.0, -0.32413333333333316, 0.0, 0.0, 0.0, 0.0), (20, 0.0, -0.36113333333333364, 0.0, 0.0, 0.0, 0.0), (21, 0.0, -0.40013333333333306, 0.0, 0.0, 0.0, 0.0), (22, 0.0, -0.4411333333333331, 0.0, 0.0, 0.0, 0.0), (23, 0.0, -0.4841333333333338, 0.0, 0.0, 0.0, 0.0), (24, 0.0, -0.529133333333333, 0.0, 0.0, 0.0, 0.0), (25, 0.0, -0.5761333333333329, 0.0, 0.0, 0.0, 0.0), (26, 0.0, -0.6251333333333339, 0.0, 0.0, 0.0, 0.0), (27, 0.0, -0.676133333333334, 0.0, 0.0, 0.0, 0.0), (28, 0.0, -0.7291333333333341, 0.0, 0.0, 0.0, 0.0), (29, 0.0, -0.7841333333333319, 0.0, 0.0, 0.0, 0.0), (30, 0.0, -0.8411333333333317, 0.0, 0.0, 0.0, 0.0), (31, 0.0, -0.9001333333333341, 0.0, 0.0, 0.0, 0.0), (32, 0.0, -0.9611333333333341, 0.0, 0.0, 0.0, 0.0), (33, 0.0, -1.0241333333333342, 0.0, 0.0, 0.0, 0.0), (34, 0.0, -1.0891333333333342, 0.0, 0.0, 0.0, 0.0), (35, 0.0, -1.1561333333333308, 0.0, 0.0, 0.0, 0.0), (36, 0.0, -1.2251333333333312, 0.0, 0.0, 0.0, 0.0), (37, 0.0, -1.2961333333333345, 0.0, 0.0, 0.0, 0.0), (38, 0.0, -1.3691333333333346, 0.0, 0.0, 0.0, 0.0), (39, 0.0, -1.4441333333333346, 0.0, 0.0, 0.0, 0.0), (40, 0.0, -1.5211333333333348, 0.0, 0.0, 0.0, 0.0), (41, 0.0, -1.6001333333333303, 0.0, 0.0, 0.0, 0.0), (42, 0.0, -1.6811333333333305, 0.0, 0.0, 0.0, 0.0), (43, 0.0, -1.7641333333333349, 0.0, 0.0, 0.0, 0.0), (44, 0.0, -1.8491333333333348, 0.0, 0.0, 0.0, 0.0), (45, 0.0, -1.936133333333335, 0.0, 0.0, 0.0, 0.0), (46, 0.0, -2.025133333333335, 0.0, 0.0, 0.0, 0.0), (47, 0.0, -2.1161333333333294, 0.0, 0.0, 0.0, 0.0), (48, 0.0, -2.2091333333333294, 0.0, 0.0, 0.0, 0.0), (49, 0.0, -2.3041333333333354, 0.0, 0.0, 0.0, 0.0), (50, 0.0, -2.401133333333336, 0.0, 0.0, 0.0, 0.0), (51, 0.0, -2.500133333333335, 0.0, 0.0, 0.0, 0.0), (52, 0.0, -2.601133333333336, 0.0, 0.0, 0.0, 0.0), (53, 0.0, -2.7041333333333357, 0.0, 0.0, 0.0, 0.0), (54, 0.0, -2.809133333333336, 0.0, 0.0, 0.0, 0.0), (55, 0.0, -2.9161333333333364, 0.0, 0.0, 0.0, 0.0), (56, 0.0, -3.0251333333333363, 0.0, 0.0, 0.0, 0.0), (57, 0.0, -3.1361333333333192, 0.0, 0.0, 0.0, 0.0), (58, 0.0, -3.249133333333318, 0.0, 0.0, 0.0, 0.0), (59, 0.0, -3.3641333333333368, 0.0, 0.0, 0.0, 0.0), (60, 0.0, -3.481133333333336, 0.0, 0.0, 0.0, 0.0), (61, 0.0, -3.6001333333333365, 0.0, 0.0, 0.0, 0.0), (62, 0.0, -3.7211333333333365, 0.0, 0.0, 0.0, 0.0), (63, 0.0, -3.8441333333333363, 0.0, 0.0, 0.0, 0.0), (64, 0.0, -3.9691333333333363, 0.0, 0.0, 0.0, 0.0), (65, 0.0, -4.096133333333338, 0.0, 0.0, 0.0, 0.0), (66, 0.0, -4.225133333333337, 0.0, 0.0, 0.0, 0.0), (67, 0.0, -4.356133333333337, 0.0, 0.0, 0.0, 0.0), (68, 0.0, -4.489133333333337, 0.0, 0.0, 0.0, 0.0), (69, 0.0, -4.624133333333312, 0.0, 0.0, 0.0, 0.0), (70, 0.0, -4.761133333333312, 0.0, 0.0, 0.0, 0.0), (71, 0.0, -4.900133333333338, 0.0, 0.0, 0.0, 0.0), (72, 0.0, -5.041133333333338, 0.0, 0.0, 0.0, 0.0), (73, 0.0, -5.184133333333339, 0.0, 0.0, 0.0, 0.0), (74, 0.0, -5.329133333333338, 0.0, 0.0, 0.0, 0.0), (75, 0.0, -5.476133333333339, 0.0, 0.0, 0.0, 0.0), (76, 0.0, -5.6251333333333395, 0.0, 0.0, 0.0, 0.0), (77, 0.0, -5.776133333333338, 0.0, 0.0, 0.0, 0.0), (78, 0.0, -5.929133333333338, 0.0, 0.0, 0.0, 0.0), (79, 0.0, -6.084133333333339, 0.0, 0.0, 0.0, 0.0), (80, 0.0, -6.241133333333339, 0.0, 0.0, 0.0, 0.0), (81, 0.0, -6.40013333333334, 0.0, 0.0, 0.0, 0.0), (82, 0.0, -6.561133333333302, 0.0, 0.0, 0.0, 0.0), (83, 0.0, -6.724133333333302, 0.0, 0.0, 0.0, 0.0), (84, 0.0, -6.88913333333334, 0.0, 0.0, 0.0, 0.0), (85, 0.0, -7.0561333333333405, 0.0, 0.0, 0.0, 0.0), (86, 0.0, -7.225133333333341, 0.0, 0.0, 0.0, 0.0), (87, 0.0, -7.39613333333334, 0.0, 0.0, 0.0, 0.0), (88, 0.0, -7.56913333333334, 0.0, 0.0, 0.0, 0.0), (89, 0.0, -7.744133333333339, 0.0, 0.0, 0.0, 0.0), (90, 0.0, -7.921133333333341, 0.0, 0.0, 0.0, 0.0), (91, 0.0, -8.10013333333334, 0.0, 0.0, 0.0, 0.0), (92, 0.0, -8.28113333333334, 0.0, 0.0, 0.0, 0.0), (93, 0.0, -8.464133333333342, 0.0, 0.0, 0.0, 0.0), (94, 0.0, -8.649133333333292, 0.0, 0.0, 0.0, 0.0), (95, 0.0, -8.836133333333294, 0.0, 0.0, 0.0, 0.0), (96, 0.0, -9.025133333333342, 0.0, 0.0, 0.0, 0.0), (97, 0.0, -9.216133333333342, 0.0, 0.0, 0.0, 0.0), (98, 0.0, -9.409133333333342, 0.0, 0.0, 0.0, 0.0), (99, 0.0, -9.60413333333334, 0.0, 0.0, 0.0, 0.0), (100, 0.0, -9.801133333333343, 0.0, 0.0, 0.0, 0.0), (101, 0.0, -4.970066666666672, 0.0, 0.0, 0.0, 0.0)]
OUT_CSV_DIR = r"C:\\Users\\s1834431\\Programs\\fem_model\\post_processing\\validation_visualisers\\abaqus_results\\job_0007_n100"

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
jobName = "job_0007_n100_abaqus"
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

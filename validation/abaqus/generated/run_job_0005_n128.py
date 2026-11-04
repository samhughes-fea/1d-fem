# -*- coding: utf-8 -*-
# Auto-generated Abaqus CAE script for job: job_0005_n128
# Run with project Python (abqpy): python this_file; abqpy saveAs() launches Abaqus.
# Outputs CSV to: C:\Users\s1834431\Programs\fem_model\post_processing\validation_visualisers\abaqus_results\job_0005_n128

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
COORDS = [[0.0, 0.0, 0.0], [0.015625, 0.0, 0.0], [0.03125, 0.0, 0.0], [0.046875, 0.0, 0.0], [0.0625, 0.0, 0.0], [0.078125, 0.0, 0.0], [0.09375, 0.0, 0.0], [0.109375, 0.0, 0.0], [0.125, 0.0, 0.0], [0.140625, 0.0, 0.0], [0.15625, 0.0, 0.0], [0.171875, 0.0, 0.0], [0.1875, 0.0, 0.0], [0.203125, 0.0, 0.0], [0.21875, 0.0, 0.0], [0.234375, 0.0, 0.0], [0.25, 0.0, 0.0], [0.265625, 0.0, 0.0], [0.28125, 0.0, 0.0], [0.296875, 0.0, 0.0], [0.3125, 0.0, 0.0], [0.328125, 0.0, 0.0], [0.34375, 0.0, 0.0], [0.359375, 0.0, 0.0], [0.375, 0.0, 0.0], [0.390625, 0.0, 0.0], [0.40625, 0.0, 0.0], [0.421875, 0.0, 0.0], [0.4375, 0.0, 0.0], [0.453125, 0.0, 0.0], [0.46875, 0.0, 0.0], [0.484375, 0.0, 0.0], [0.5, 0.0, 0.0], [0.515625, 0.0, 0.0], [0.53125, 0.0, 0.0], [0.546875, 0.0, 0.0], [0.5625, 0.0, 0.0], [0.578125, 0.0, 0.0], [0.59375, 0.0, 0.0], [0.609375, 0.0, 0.0], [0.625, 0.0, 0.0], [0.640625, 0.0, 0.0], [0.65625, 0.0, 0.0], [0.671875, 0.0, 0.0], [0.6875, 0.0, 0.0], [0.703125, 0.0, 0.0], [0.71875, 0.0, 0.0], [0.734375, 0.0, 0.0], [0.75, 0.0, 0.0], [0.765625, 0.0, 0.0], [0.78125, 0.0, 0.0], [0.796875, 0.0, 0.0], [0.8125, 0.0, 0.0], [0.828125, 0.0, 0.0], [0.84375, 0.0, 0.0], [0.859375, 0.0, 0.0], [0.875, 0.0, 0.0], [0.890625, 0.0, 0.0], [0.90625, 0.0, 0.0], [0.921875, 0.0, 0.0], [0.9375, 0.0, 0.0], [0.953125, 0.0, 0.0], [0.96875, 0.0, 0.0], [0.984375, 0.0, 0.0], [1.0, 0.0, 0.0], [1.015625, 0.0, 0.0], [1.03125, 0.0, 0.0], [1.046875, 0.0, 0.0], [1.0625, 0.0, 0.0], [1.078125, 0.0, 0.0], [1.09375, 0.0, 0.0], [1.109375, 0.0, 0.0], [1.125, 0.0, 0.0], [1.140625, 0.0, 0.0], [1.15625, 0.0, 0.0], [1.171875, 0.0, 0.0], [1.1875, 0.0, 0.0], [1.203125, 0.0, 0.0], [1.21875, 0.0, 0.0], [1.234375, 0.0, 0.0], [1.25, 0.0, 0.0], [1.265625, 0.0, 0.0], [1.28125, 0.0, 0.0], [1.296875, 0.0, 0.0], [1.3125, 0.0, 0.0], [1.328125, 0.0, 0.0], [1.34375, 0.0, 0.0], [1.359375, 0.0, 0.0], [1.375, 0.0, 0.0], [1.390625, 0.0, 0.0], [1.40625, 0.0, 0.0], [1.421875, 0.0, 0.0], [1.4375, 0.0, 0.0], [1.453125, 0.0, 0.0], [1.46875, 0.0, 0.0], [1.484375, 0.0, 0.0], [1.5, 0.0, 0.0], [1.515625, 0.0, 0.0], [1.53125, 0.0, 0.0], [1.546875, 0.0, 0.0], [1.5625, 0.0, 0.0], [1.578125, 0.0, 0.0], [1.59375, 0.0, 0.0], [1.609375, 0.0, 0.0], [1.625, 0.0, 0.0], [1.640625, 0.0, 0.0], [1.65625, 0.0, 0.0], [1.671875, 0.0, 0.0], [1.6875, 0.0, 0.0], [1.703125, 0.0, 0.0], [1.71875, 0.0, 0.0], [1.734375, 0.0, 0.0], [1.75, 0.0, 0.0], [1.765625, 0.0, 0.0], [1.78125, 0.0, 0.0], [1.796875, 0.0, 0.0], [1.8125, 0.0, 0.0], [1.828125, 0.0, 0.0], [1.84375, 0.0, 0.0], [1.859375, 0.0, 0.0], [1.875, 0.0, 0.0], [1.890625, 0.0, 0.0], [1.90625, 0.0, 0.0], [1.921875, 0.0, 0.0], [1.9375, 0.0, 0.0], [1.953125, 0.0, 0.0], [1.96875, 0.0, 0.0], [1.984375, 0.0, 0.0], [2.0, 0.0, 0.0]]
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
DISTRIBUTED_LOADS = [[0.0, 0.0, 0.0, 0.0, -0.0, 0.0, 0.0, 0.0, 0.0], [0.015625, 0.0, 0.0, 0.0, -0.030518, 0.0, 0.0, 0.0, 0.0], [0.03125, 0.0, 0.0, 0.0, -0.12207, 0.0, 0.0, 0.0, 0.0], [0.046875, 0.0, 0.0, 0.0, -0.274658, 0.0, 0.0, 0.0, 0.0], [0.0625, 0.0, 0.0, 0.0, -0.488281, 0.0, 0.0, 0.0, 0.0], [0.078125, 0.0, 0.0, 0.0, -0.762939, 0.0, 0.0, 0.0, 0.0], [0.09375, 0.0, 0.0, 0.0, -1.098633, 0.0, 0.0, 0.0, 0.0], [0.109375, 0.0, 0.0, 0.0, -1.495361, 0.0, 0.0, 0.0, 0.0], [0.125, 0.0, 0.0, 0.0, -1.953125, 0.0, 0.0, 0.0, 0.0], [0.140625, 0.0, 0.0, 0.0, -2.471924, 0.0, 0.0, 0.0, 0.0], [0.15625, 0.0, 0.0, 0.0, -3.051758, 0.0, 0.0, 0.0, 0.0], [0.171875, 0.0, 0.0, 0.0, -3.692627, 0.0, 0.0, 0.0, 0.0], [0.1875, 0.0, 0.0, 0.0, -4.394531, 0.0, 0.0, 0.0, 0.0], [0.203125, 0.0, 0.0, 0.0, -5.157471, 0.0, 0.0, 0.0, 0.0], [0.21875, 0.0, 0.0, 0.0, -5.981445, 0.0, 0.0, 0.0, 0.0], [0.234375, 0.0, 0.0, 0.0, -6.866455, 0.0, 0.0, 0.0, 0.0], [0.25, 0.0, 0.0, 0.0, -7.8125, 0.0, 0.0, 0.0, 0.0], [0.265625, 0.0, 0.0, 0.0, -8.81958, 0.0, 0.0, 0.0, 0.0], [0.28125, 0.0, 0.0, 0.0, -9.887695, 0.0, 0.0, 0.0, 0.0], [0.296875, 0.0, 0.0, 0.0, -11.016846, 0.0, 0.0, 0.0, 0.0], [0.3125, 0.0, 0.0, 0.0, -12.207031, 0.0, 0.0, 0.0, 0.0], [0.328125, 0.0, 0.0, 0.0, -13.458252, 0.0, 0.0, 0.0, 0.0], [0.34375, 0.0, 0.0, 0.0, -14.770508, 0.0, 0.0, 0.0, 0.0], [0.359375, 0.0, 0.0, 0.0, -16.143799, 0.0, 0.0, 0.0, 0.0], [0.375, 0.0, 0.0, 0.0, -17.578125, 0.0, 0.0, 0.0, 0.0], [0.390625, 0.0, 0.0, 0.0, -19.073486, 0.0, 0.0, 0.0, 0.0], [0.40625, 0.0, 0.0, 0.0, -20.629883, 0.0, 0.0, 0.0, 0.0], [0.421875, 0.0, 0.0, 0.0, -22.247314, 0.0, 0.0, 0.0, 0.0], [0.4375, 0.0, 0.0, 0.0, -23.925781, 0.0, 0.0, 0.0, 0.0], [0.453125, 0.0, 0.0, 0.0, -25.665283, 0.0, 0.0, 0.0, 0.0], [0.46875, 0.0, 0.0, 0.0, -27.46582, 0.0, 0.0, 0.0, 0.0], [0.484375, 0.0, 0.0, 0.0, -29.327393, 0.0, 0.0, 0.0, 0.0], [0.5, 0.0, 0.0, 0.0, -31.25, 0.0, 0.0, 0.0, 0.0], [0.515625, 0.0, 0.0, 0.0, -33.233643, 0.0, 0.0, 0.0, 0.0], [0.53125, 0.0, 0.0, 0.0, -35.27832, 0.0, 0.0, 0.0, 0.0], [0.546875, 0.0, 0.0, 0.0, -37.384033, 0.0, 0.0, 0.0, 0.0], [0.5625, 0.0, 0.0, 0.0, -39.550781, 0.0, 0.0, 0.0, 0.0], [0.578125, 0.0, 0.0, 0.0, -41.778564, 0.0, 0.0, 0.0, 0.0], [0.59375, 0.0, 0.0, 0.0, -44.067383, 0.0, 0.0, 0.0, 0.0], [0.609375, 0.0, 0.0, 0.0, -46.417236, 0.0, 0.0, 0.0, 0.0], [0.625, 0.0, 0.0, 0.0, -48.828125, 0.0, 0.0, 0.0, 0.0], [0.640625, 0.0, 0.0, 0.0, -51.300049, 0.0, 0.0, 0.0, 0.0], [0.65625, 0.0, 0.0, 0.0, -53.833008, 0.0, 0.0, 0.0, 0.0], [0.671875, 0.0, 0.0, 0.0, -56.427002, 0.0, 0.0, 0.0, 0.0], [0.6875, 0.0, 0.0, 0.0, -59.082031, 0.0, 0.0, 0.0, 0.0], [0.703125, 0.0, 0.0, 0.0, -61.798096, 0.0, 0.0, 0.0, 0.0], [0.71875, 0.0, 0.0, 0.0, -64.575195, 0.0, 0.0, 0.0, 0.0], [0.734375, 0.0, 0.0, 0.0, -67.41333, 0.0, 0.0, 0.0, 0.0], [0.75, 0.0, 0.0, 0.0, -70.3125, 0.0, 0.0, 0.0, 0.0], [0.765625, 0.0, 0.0, 0.0, -73.272705, 0.0, 0.0, 0.0, 0.0], [0.78125, 0.0, 0.0, 0.0, -76.293945, 0.0, 0.0, 0.0, 0.0], [0.796875, 0.0, 0.0, 0.0, -79.376221, 0.0, 0.0, 0.0, 0.0], [0.8125, 0.0, 0.0, 0.0, -82.519531, 0.0, 0.0, 0.0, 0.0], [0.828125, 0.0, 0.0, 0.0, -85.723877, 0.0, 0.0, 0.0, 0.0], [0.84375, 0.0, 0.0, 0.0, -88.989258, 0.0, 0.0, 0.0, 0.0], [0.859375, 0.0, 0.0, 0.0, -92.315674, 0.0, 0.0, 0.0, 0.0], [0.875, 0.0, 0.0, 0.0, -95.703125, 0.0, 0.0, 0.0, 0.0], [0.890625, 0.0, 0.0, 0.0, -99.151611, 0.0, 0.0, 0.0, 0.0], [0.90625, 0.0, 0.0, 0.0, -102.661133, 0.0, 0.0, 0.0, 0.0], [0.921875, 0.0, 0.0, 0.0, -106.231689, 0.0, 0.0, 0.0, 0.0], [0.9375, 0.0, 0.0, 0.0, -109.863281, 0.0, 0.0, 0.0, 0.0], [0.953125, 0.0, 0.0, 0.0, -113.555908, 0.0, 0.0, 0.0, 0.0], [0.96875, 0.0, 0.0, 0.0, -117.30957, 0.0, 0.0, 0.0, 0.0], [0.984375, 0.0, 0.0, 0.0, -121.124268, 0.0, 0.0, 0.0, 0.0], [1.0, 0.0, 0.0, 0.0, -125.0, 0.0, 0.0, 0.0, 0.0], [1.015625, 0.0, 0.0, 0.0, -128.936768, 0.0, 0.0, 0.0, 0.0], [1.03125, 0.0, 0.0, 0.0, -132.93457, 0.0, 0.0, 0.0, 0.0], [1.046875, 0.0, 0.0, 0.0, -136.993408, 0.0, 0.0, 0.0, 0.0], [1.0625, 0.0, 0.0, 0.0, -141.113281, 0.0, 0.0, 0.0, 0.0], [1.078125, 0.0, 0.0, 0.0, -145.294189, 0.0, 0.0, 0.0, 0.0], [1.09375, 0.0, 0.0, 0.0, -149.536133, 0.0, 0.0, 0.0, 0.0], [1.109375, 0.0, 0.0, 0.0, -153.839111, 0.0, 0.0, 0.0, 0.0], [1.125, 0.0, 0.0, 0.0, -158.203125, 0.0, 0.0, 0.0, 0.0], [1.140625, 0.0, 0.0, 0.0, -162.628174, 0.0, 0.0, 0.0, 0.0], [1.15625, 0.0, 0.0, 0.0, -167.114258, 0.0, 0.0, 0.0, 0.0], [1.171875, 0.0, 0.0, 0.0, -171.661377, 0.0, 0.0, 0.0, 0.0], [1.1875, 0.0, 0.0, 0.0, -176.269531, 0.0, 0.0, 0.0, 0.0], [1.203125, 0.0, 0.0, 0.0, -180.938721, 0.0, 0.0, 0.0, 0.0], [1.21875, 0.0, 0.0, 0.0, -185.668945, 0.0, 0.0, 0.0, 0.0], [1.234375, 0.0, 0.0, 0.0, -190.460205, 0.0, 0.0, 0.0, 0.0], [1.25, 0.0, 0.0, 0.0, -195.3125, 0.0, 0.0, 0.0, 0.0], [1.265625, 0.0, 0.0, 0.0, -200.22583, 0.0, 0.0, 0.0, 0.0], [1.28125, 0.0, 0.0, 0.0, -205.200195, 0.0, 0.0, 0.0, 0.0], [1.296875, 0.0, 0.0, 0.0, -210.235596, 0.0, 0.0, 0.0, 0.0], [1.3125, 0.0, 0.0, 0.0, -215.332031, 0.0, 0.0, 0.0, 0.0], [1.328125, 0.0, 0.0, 0.0, -220.489502, 0.0, 0.0, 0.0, 0.0], [1.34375, 0.0, 0.0, 0.0, -225.708008, 0.0, 0.0, 0.0, 0.0], [1.359375, 0.0, 0.0, 0.0, -230.987549, 0.0, 0.0, 0.0, 0.0], [1.375, 0.0, 0.0, 0.0, -236.328125, 0.0, 0.0, 0.0, 0.0], [1.390625, 0.0, 0.0, 0.0, -241.729736, 0.0, 0.0, 0.0, 0.0], [1.40625, 0.0, 0.0, 0.0, -247.192383, 0.0, 0.0, 0.0, 0.0], [1.421875, 0.0, 0.0, 0.0, -252.716064, 0.0, 0.0, 0.0, 0.0], [1.4375, 0.0, 0.0, 0.0, -258.300781, 0.0, 0.0, 0.0, 0.0], [1.453125, 0.0, 0.0, 0.0, -263.946533, 0.0, 0.0, 0.0, 0.0], [1.46875, 0.0, 0.0, 0.0, -269.65332, 0.0, 0.0, 0.0, 0.0], [1.484375, 0.0, 0.0, 0.0, -275.421143, 0.0, 0.0, 0.0, 0.0], [1.5, 0.0, 0.0, 0.0, -281.25, 0.0, 0.0, 0.0, 0.0], [1.515625, 0.0, 0.0, 0.0, -287.139893, 0.0, 0.0, 0.0, 0.0], [1.53125, 0.0, 0.0, 0.0, -293.09082, 0.0, 0.0, 0.0, 0.0], [1.546875, 0.0, 0.0, 0.0, -299.102783, 0.0, 0.0, 0.0, 0.0], [1.5625, 0.0, 0.0, 0.0, -305.175781, 0.0, 0.0, 0.0, 0.0], [1.578125, 0.0, 0.0, 0.0, -311.309814, 0.0, 0.0, 0.0, 0.0], [1.59375, 0.0, 0.0, 0.0, -317.504883, 0.0, 0.0, 0.0, 0.0], [1.609375, 0.0, 0.0, 0.0, -323.760986, 0.0, 0.0, 0.0, 0.0], [1.625, 0.0, 0.0, 0.0, -330.078125, 0.0, 0.0, 0.0, 0.0], [1.640625, 0.0, 0.0, 0.0, -336.456299, 0.0, 0.0, 0.0, 0.0], [1.65625, 0.0, 0.0, 0.0, -342.895508, 0.0, 0.0, 0.0, 0.0], [1.671875, 0.0, 0.0, 0.0, -349.395752, 0.0, 0.0, 0.0, 0.0], [1.6875, 0.0, 0.0, 0.0, -355.957031, 0.0, 0.0, 0.0, 0.0], [1.703125, 0.0, 0.0, 0.0, -362.579346, 0.0, 0.0, 0.0, 0.0], [1.71875, 0.0, 0.0, 0.0, -369.262695, 0.0, 0.0, 0.0, 0.0], [1.734375, 0.0, 0.0, 0.0, -376.00708, 0.0, 0.0, 0.0, 0.0], [1.75, 0.0, 0.0, 0.0, -382.8125, 0.0, 0.0, 0.0, 0.0], [1.765625, 0.0, 0.0, 0.0, -389.678955, 0.0, 0.0, 0.0, 0.0], [1.78125, 0.0, 0.0, 0.0, -396.606445, 0.0, 0.0, 0.0, 0.0], [1.796875, 0.0, 0.0, 0.0, -403.594971, 0.0, 0.0, 0.0, 0.0], [1.8125, 0.0, 0.0, 0.0, -410.644531, 0.0, 0.0, 0.0, 0.0], [1.828125, 0.0, 0.0, 0.0, -417.755127, 0.0, 0.0, 0.0, 0.0], [1.84375, 0.0, 0.0, 0.0, -424.926758, 0.0, 0.0, 0.0, 0.0], [1.859375, 0.0, 0.0, 0.0, -432.159424, 0.0, 0.0, 0.0, 0.0], [1.875, 0.0, 0.0, 0.0, -439.453125, 0.0, 0.0, 0.0, 0.0], [1.890625, 0.0, 0.0, 0.0, -446.807861, 0.0, 0.0, 0.0, 0.0], [1.90625, 0.0, 0.0, 0.0, -454.223633, 0.0, 0.0, 0.0, 0.0], [1.921875, 0.0, 0.0, 0.0, -461.700439, 0.0, 0.0, 0.0, 0.0], [1.9375, 0.0, 0.0, 0.0, -469.238281, 0.0, 0.0, 0.0, 0.0], [1.953125, 0.0, 0.0, 0.0, -476.837158, 0.0, 0.0, 0.0, 0.0], [1.96875, 0.0, 0.0, 0.0, -484.49707, 0.0, 0.0, 0.0, 0.0], [1.984375, 0.0, 0.0, 0.0, -492.218018, 0.0, 0.0, 0.0, 0.0], [2.0, 0.0, 0.0, 0.0, -500.0, 0.0, 0.0, 0.0, 0.0]]
DISTRIBUTED_EQUIVALENT_NODAL = [(1, 0.0, -3.1791881819292867e-05, 0.0, 0.0, 0.0, 0.0), (2, 0.0, -0.0005404208333333333, 0.0, 0.0, 0.0, 0.0), (3, 0.0, -0.0019709233294277705, 0.0, 0.0, 0.0, 0.0), (4, 0.0, -0.004355109208330586, 0.0, 0.0, 0.0, 0.0), (5, 0.0, -0.007692968509124888, 0.0, 0.0, 0.0, 0.0), (6, 0.0, -0.01198450167704487, 0.0, 0.0, 0.0, 0.0), (7, 0.0, -0.017229716892070646, 0.0, 0.0, 0.0, 0.0), (8, 0.0, -0.023428595442172566, 0.0, 0.0, 0.0, 0.0), (9, 0.0, -0.030581155948614092, 0.0, 0.0, 0.0, 0.0), (10, 0.0, -0.03868739068524607, 0.0, 0.0, 0.0, 0.0), (11, 0.0, -0.04774729693540165, 0.0, 0.0, 0.0, 0.0), (12, 0.0, -0.057760874698147384, 0.0, 0.0, 0.0, 0.0), (13, 0.0, -0.06872812669388387, 0.0, 0.0, 0.0, 0.0), (14, 0.0, -0.08064906063569222, 0.0, 0.0, 0.0, 0.0), (15, 0.0, -0.09352365795084729, 0.0, 0.0, 0.0, 0.0), (16, 0.0, -0.10735193717029362, 0.0, 0.0, 0.0, 0.0), (17, 0.0, -0.12213389078985319, 0.0, 0.0, 0.0, 0.0), (18, 0.0, -0.13786951529529365, 0.0, 0.0, 0.0, 0.0), (19, 0.0, -0.15455881420084736, 0.0, 0.0, 0.0, 0.0), (20, 0.0, -0.17220179501069213, 0.0, 0.0, 0.0, 0.0), (21, 0.0, -0.19079843919388434, 0.0, 0.0, 0.0, 0.0), (22, 0.0, -0.21034876532314573, 0.0, 0.0, 0.0, 0.0), (23, 0.0, -0.23085276568540802, 0.0, 0.0, 0.0, 0.0), (24, 0.0, -0.2523104375602224, 0.0, 0.0, 0.0, 0.0), (25, 0.0, -0.2747217809487028, 0.0, 0.0, 0.0, 0.0), (26, 0.0, -0.29808679856684206, 0.0, 0.0, 0.0, 0.0), (27, 0.0, -0.3224054981433043, 0.0, 0.0, 0.0, 0.0), (28, 0.0, -0.3476778610474409, 0.0, 0.0, 0.0, 0.0), (29, 0.0, -0.37390390602630735, 0.0, 0.0, 0.0, 0.0), (30, 0.0, -0.401083624769205, 0.0, 0.0, 0.0, 0.0), (31, 0.0, -0.42921701731874773, 0.0, 0.0, 0.0, 0.0), (32, 0.0, -0.45830409181517934, 0.0, 0.0, 0.0, 0.0), (33, 0.0, -0.4883448301080354, 0.0, 0.0, 0.0, 0.0), (34, 0.0, -0.5193392480651793, 0.0, 0.0, 0.0, 0.0), (35, 0.0, -0.5512873298187477, 0.0, 0.0, 0.0, 0.0), (36, 0.0, -0.5841890935192051, 0.0, 0.0, 0.0, 0.0), (37, 0.0, -0.6180445310263074, 0.0, 0.0, 0.0, 0.0), (38, 0.0, -0.652853642297441, 0.0, 0.0, 0.0, 0.0), (39, 0.0, -0.6886164356433043, 0.0, 0.0, 0.0, 0.0), (40, 0.0, -0.7253328923168421, 0.0, 0.0, 0.0, 0.0), (41, 0.0, -0.7630030309487028, 0.0, 0.0, 0.0, 0.0), (42, 0.0, -0.8016268438102224, 0.0, 0.0, 0.0, 0.0), (43, 0.0, -0.8412043281854082, 0.0, 0.0, 0.0, 0.0), (44, 0.0, -0.8817354840731457, 0.0, 0.0, 0.0, 0.0), (45, 0.0, -0.9232203141938844, 0.0, 0.0, 0.0, 0.0), (46, 0.0, -0.9656588262606922, 0.0, 0.0, 0.0, 0.0), (47, 0.0, -1.0090510017008474, 0.0, 0.0, 0.0, 0.0), (48, 0.0, -1.0533968590452938, 0.0, 0.0, 0.0, 0.0), (49, 0.0, -1.0986963907898535, 0.0, 0.0, 0.0, 0.0), (50, 0.0, -1.1449495934202938, 0.0, 0.0, 0.0, 0.0), (51, 0.0, -1.1921564704508474, 0.0, 0.0, 0.0, 0.0), (52, 0.0, -1.2403170293856924, 0.0, 0.0, 0.0, 0.0), (53, 0.0, -1.2894312516938844, 0.0, 0.0, 0.0, 0.0), (54, 0.0, -1.339499155948146, 0.0, 0.0, 0.0, 0.0), (55, 0.0, -1.390520734435408, 0.0, 0.0, 0.0, 0.0), (56, 0.0, -1.4424959844352228, 0.0, 0.0, 0.0, 0.0), (57, 0.0, -1.495424905948703, 0.0, 0.0, 0.0, 0.0), (58, 0.0, -1.5493075016918423, 0.0, 0.0, 0.0, 0.0), (59, 0.0, -1.6041437793933047, 0.0, 0.0, 0.0, 0.0), (60, 0.0, -1.6599337204224411, 0.0, 0.0, 0.0, 0.0), (61, 0.0, -1.7166773435263079, 0.0, 0.0, 0.0, 0.0), (62, 0.0, -1.7743746403942053, 0.0, 0.0, 0.0, 0.0), (63, 0.0, -1.8330256110687477, 0.0, 0.0, 0.0, 0.0), (64, 0.0, -1.8926302636901795, 0.0, 0.0, 0.0, 0.0), (65, 0.0, -1.9531885801080358, 0.0, 0.0, 0.0, 0.0), (66, 0.0, -2.0147005761901795, 0.0, 0.0, 0.0, 0.0), (67, 0.0, -2.0771662360687486, 0.0, 0.0, 0.0, 0.0), (68, 0.0, -2.1405855778942056, 0.0, 0.0, 0.0, 0.0), (69, 0.0, -2.2049585935263085, 0.0, 0.0, 0.0, 0.0), (70, 0.0, -2.270285282922441, 0.0, 0.0, 0.0, 0.0), (71, 0.0, -2.336565654393304, 0.0, 0.0, 0.0, 0.0), (72, 0.0, -2.4037996891918425, 0.0, 0.0, 0.0, 0.0), (73, 0.0, -2.4719874059487026, 0.0, 0.0, 0.0, 0.0), (74, 0.0, -2.5411287969352228, 0.0, 0.0, 0.0, 0.0), (75, 0.0, -2.611223859435409, 0.0, 0.0, 0.0, 0.0), (76, 0.0, -2.6822725934481455, 0.0, 0.0, 0.0, 0.0), (77, 0.0, -2.754275001693885, 0.0, 0.0, 0.0, 0.0), (78, 0.0, -2.827231091885692, 0.0, 0.0, 0.0, 0.0), (79, 0.0, -2.901140845450848, 0.0, 0.0, 0.0, 0.0), (80, 0.0, -2.9760042809202947, 0.0, 0.0, 0.0, 0.0), (81, 0.0, -3.051821390789854, 0.0, 0.0, 0.0, 0.0), (82, 0.0, -3.1285921715452947, 0.0, 0.0, 0.0, 0.0), (83, 0.0, -3.2063166267008483, 0.0, 0.0, 0.0, 0.0), (84, 0.0, -3.284994763760693, 0.0, 0.0, 0.0, 0.0), (85, 0.0, -3.364626564193885, 0.0, 0.0, 0.0, 0.0), (86, 0.0, -3.445212046573146, 0.0, 0.0, 0.0, 0.0), (87, 0.0, -3.526751203185409, 0.0, 0.0, 0.0, 0.0), (88, 0.0, -3.6092440313102228, 0.0, 0.0, 0.0, 0.0), (89, 0.0, -3.6926905309487035, 0.0, 0.0, 0.0, 0.0), (90, 0.0, -3.7770907048168425, 0.0, 0.0, 0.0, 0.0), (91, 0.0, -3.862444560643305, 0.0, 0.0, 0.0, 0.0), (92, 0.0, -3.948752079797442, 0.0, 0.0, 0.0, 0.0), (93, 0.0, -4.0360132810263085, 0.0, 0.0, 0.0, 0.0), (94, 0.0, -4.124228156019206, 0.0, 0.0, 0.0, 0.0), (95, 0.0, -4.213396704818749, 0.0, 0.0, 0.0, 0.0), (96, 0.0, -4.3035189355651795, 0.0, 0.0, 0.0, 0.0), (97, 0.0, -4.394594830108037, 0.0, 0.0, 0.0, 0.0), (98, 0.0, -4.4866244043151795, 0.0, 0.0, 0.0, 0.0), (99, 0.0, -4.579607642318749, 0.0, 0.0, 0.0, 0.0), (100, 0.0, -4.673544562269205, 0.0, 0.0, 0.0, 0.0), (101, 0.0, -4.7684351560263085, 0.0, 0.0, 0.0, 0.0), (102, 0.0, -4.864279423547442, 0.0, 0.0, 0.0, 0.0), (103, 0.0, -4.961077373143306, 0.0, 0.0, 0.0, 0.0), (104, 0.0, -5.058828986066843, 0.0, 0.0, 0.0, 0.0), (105, 0.0, -5.157534280948703, 0.0, 0.0, 0.0, 0.0), (106, 0.0, -5.257193250060223, 0.0, 0.0, 0.0, 0.0), (107, 0.0, -5.357805890685409, 0.0, 0.0, 0.0, 0.0), (108, 0.0, -5.459372202823147, 0.0, 0.0, 0.0, 0.0), (109, 0.0, -5.561892189193884, 0.0, 0.0, 0.0, 0.0), (110, 0.0, -5.665365857510693, 0.0, 0.0, 0.0, 0.0), (111, 0.0, -5.769793189200849, 0.0, 0.0, 0.0, 0.0), (112, 0.0, -5.875174202795295, 0.0, 0.0, 0.0, 0.0), (113, 0.0, -5.981508890789854, 0.0, 0.0, 0.0, 0.0), (114, 0.0, -6.088797249670295, 0.0, 0.0, 0.0, 0.0), (115, 0.0, -6.19703928295085, 0.0, 0.0, 0.0, 0.0), (116, 0.0, -6.306234998135694, 0.0, 0.0, 0.0, 0.0), (117, 0.0, -6.416384376693886, 0.0, 0.0, 0.0, 0.0), (118, 0.0, -6.52748743719815, 0.0, 0.0, 0.0, 0.0), (119, 0.0, -6.639544171935403, 0.0, 0.0, 0.0, 0.0), (120, 0.0, -6.752554578185247, 0.0, 0.0, 0.0, 0.0), (121, 0.0, -6.8665186559486155, 0.0, 0.0, 0.0, 0.0), (122, 0.0, -6.981436407942174, 0.0, 0.0, 0.0, 0.0), (123, 0.0, -7.097307841892071, 0.0, 0.0, 0.0, 0.0), (124, 0.0, -7.214132939177047, 0.0, 0.0, 0.0, 0.0), (125, 0.0, -7.331911718509126, 0.0, 0.0, 0.0, 0.0), (126, 0.0, -7.450644171708332, 0.0, 0.0, 0.0, 0.0), (127, 0.0, -7.570330298329427, 0.0, 0.0, 0.0, 0.0), (128, 0.0, -7.690970108333334, 0.0, 0.0, 0.0, 0.0), (129, 0.0, -3.8879712450068196, 0.0, 0.0, 0.0, 0.0)]
OUT_CSV_DIR = r"C:\\Users\\s1834431\\Programs\\fem_model\\post_processing\\validation_visualisers\\abaqus_results\\job_0005_n128"

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
jobName = "job_0005_n128_abaqus"
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

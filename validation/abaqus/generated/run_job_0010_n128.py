# -*- coding: utf-8 -*-
# Auto-generated Abaqus CAE script for job: job_0010_n128
# Run with project Python (abqpy): python this_file; abqpy saveAs() launches Abaqus.
# Outputs CSV to: C:\Users\s1834431\Programs\fem_model\post_processing\validation_visualisers\abaqus_results\job_0010_n128

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
DISTRIBUTED_LOADS = [[0.0, 0.0, 0.0, 0.0, -0.0, 0.0, 0.0, 0.0, 0.0], [0.015625, 0.0, 0.0, 0.0, -3.90625, 0.0, 0.0, 0.0, 0.0], [0.03125, 0.0, 0.0, 0.0, -7.8125, 0.0, 0.0, 0.0, 0.0], [0.046875, 0.0, 0.0, 0.0, -11.71875, 0.0, 0.0, 0.0, 0.0], [0.0625, 0.0, 0.0, 0.0, -15.625, 0.0, 0.0, 0.0, 0.0], [0.078125, 0.0, 0.0, 0.0, -19.53125, 0.0, 0.0, 0.0, 0.0], [0.09375, 0.0, 0.0, 0.0, -23.4375, 0.0, 0.0, 0.0, 0.0], [0.109375, 0.0, 0.0, 0.0, -27.34375, 0.0, 0.0, 0.0, 0.0], [0.125, 0.0, 0.0, 0.0, -31.25, 0.0, 0.0, 0.0, 0.0], [0.140625, 0.0, 0.0, 0.0, -35.15625, 0.0, 0.0, 0.0, 0.0], [0.15625, 0.0, 0.0, 0.0, -39.0625, 0.0, 0.0, 0.0, 0.0], [0.171875, 0.0, 0.0, 0.0, -42.96875, 0.0, 0.0, 0.0, 0.0], [0.1875, 0.0, 0.0, 0.0, -46.875, 0.0, 0.0, 0.0, 0.0], [0.203125, 0.0, 0.0, 0.0, -50.78125, 0.0, 0.0, 0.0, 0.0], [0.21875, 0.0, 0.0, 0.0, -54.6875, 0.0, 0.0, 0.0, 0.0], [0.234375, 0.0, 0.0, 0.0, -58.59375, 0.0, 0.0, 0.0, 0.0], [0.25, 0.0, 0.0, 0.0, -62.5, 0.0, 0.0, 0.0, 0.0], [0.265625, 0.0, 0.0, 0.0, -66.40625, 0.0, 0.0, 0.0, 0.0], [0.28125, 0.0, 0.0, 0.0, -70.3125, 0.0, 0.0, 0.0, 0.0], [0.296875, 0.0, 0.0, 0.0, -74.21875, 0.0, 0.0, 0.0, 0.0], [0.3125, 0.0, 0.0, 0.0, -78.125, 0.0, 0.0, 0.0, 0.0], [0.328125, 0.0, 0.0, 0.0, -82.03125, 0.0, 0.0, 0.0, 0.0], [0.34375, 0.0, 0.0, 0.0, -85.9375, 0.0, 0.0, 0.0, 0.0], [0.359375, 0.0, 0.0, 0.0, -89.84375, 0.0, 0.0, 0.0, 0.0], [0.375, 0.0, 0.0, 0.0, -93.75, 0.0, 0.0, 0.0, 0.0], [0.390625, 0.0, 0.0, 0.0, -97.65625, 0.0, 0.0, 0.0, 0.0], [0.40625, 0.0, 0.0, 0.0, -101.5625, 0.0, 0.0, 0.0, 0.0], [0.421875, 0.0, 0.0, 0.0, -105.46875, 0.0, 0.0, 0.0, 0.0], [0.4375, 0.0, 0.0, 0.0, -109.375, 0.0, 0.0, 0.0, 0.0], [0.453125, 0.0, 0.0, 0.0, -113.28125, 0.0, 0.0, 0.0, 0.0], [0.46875, 0.0, 0.0, 0.0, -117.1875, 0.0, 0.0, 0.0, 0.0], [0.484375, 0.0, 0.0, 0.0, -121.09375, 0.0, 0.0, 0.0, 0.0], [0.5, 0.0, 0.0, 0.0, -125.0, 0.0, 0.0, 0.0, 0.0], [0.515625, 0.0, 0.0, 0.0, -128.90625, 0.0, 0.0, 0.0, 0.0], [0.53125, 0.0, 0.0, 0.0, -132.8125, 0.0, 0.0, 0.0, 0.0], [0.546875, 0.0, 0.0, 0.0, -136.71875, 0.0, 0.0, 0.0, 0.0], [0.5625, 0.0, 0.0, 0.0, -140.625, 0.0, 0.0, 0.0, 0.0], [0.578125, 0.0, 0.0, 0.0, -144.53125, 0.0, 0.0, 0.0, 0.0], [0.59375, 0.0, 0.0, 0.0, -148.4375, 0.0, 0.0, 0.0, 0.0], [0.609375, 0.0, 0.0, 0.0, -152.34375, 0.0, 0.0, 0.0, 0.0], [0.625, 0.0, 0.0, 0.0, -156.25, 0.0, 0.0, 0.0, 0.0], [0.640625, 0.0, 0.0, 0.0, -160.15625, 0.0, 0.0, 0.0, 0.0], [0.65625, 0.0, 0.0, 0.0, -164.0625, 0.0, 0.0, 0.0, 0.0], [0.671875, 0.0, 0.0, 0.0, -167.96875, 0.0, 0.0, 0.0, 0.0], [0.6875, 0.0, 0.0, 0.0, -171.875, 0.0, 0.0, 0.0, 0.0], [0.703125, 0.0, 0.0, 0.0, -175.78125, 0.0, 0.0, 0.0, 0.0], [0.71875, 0.0, 0.0, 0.0, -179.6875, 0.0, 0.0, 0.0, 0.0], [0.734375, 0.0, 0.0, 0.0, -183.59375, 0.0, 0.0, 0.0, 0.0], [0.75, 0.0, 0.0, 0.0, -187.5, 0.0, 0.0, 0.0, 0.0], [0.765625, 0.0, 0.0, 0.0, -191.40625, 0.0, 0.0, 0.0, 0.0], [0.78125, 0.0, 0.0, 0.0, -195.3125, 0.0, 0.0, 0.0, 0.0], [0.796875, 0.0, 0.0, 0.0, -199.21875, 0.0, 0.0, 0.0, 0.0], [0.8125, 0.0, 0.0, 0.0, -203.125, 0.0, 0.0, 0.0, 0.0], [0.828125, 0.0, 0.0, 0.0, -207.03125, 0.0, 0.0, 0.0, 0.0], [0.84375, 0.0, 0.0, 0.0, -210.9375, 0.0, 0.0, 0.0, 0.0], [0.859375, 0.0, 0.0, 0.0, -214.84375, 0.0, 0.0, 0.0, 0.0], [0.875, 0.0, 0.0, 0.0, -218.75, 0.0, 0.0, 0.0, 0.0], [0.890625, 0.0, 0.0, 0.0, -222.65625, 0.0, 0.0, 0.0, 0.0], [0.90625, 0.0, 0.0, 0.0, -226.5625, 0.0, 0.0, 0.0, 0.0], [0.921875, 0.0, 0.0, 0.0, -230.46875, 0.0, 0.0, 0.0, 0.0], [0.9375, 0.0, 0.0, 0.0, -234.375, 0.0, 0.0, 0.0, 0.0], [0.953125, 0.0, 0.0, 0.0, -238.28125, 0.0, 0.0, 0.0, 0.0], [0.96875, 0.0, 0.0, 0.0, -242.1875, 0.0, 0.0, 0.0, 0.0], [0.984375, 0.0, 0.0, 0.0, -246.09375, 0.0, 0.0, 0.0, 0.0], [1.0, 0.0, 0.0, 0.0, -250.0, 0.0, 0.0, 0.0, 0.0], [1.015625, 0.0, 0.0, 0.0, -253.90625, 0.0, 0.0, 0.0, 0.0], [1.03125, 0.0, 0.0, 0.0, -257.8125, 0.0, 0.0, 0.0, 0.0], [1.046875, 0.0, 0.0, 0.0, -261.71875, 0.0, 0.0, 0.0, 0.0], [1.0625, 0.0, 0.0, 0.0, -265.625, 0.0, 0.0, 0.0, 0.0], [1.078125, 0.0, 0.0, 0.0, -269.53125, 0.0, 0.0, 0.0, 0.0], [1.09375, 0.0, 0.0, 0.0, -273.4375, 0.0, 0.0, 0.0, 0.0], [1.109375, 0.0, 0.0, 0.0, -277.34375, 0.0, 0.0, 0.0, 0.0], [1.125, 0.0, 0.0, 0.0, -281.25, 0.0, 0.0, 0.0, 0.0], [1.140625, 0.0, 0.0, 0.0, -285.15625, 0.0, 0.0, 0.0, 0.0], [1.15625, 0.0, 0.0, 0.0, -289.0625, 0.0, 0.0, 0.0, 0.0], [1.171875, 0.0, 0.0, 0.0, -292.96875, 0.0, 0.0, 0.0, 0.0], [1.1875, 0.0, 0.0, 0.0, -296.875, 0.0, 0.0, 0.0, 0.0], [1.203125, 0.0, 0.0, 0.0, -300.78125, 0.0, 0.0, 0.0, 0.0], [1.21875, 0.0, 0.0, 0.0, -304.6875, 0.0, 0.0, 0.0, 0.0], [1.234375, 0.0, 0.0, 0.0, -308.59375, 0.0, 0.0, 0.0, 0.0], [1.25, 0.0, 0.0, 0.0, -312.5, 0.0, 0.0, 0.0, 0.0], [1.265625, 0.0, 0.0, 0.0, -316.40625, 0.0, 0.0, 0.0, 0.0], [1.28125, 0.0, 0.0, 0.0, -320.3125, 0.0, 0.0, 0.0, 0.0], [1.296875, 0.0, 0.0, 0.0, -324.21875, 0.0, 0.0, 0.0, 0.0], [1.3125, 0.0, 0.0, 0.0, -328.125, 0.0, 0.0, 0.0, 0.0], [1.328125, 0.0, 0.0, 0.0, -332.03125, 0.0, 0.0, 0.0, 0.0], [1.34375, 0.0, 0.0, 0.0, -335.9375, 0.0, 0.0, 0.0, 0.0], [1.359375, 0.0, 0.0, 0.0, -339.84375, 0.0, 0.0, 0.0, 0.0], [1.375, 0.0, 0.0, 0.0, -343.75, 0.0, 0.0, 0.0, 0.0], [1.390625, 0.0, 0.0, 0.0, -347.65625, 0.0, 0.0, 0.0, 0.0], [1.40625, 0.0, 0.0, 0.0, -351.5625, 0.0, 0.0, 0.0, 0.0], [1.421875, 0.0, 0.0, 0.0, -355.46875, 0.0, 0.0, 0.0, 0.0], [1.4375, 0.0, 0.0, 0.0, -359.375, 0.0, 0.0, 0.0, 0.0], [1.453125, 0.0, 0.0, 0.0, -363.28125, 0.0, 0.0, 0.0, 0.0], [1.46875, 0.0, 0.0, 0.0, -367.1875, 0.0, 0.0, 0.0, 0.0], [1.484375, 0.0, 0.0, 0.0, -371.09375, 0.0, 0.0, 0.0, 0.0], [1.5, 0.0, 0.0, 0.0, -375.0, 0.0, 0.0, 0.0, 0.0], [1.515625, 0.0, 0.0, 0.0, -378.90625, 0.0, 0.0, 0.0, 0.0], [1.53125, 0.0, 0.0, 0.0, -382.8125, 0.0, 0.0, 0.0, 0.0], [1.546875, 0.0, 0.0, 0.0, -386.71875, 0.0, 0.0, 0.0, 0.0], [1.5625, 0.0, 0.0, 0.0, -390.625, 0.0, 0.0, 0.0, 0.0], [1.578125, 0.0, 0.0, 0.0, -394.53125, 0.0, 0.0, 0.0, 0.0], [1.59375, 0.0, 0.0, 0.0, -398.4375, 0.0, 0.0, 0.0, 0.0], [1.609375, 0.0, 0.0, 0.0, -402.34375, 0.0, 0.0, 0.0, 0.0], [1.625, 0.0, 0.0, 0.0, -406.25, 0.0, 0.0, 0.0, 0.0], [1.640625, 0.0, 0.0, 0.0, -410.15625, 0.0, 0.0, 0.0, 0.0], [1.65625, 0.0, 0.0, 0.0, -414.0625, 0.0, 0.0, 0.0, 0.0], [1.671875, 0.0, 0.0, 0.0, -417.96875, 0.0, 0.0, 0.0, 0.0], [1.6875, 0.0, 0.0, 0.0, -421.875, 0.0, 0.0, 0.0, 0.0], [1.703125, 0.0, 0.0, 0.0, -425.78125, 0.0, 0.0, 0.0, 0.0], [1.71875, 0.0, 0.0, 0.0, -429.6875, 0.0, 0.0, 0.0, 0.0], [1.734375, 0.0, 0.0, 0.0, -433.59375, 0.0, 0.0, 0.0, 0.0], [1.75, 0.0, 0.0, 0.0, -437.5, 0.0, 0.0, 0.0, 0.0], [1.765625, 0.0, 0.0, 0.0, -441.40625, 0.0, 0.0, 0.0, 0.0], [1.78125, 0.0, 0.0, 0.0, -445.3125, 0.0, 0.0, 0.0, 0.0], [1.796875, 0.0, 0.0, 0.0, -449.21875, 0.0, 0.0, 0.0, 0.0], [1.8125, 0.0, 0.0, 0.0, -453.125, 0.0, 0.0, 0.0, 0.0], [1.828125, 0.0, 0.0, 0.0, -457.03125, 0.0, 0.0, 0.0, 0.0], [1.84375, 0.0, 0.0, 0.0, -460.9375, 0.0, 0.0, 0.0, 0.0], [1.859375, 0.0, 0.0, 0.0, -464.84375, 0.0, 0.0, 0.0, 0.0], [1.875, 0.0, 0.0, 0.0, -468.75, 0.0, 0.0, 0.0, 0.0], [1.890625, 0.0, 0.0, 0.0, -472.65625, 0.0, 0.0, 0.0, 0.0], [1.90625, 0.0, 0.0, 0.0, -476.5625, 0.0, 0.0, 0.0, 0.0], [1.921875, 0.0, 0.0, 0.0, -480.46875, 0.0, 0.0, 0.0, 0.0], [1.9375, 0.0, 0.0, 0.0, -484.375, 0.0, 0.0, 0.0, 0.0], [1.953125, 0.0, 0.0, 0.0, -488.28125, 0.0, 0.0, 0.0, 0.0], [1.96875, 0.0, 0.0, 0.0, -492.1875, 0.0, 0.0, 0.0, 0.0], [1.984375, 0.0, 0.0, 0.0, -496.09375, 0.0, 0.0, 0.0, 0.0], [2.0, 0.0, 0.0, 0.0, -500.0, 0.0, 0.0, 0.0, 0.0]]
DISTRIBUTED_EQUIVALENT_NODAL = [(1, 0.0, -0.009155273437500003, 0.0, 0.0, 0.0, 0.0), (2, 0.0, -0.061035156250000014, 0.0, 0.0, 0.0, 0.0), (3, 0.0, -0.12207031250000003, 0.0, 0.0, 0.0, 0.0), (4, 0.0, -0.18310546875000003, 0.0, 0.0, 0.0, 0.0), (5, 0.0, -0.24414062500000003, 0.0, 0.0, 0.0, 0.0), (6, 0.0, -0.30517578125, 0.0, 0.0, 0.0, 0.0), (7, 0.0, -0.36621093750000006, 0.0, 0.0, 0.0, 0.0), (8, 0.0, -0.42724609375, 0.0, 0.0, 0.0, 0.0), (9, 0.0, -0.48828125000000006, 0.0, 0.0, 0.0, 0.0), (10, 0.0, -0.54931640625, 0.0, 0.0, 0.0, 0.0), (11, 0.0, -0.6103515625, 0.0, 0.0, 0.0, 0.0), (12, 0.0, -0.6713867187500002, 0.0, 0.0, 0.0, 0.0), (13, 0.0, -0.7324218750000002, 0.0, 0.0, 0.0, 0.0), (14, 0.0, -0.7934570312500001, 0.0, 0.0, 0.0, 0.0), (15, 0.0, -0.8544921875000002, 0.0, 0.0, 0.0, 0.0), (16, 0.0, -0.9155273437500002, 0.0, 0.0, 0.0, 0.0), (17, 0.0, -0.9765625000000002, 0.0, 0.0, 0.0, 0.0), (18, 0.0, -1.0375976562500004, 0.0, 0.0, 0.0, 0.0), (19, 0.0, -1.0986328125000002, 0.0, 0.0, 0.0, 0.0), (20, 0.0, -1.15966796875, 0.0, 0.0, 0.0, 0.0), (21, 0.0, -1.2207031250000004, 0.0, 0.0, 0.0, 0.0), (22, 0.0, -1.2817382812500002, 0.0, 0.0, 0.0, 0.0), (23, 0.0, -1.3427734375000004, 0.0, 0.0, 0.0, 0.0), (24, 0.0, -1.4038085937500004, 0.0, 0.0, 0.0, 0.0), (25, 0.0, -1.46484375, 0.0, 0.0, 0.0, 0.0), (26, 0.0, -1.5258789062500004, 0.0, 0.0, 0.0, 0.0), (27, 0.0, -1.5869140625000004, 0.0, 0.0, 0.0, 0.0), (28, 0.0, -1.64794921875, 0.0, 0.0, 0.0, 0.0), (29, 0.0, -1.7089843750000004, 0.0, 0.0, 0.0, 0.0), (30, 0.0, -1.7700195312500004, 0.0, 0.0, 0.0, 0.0), (31, 0.0, -1.8310546875000004, 0.0, 0.0, 0.0, 0.0), (32, 0.0, -1.8920898437500004, 0.0, 0.0, 0.0, 0.0), (33, 0.0, -1.9531250000000004, 0.0, 0.0, 0.0, 0.0), (34, 0.0, -2.01416015625, 0.0, 0.0, 0.0, 0.0), (35, 0.0, -2.0751953125000004, 0.0, 0.0, 0.0, 0.0), (36, 0.0, -2.1362304687500004, 0.0, 0.0, 0.0, 0.0), (37, 0.0, -2.1972656250000004, 0.0, 0.0, 0.0, 0.0), (38, 0.0, -2.25830078125, 0.0, 0.0, 0.0, 0.0), (39, 0.0, -2.3193359375, 0.0, 0.0, 0.0, 0.0), (40, 0.0, -2.3803710937500004, 0.0, 0.0, 0.0, 0.0), (41, 0.0, -2.4414062500000004, 0.0, 0.0, 0.0, 0.0), (42, 0.0, -2.5024414062500004, 0.0, 0.0, 0.0, 0.0), (43, 0.0, -2.563476562500001, 0.0, 0.0, 0.0, 0.0), (44, 0.0, -2.624511718750001, 0.0, 0.0, 0.0, 0.0), (45, 0.0, -2.685546875, 0.0, 0.0, 0.0, 0.0), (46, 0.0, -2.7465820312500004, 0.0, 0.0, 0.0, 0.0), (47, 0.0, -2.8076171875000004, 0.0, 0.0, 0.0, 0.0), (48, 0.0, -2.8686523437500004, 0.0, 0.0, 0.0, 0.0), (49, 0.0, -2.929687500000001, 0.0, 0.0, 0.0, 0.0), (50, 0.0, -2.990722656250001, 0.0, 0.0, 0.0, 0.0), (51, 0.0, -3.0517578125, 0.0, 0.0, 0.0, 0.0), (52, 0.0, -3.1127929687500004, 0.0, 0.0, 0.0, 0.0), (53, 0.0, -3.173828125, 0.0, 0.0, 0.0, 0.0), (54, 0.0, -3.2348632812500004, 0.0, 0.0, 0.0, 0.0), (55, 0.0, -3.2958984375000004, 0.0, 0.0, 0.0, 0.0), (56, 0.0, -3.35693359375, 0.0, 0.0, 0.0, 0.0), (57, 0.0, -3.417968750000001, 0.0, 0.0, 0.0, 0.0), (58, 0.0, -3.47900390625, 0.0, 0.0, 0.0, 0.0), (59, 0.0, -3.5400390625, 0.0, 0.0, 0.0, 0.0), (60, 0.0, -3.6010742187500004, 0.0, 0.0, 0.0, 0.0), (61, 0.0, -3.662109375000001, 0.0, 0.0, 0.0, 0.0), (62, 0.0, -3.723144531250001, 0.0, 0.0, 0.0, 0.0), (63, 0.0, -3.7841796875, 0.0, 0.0, 0.0, 0.0), (64, 0.0, -3.84521484375, 0.0, 0.0, 0.0, 0.0), (65, 0.0, -3.9062500000000004, 0.0, 0.0, 0.0, 0.0), (66, 0.0, -3.967285156250001, 0.0, 0.0, 0.0, 0.0), (67, 0.0, -4.028320312500001, 0.0, 0.0, 0.0, 0.0), (68, 0.0, -4.089355468750002, 0.0, 0.0, 0.0, 0.0), (69, 0.0, -4.150390625000002, 0.0, 0.0, 0.0, 0.0), (70, 0.0, -4.211425781250002, 0.0, 0.0, 0.0, 0.0), (71, 0.0, -4.272460937500001, 0.0, 0.0, 0.0, 0.0), (72, 0.0, -4.333496093750002, 0.0, 0.0, 0.0, 0.0), (73, 0.0, -4.394531250000002, 0.0, 0.0, 0.0, 0.0), (74, 0.0, -4.455566406250001, 0.0, 0.0, 0.0, 0.0), (75, 0.0, -4.516601562500002, 0.0, 0.0, 0.0, 0.0), (76, 0.0, -4.577636718750002, 0.0, 0.0, 0.0, 0.0), (77, 0.0, -4.638671875000002, 0.0, 0.0, 0.0, 0.0), (78, 0.0, -4.69970703125, 0.0, 0.0, 0.0, 0.0), (79, 0.0, -4.760742187500002, 0.0, 0.0, 0.0, 0.0), (80, 0.0, -4.821777343750002, 0.0, 0.0, 0.0, 0.0), (81, 0.0, -4.882812500000001, 0.0, 0.0, 0.0, 0.0), (82, 0.0, -4.943847656250002, 0.0, 0.0, 0.0, 0.0), (83, 0.0, -5.004882812500002, 0.0, 0.0, 0.0, 0.0), (84, 0.0, -5.065917968750002, 0.0, 0.0, 0.0, 0.0), (85, 0.0, -5.126953125000002, 0.0, 0.0, 0.0, 0.0), (86, 0.0, -5.187988281250002, 0.0, 0.0, 0.0, 0.0), (87, 0.0, -5.249023437500002, 0.0, 0.0, 0.0, 0.0), (88, 0.0, -5.310058593750001, 0.0, 0.0, 0.0, 0.0), (89, 0.0, -5.371093750000002, 0.0, 0.0, 0.0, 0.0), (90, 0.0, -5.432128906250002, 0.0, 0.0, 0.0, 0.0), (91, 0.0, -5.493164062500002, 0.0, 0.0, 0.0, 0.0), (92, 0.0, -5.554199218750002, 0.0, 0.0, 0.0, 0.0), (93, 0.0, -5.615234375, 0.0, 0.0, 0.0, 0.0), (94, 0.0, -5.676269531250002, 0.0, 0.0, 0.0, 0.0), (95, 0.0, -5.737304687500002, 0.0, 0.0, 0.0, 0.0), (96, 0.0, -5.798339843750002, 0.0, 0.0, 0.0, 0.0), (97, 0.0, -5.859375000000002, 0.0, 0.0, 0.0, 0.0), (98, 0.0, -5.92041015625, 0.0, 0.0, 0.0, 0.0), (99, 0.0, -5.9814453125, 0.0, 0.0, 0.0, 0.0), (100, 0.0, -6.04248046875, 0.0, 0.0, 0.0, 0.0), (101, 0.0, -6.103515625000002, 0.0, 0.0, 0.0, 0.0), (102, 0.0, -6.164550781250002, 0.0, 0.0, 0.0, 0.0), (103, 0.0, -6.225585937500002, 0.0, 0.0, 0.0, 0.0), (104, 0.0, -6.286621093750002, 0.0, 0.0, 0.0, 0.0), (105, 0.0, -6.347656250000002, 0.0, 0.0, 0.0, 0.0), (106, 0.0, -6.408691406250002, 0.0, 0.0, 0.0, 0.0), (107, 0.0, -6.469726562500001, 0.0, 0.0, 0.0, 0.0), (108, 0.0, -6.530761718750002, 0.0, 0.0, 0.0, 0.0), (109, 0.0, -6.591796875000002, 0.0, 0.0, 0.0, 0.0), (110, 0.0, -6.652832031250002, 0.0, 0.0, 0.0, 0.0), (111, 0.0, -6.713867187500002, 0.0, 0.0, 0.0, 0.0), (112, 0.0, -6.774902343750002, 0.0, 0.0, 0.0, 0.0), (113, 0.0, -6.835937500000003, 0.0, 0.0, 0.0, 0.0), (114, 0.0, -6.896972656250002, 0.0, 0.0, 0.0, 0.0), (115, 0.0, -6.958007812500002, 0.0, 0.0, 0.0, 0.0), (116, 0.0, -7.019042968750001, 0.0, 0.0, 0.0, 0.0), (117, 0.0, -7.080078125000002, 0.0, 0.0, 0.0, 0.0), (118, 0.0, -7.141113281250002, 0.0, 0.0, 0.0, 0.0), (119, 0.0, -7.202148437500002, 0.0, 0.0, 0.0, 0.0), (120, 0.0, -7.263183593750002, 0.0, 0.0, 0.0, 0.0), (121, 0.0, -7.324218750000003, 0.0, 0.0, 0.0, 0.0), (122, 0.0, -7.385253906250002, 0.0, 0.0, 0.0, 0.0), (123, 0.0, -7.446289062500002, 0.0, 0.0, 0.0, 0.0), (124, 0.0, -7.507324218750002, 0.0, 0.0, 0.0, 0.0), (125, 0.0, -7.568359375000002, 0.0, 0.0, 0.0, 0.0), (126, 0.0, -7.629394531250001, 0.0, 0.0, 0.0, 0.0), (127, 0.0, -7.690429687500001, 0.0, 0.0, 0.0, 0.0), (128, 0.0, -7.751464843750002, 0.0, 0.0, 0.0, 0.0), (129, 0.0, -3.897094726562501, 0.0, 0.0, 0.0, 0.0)]
OUT_CSV_DIR = r"C:\\Users\\s1834431\\Programs\\fem_model\\post_processing\\validation_visualisers\\abaqus_results\\job_0010_n128"

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
jobName = "job_0010_n128_abaqus"
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

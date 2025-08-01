# blade_structural_model\internal_force\blade_load_profile_reader.py

from pathlib import Path
import pandas as pd
import numpy as np

class LoadProfileReader:
    def __init__(self, csv_path: Path, R: float, L: float):
        self.csv_path = csv_path
        self.R = R
        self.L = L
        self.offset = R - L
        self.data = {}
        self.rR = None
        self.rR_locii = None
        self.xL = None
        self.xL_locii = None
        self.dx = None

    def read(self):
        df = pd.read_csv(self.csv_path)
        df.columns = df.columns.str.strip("[]")

        x = df["x"].to_numpy()
        self.dx = float(np.mean(np.diff(x)))
        self.rR = (x + self.offset) / self.R
        self.xL = x / self.L

        self.data = {
            "F_y": df["F_y"].to_numpy(),
            "F_z": -df["F_z"].to_numpy(),
            "M_x": df["M_x"].to_numpy(),
        }

        self.rR_locii = {
            "rotor_origin": 0.0,
            "blade_root": self.offset / self.R,
            "blade_tip": 1.0,
        }

        self.xL_locii = {
            "blade_root": 0.0,
            "blade_tip": 1.0,
        }

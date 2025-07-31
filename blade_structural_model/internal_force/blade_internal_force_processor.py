import numpy as np

class InternalForceProcessor:
    def __init__(self, loads: dict, dx: float):
        self.loads = loads
        self.dx = dx
        self.results = {}

    def cumtrapz_uniform(self, y: np.ndarray) -> np.ndarray:
        return np.concatenate(([0.0], np.cumsum((y[:-1] + y[1:]) * 0.5 * self.dx)))

    def compute(self):
        qy = self.loads["F_y"][::-1]
        qz = self.loads["F_z"][::-1]
        mx = self.loads["M_x"][::-1]

        Vy = self.cumtrapz_uniform(qy)[::-1]
        Vz = self.cumtrapz_uniform(qz)[::-1]
        Mz = -self.cumtrapz_uniform(Vy[::-1])[::-1]
        My = -self.cumtrapz_uniform(Vz[::-1])[::-1]
        T  = self.cumtrapz_uniform(mx)[::-1]

        self.results = {
            "V_y": Vy,
            "M_z": Mz,
            "V_z": Vz,
            "M_y": My,
            "T": T,
        }

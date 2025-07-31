from dataclasses import dataclass
import numpy as np
import pandas as pd

@dataclass
class BladeInternalForceRow:
    tsr: str
    r_over_R: np.ndarray

    # Distributed Loads (inputs)
    f_y: np.ndarray      # Flapwise force [N/m]
    f_z: np.ndarray      # Edgewise force [N/m]
    m_x: np.ndarray      # Torsional moment [Nm/m]

    # Internal Resultants (outputs)
    V_y: np.ndarray      # Flapwise shear [N]
    M_z: np.ndarray      # Flapwise bending moment [Nm]
    V_z: np.ndarray      # Edgewise shear [N]
    M_y: np.ndarray      # Edgewise bending moment [Nm]
    T: np.ndarray        # Torsional resultant [Nm]

    def to_dataframe(self) -> pd.DataFrame:
        return pd.DataFrame({
            "r/R": self.r_over_R,
            "f_y [N/m]": self.f_y,
            "V_y [N]": self.V_y,
            "M_z [Nm]": self.M_z,
            "f_z [N/m]": self.f_z,
            "V_z [N]": self.V_z,
            "M_y [Nm]": self.M_y,
            "m_x [Nm/m]": self.m_x,
            "T [Nm]": self.T,
        })
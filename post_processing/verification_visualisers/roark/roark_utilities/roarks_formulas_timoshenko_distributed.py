# post_processing/verification_visualisers/roark_utilities/roarks_formulas_timoshenko_distributed.py
"""
Timoshenko beam formulas for cantilever with distributed loads.

Same load types as Euler–Bernoulli (Roark): UDL, triangular, parabolic.
V and M unchanged; θ = ∫ M/(EI) dx unchanged. Deflection adds shear:
  du/dx = θ + V/(k_s A G)  =>  u = u_EB + (1/(k_s A G)) ∫_0^x V(ξ) dξ.

Closed-form ∫_0^x V(ξ) dξ:
  UDL:       -w(Lx - x²/2)
  Triangular: -w/(2L)(L²x - x³/3)
  Parabolic:  -w/(3L²)(L³x - x⁴/4)
"""

import numpy as np

from roarks_formulas_euler_bernoulli_distributed import RoarksFormulaeDistributedLoad  # type: ignore


DEFAULT_K_S: float = 5.0 / 6.0


class TimoshenkoFormulaeDistributedLoad(RoarksFormulaeDistributedLoad):
    """
    Timoshenko beam response for distributed loads.
    Extends Roark (E–B) with shear deflection; requires A, G, and optional k_s.
    """

    def __init__(
        self,
        L: float,
        E: float,
        I: float,
        w: float,
        A: float,
        G: float,
        load_type: str,
        k_s: float = DEFAULT_K_S,
    ):
        super().__init__(L, E, I, w, load_type)
        self.A = A
        self.G = G
        self.k_s = k_s
        if self.A <= 0 or self.G <= 0 or not (0 < self.k_s <= 1):
            raise ValueError("A, G must be positive; k_s in (0, 1]")

    def _integral_V(self, x: np.ndarray) -> np.ndarray:
        """∫_0^x V(ξ) dξ (closed form)."""
        if self.load_type == "udl":
            return -self.w * (self.L * x - x**2 / 2)
        elif self.load_type == "triangular":
            return -self.w / (2 * self.L) * (self.L**2 * x - x**3 / 3)
        else:  # parabolic
            return -self.w / (3 * self.L**2) * (self.L**3 * x - x**4 / 4)

    def deflection(self, x: np.ndarray) -> np.ndarray:
        """Deflection u_y(x): E–B deflection plus shear contribution."""
        if not np.all(np.diff(x) >= 0):
            raise ValueError("Input array x must be sorted in ascending order")
        u_eb = super().deflection(x)
        integral_V = self._integral_V(x)
        u_shear = integral_V / (self.k_s * self.A * self.G)
        return u_eb + u_shear

    def response(self, x: np.ndarray) -> dict:
        return {
            "intensity": self.intensity(x),
            "shear": self.shear(x),
            "moment": self.moment(x),
            "rotation": self.rotation(x),
            "deflection": self.deflection(x),
        }


def timoshenko_distributed_load_response(
    x: np.ndarray,
    L: float,
    E: float,
    I: float,
    w: float,
    A: float,
    G: float,
    load_type: str,
    k_s: float = DEFAULT_K_S,
) -> dict:
    """
    Convenience wrapper: Timoshenko distributed-load response.
    load_type: 'udl', 'triangular', or 'parabolic'.
    Returns dict: {intensity, shear, moment, rotation, deflection}.
    """
    solver = TimoshenkoFormulaeDistributedLoad(L, E, I, w, A, G, load_type, k_s=k_s)
    return solver.response(x)

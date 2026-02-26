# post_processing/verification_visualisers/roark_utilities/roarks_formulas_timoshenko_point.py
"""
Timoshenko beam formulas for cantilever with point load.

Same load cases as Euler–Bernoulli (Roark): end, mid-span, quarter-span.
V and M are unchanged; θ = ∫ M/(EI) dx is unchanged. Deflection adds shear:
  du/dx = θ + V/(k_s A G)  =>  u = u_EB + ∫_0^x V(ξ)/(k_s A G) dξ.

Point load P at x = a: V(x) = -P for x < a, 0 for x ≥ a.
  => u_shear(x) = -P·min(x,a)/(k_s A G)  [in physical sign: adds downward deflection].
Convention: returned deflection/rotation match roarks_formulas_euler_bernoulli_point (same sign as FEM).
"""

import numpy as np

from roarks_formulas_euler_bernoulli_point import RoarksFormulaePointLoad  # type: ignore


# Default shear correction factor (e.g. 5/6 for rectangular cross-section)
DEFAULT_K_S: float = 5.0 / 6.0


class TimoshenkoFormulaePointLoad(RoarksFormulaePointLoad):
    """
    Timoshenko beam response for cantilever with point load.
    Extends Roark (E–B) with shear deflection; requires A, G, and optional k_s.
    """

    def __init__(self, L, E, I, P, A, G, k_s=DEFAULT_K_S, a=None, load_type=None):
        super().__init__(L, E, I, P, a=a, load_type=load_type)
        self.A = A
        self.G = G
        self.k_s = k_s
        if self.A <= 0 or self.G <= 0 or not (0 < self.k_s <= 1):
            raise ValueError("A, G must be positive; k_s in (0, 1]")

    def deflection(self, x):
        """Deflection u_y(x): E–B deflection plus shear contribution (same sign convention as Roark)."""
        u_eb = super().deflection(x)  # already sign-flipped to match FEM
        # Shear: u_shear_physical = P·min(x,a)/(k_s A G); we return -u_physical => subtract it from E–B return
        x_safe = np.where(x < self.a, x, self.a)
        u_shear_term = self.P * x_safe / (self.k_s * self.A * self.G)
        return u_eb - u_shear_term

    # rotation, shear, moment, intensity unchanged from E–B
    def response(self, x):
        return {
            "intensity": self.intensity(x),
            "shear": self.shear(x),
            "moment": self.moment(x),
            "rotation": self.rotation(x),
            "deflection": self.deflection(x),
        }


def timoshenko_point_load_response(x, L, E, I, P, A, G, k_s=DEFAULT_K_S, load_type="end"):
    """
    Convenience wrapper: Timoshenko point-load response.
    load_type: 'end', 'mid', or 'quarter'.
    Returns dict: {intensity, shear, moment, rotation, deflection}.
    """
    solver = TimoshenkoFormulaePointLoad(L, E, I, P, A, G, k_s=k_s, load_type=load_type)
    return solver.response(x)

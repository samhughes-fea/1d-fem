# post_processing/verification_visualisers/roark_utilities/roarks_formulas_euler_bernoulli_point.py

import numpy as np

class RoarksFormulaePointLoad:
    """
    Calculates beam responses for a cantilever with point load using Roark's formulas.
    
    Parameters:
        L (float): Beam length
        E (float): Young's modulus
        I (float): Area moment of inertia
        P (float): Point load magnitude (positive downward)
        a (float, optional): Load position from fixed end (0 ≤ a ≤ L)
        load_type (str, optional): Predefined load position ('end', 'mid', 'quarter')
    """
    def __init__(self, L, E, I, P, a=None, load_type=None):
        self.L = L
        self.E = E
        self.I = I
        self.P = P
        
        if a is not None:
            self.a = a
        elif load_type is not None:
            self.set_load_type(load_type)
        else:
            raise ValueError("Must specify either 'a' or 'load_type'")
            
        self.validate_parameters()
        
    def validate_parameters(self):
        """Check for valid physical parameters and load position"""
        if self.L <= 0:
            raise ValueError("Beam length L must be positive")
        if self.E <= 0:
            raise ValueError("Young's modulus E must be positive")
        if self.I <= 0:
            raise ValueError("Moment of inertia I must be positive")
        if not (0 <= self.a <= self.L):
            raise ValueError(f"Load position a={self.a} must be in [0, L={self.L}]")

    def set_load_type(self, load_type):
        """Set load position using predefined types"""
        if load_type == "end":
            self.a = self.L
        elif load_type == "mid":
            self.a = self.L / 2
        elif load_type == "quarter":
            self.a = self.L / 4
        else:
            raise ValueError("load_type must be 'end', 'mid', or 'quarter'")
    
    def intensity(self, x):
        """Point load intensity (Dirac delta approximation)"""
        q = np.zeros_like(x)
        idx = np.argmin(np.abs(x - self.a))
        q[idx] = -self.P  # Negative for downward load
        return q
    
    def shear(self, x):
        """Shear force distribution V(x)"""
        return np.where(x < self.a, -self.P, 0.0)
    
    def moment(self, x):
        """Bending moment distribution M(x)"""
        return np.where(x < self.a, -self.P * (self.a - x), 0.0)
    
    def rotation(self, x):
        """Slope/rotation θ_z(x). Sign flipped so positive load → same sign as FEM (θ_z)."""
        theta_a = -(self.P * self.a**2) / (2 * self.E * self.I)
        region1 = (x < self.a)
        theta = np.empty_like(x)
        theta[region1] = -(self.P * x[region1]) * (2*self.a - x[region1]) / (2 * self.E * self.I)
        theta[~region1] = theta_a
        return -theta

    def deflection(self, x):
        """Deflection u_y(x). Sign flipped so positive load → same sign as FEM (u_y)."""
        u_a = -(self.P * self.a**3) / (3 * self.E * self.I)
        theta_a = -(self.P * self.a**2) / (2 * self.E * self.I)
        region1 = (x < self.a)
        u = np.empty_like(x)
        u[region1] = -(self.P * x[region1]**2) * (3*self.a - x[region1]) / (6 * self.E * self.I)
        u[~region1] = u_a + theta_a * (x[~region1] - self.a)
        return -u
    
    def response(self, x):
        """
        Compute complete beam response as a dictionary of:
        {intensity, shear, moment, rotation, deflection}
        """
        return {
            "intensity": self.intensity(x),
            "shear": self.shear(x),
            "moment": self.moment(x),
            "rotation": self.rotation(x),
            "deflection": self.deflection(x)
        }


def roark_point_load_response(x, L, E, I, P, load_type):
    """
    Convenience wrapper for RoarksFormulaePointLoad.response().
    load_type: 'end', 'mid', or 'quarter'.
    Returns dict: {intensity, shear, moment, rotation, deflection}.
    """
    solver = RoarksFormulaePointLoad(L, E, I, P, load_type=load_type)
    return solver.response(x)

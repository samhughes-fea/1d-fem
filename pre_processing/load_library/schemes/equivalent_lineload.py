import numpy as np

# note this could be implemented in element formulation logic

class RigidLoadTransformer:
    def __init__(self, cop_path: np.ndarray, cop_line_load: np.ndarray, centroid_path: np.ndarray):
        """
        Parameters:
        - cop_path:       (N, 3) array of [x, y, z] for CoP line
        - cop_line_load:  (N, 9) array with columns [x, y, z, Fx, Fy, Fz, Mx, My, Mz]
        - centroid_path:  (N, 3) array of [x, y, z] for centroid line
        """
        self._validate_inputs(cop_path, cop_line_load, centroid_path)
        self.cop_path = cop_path
        self.cop_line_load = cop_line_load
        self.centroid_path = centroid_path

        # Perform transformation immediately
        self.transformed_loads = self._transform()

    def _validate_inputs(self, cop_path, cop_line_load, centroid_path):
        if cop_path.shape != centroid_path.shape or cop_path.shape[1] != 3:
            raise ValueError("cop_path and centroid_path must be (N, 3) arrays of the same shape.")
        if cop_line_load.shape[0] != cop_path.shape[0] or cop_line_load.shape[1] != 9:
            raise ValueError("cop_line_load must be of shape (N, 9) matching the path length.")

    def _transform(self) -> np.ndarray:
        """
        Internal method to transform loads from CoP path to centroid path.

        Returns:
        - transformed_loads: (N, 9) array at centroid path: [x, y, z, Fx, Fy, Fz, Mx, My, Mz]
        """
        # Extract forces and moments
        F = self.cop_line_load[:, 3:6]  # Fx, Fy, Fz
        M = self.cop_line_load[:, 6:9]  # Mx, My, Mz

        # Compute offset vector from CoP to centroid
        r = self.centroid_path - self.cop_path  # shape (N, 3)

        # Compute r × F
        r_cross_F = np.cross(r, F)

        # Add moment correction
        M_transformed = M + r_cross_F

        # Combine: new position + same forces + updated moments
        return np.hstack((self.centroid_path, F, M_transformed))
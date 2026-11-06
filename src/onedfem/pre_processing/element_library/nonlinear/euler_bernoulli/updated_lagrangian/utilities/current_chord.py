# Current-chord frame for Updated-Lagrangian beams (local x along the live chord).

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from onedfem.pre_processing.element_library.element_frame import node_dof_rotation, straight_frame


@dataclass
class CurrentChordFrameOperator:
    """Map increment DOFs between global axes and the current chord."""

    def rotation(self, node_coords: np.ndarray) -> np.ndarray:
        c = np.asarray(node_coords, dtype=np.float64)
        return straight_frame(c[0], c[1])

    def dof_map(self, node_coords: np.ndarray, n_dof: int) -> np.ndarray:
        R = self.rotation(node_coords)
        return node_dof_rotation(R, R, int(n_dof))

    def pad_sigma(self, K_mat: np.ndarray, K_sigma: np.ndarray) -> np.ndarray:
        n = K_mat.shape[0]
        if K_sigma.shape[0] == n:
            return K_sigma
        Ks = np.zeros_like(K_mat)
        m = min(n, K_sigma.shape[0])
        Ks[:m, :m] = K_sigma[:m, :m]
        return Ks


class UpdatedLagrangianChordMixin:
    """Local ``K_mat``/``K_σ`` on the current chord, rotated to global increment DOFs."""

    def _rebuild_ul_T(self) -> None:
        n = int(getattr(self, "_U_ref", np.zeros(12)).size)
        self._ul_T = self.current_chord_frame_operator.dof_map(self.node_coords, n)

    def _ul_local_increment(self, U_e: np.ndarray):
        u = self._increment(U_e)
        T = getattr(self, "_ul_T", None)
        if T is None or T.shape[0] != u.size:
            self._rebuild_ul_T()
            T = self._ul_T
        return T @ u, T

    def tangent_stiffness_matrix(self, U_e: np.ndarray) -> np.ndarray:
        u_loc, T = self._ul_local_increment(U_e)
        K_mat = self._get_K_0_current()
        K_sigma = self.current_chord_frame_operator.pad_sigma(
            K_mat, np.asarray(self.linear_geometric_stiffness_matrix(u_loc), dtype=np.float64)
        )
        K_loc = K_mat + K_sigma
        return T.T @ K_loc @ T

    def internal_force_vector(self, U_e: np.ndarray) -> np.ndarray:
        u_loc, T = self._ul_local_increment(U_e)
        K_mat = self._get_K_0_current()
        K_sigma = self.current_chord_frame_operator.pad_sigma(
            K_mat, np.asarray(self.linear_geometric_stiffness_matrix(u_loc), dtype=np.float64)
        )
        return T.T @ (K_mat @ u_loc + K_sigma @ u_loc)

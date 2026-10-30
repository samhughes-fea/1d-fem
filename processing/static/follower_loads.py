# processing/static/follower_loads.py
"""Nodal follower loads: forces and moments that turn with the node they act on.

A follower load is given in the node's undeformed orientation (``follower_load.txt``,
global components) and turns with the node's rotation ``R = exp(theta)``, ``theta`` the
node's total rotation vector (the convention of the corotational elements, whose Newton
update is additive). Its work-conjugate generalised forces are

- on the node's translations: ``R F_0``;
- on the node's rotation vector: ``J(theta)^T R M_0``, with ``J`` the tangent of the
  exponential map (``delta(R) R^T = [J delta_theta x]``), since a spatial moment does work
  on the spatial spin ``J delta_theta``.

The load depends on the displacement, so Newton needs its load stiffness
``K_L = -d g / dU`` (unsymmetric: a follower load is not conservative). Both are exact:
``d(R F_0)/d theta = -[R F_0 x] J`` in closed form, and the moment term by central
differences of the closed-form ``J^T R M_0``.

The nonlinear static runner adds ``lambda g(U)`` to the external load and
``lambda K_L(U)`` to the tangent at every Newton iteration.

Public API
----------
parse_follower_load, FollowerLoads.
"""

from __future__ import annotations

import os
import re
from typing import List, Optional, Tuple

import numpy as np
import numpy.typing as npt
import scipy.sparse as sp

__all__ = ["parse_follower_load", "FollowerLoads", "FOLLOWER_LOAD_FILE"]

FOLLOWER_LOAD_FILE = "follower_load.txt"


def parse_follower_load(path: str) -> npt.NDArray[np.float64]:
    """Rows ``[x, y, z, F_x, F_y, F_z, M_x, M_y, M_z]`` under a ``[Follower load]`` header."""
    if not os.path.exists(path):
        raise FileNotFoundError(path)
    rows: List[List[float]] = []
    in_section = False
    with open(path, encoding="utf-8") as f:
        for raw in f:
            line = raw.split("#")[0].strip()
            if not line:
                continue
            if re.match(r"^\[follower[ _]?load\]$", line, re.IGNORECASE):
                in_section = True
                continue
            if not in_section or line.startswith("["):
                continue
            parts = line.split()
            if len(parts) != 9:
                raise ValueError(f"{path}: a follower load row has 9 numbers, got {line!r}")
            rows.append([float(p) for p in parts])
    return np.array(rows, dtype=np.float64).reshape(-1, 9)


def _skew(a: np.ndarray) -> np.ndarray:
    return np.array([[0.0, -a[2], a[1]], [a[2], 0.0, -a[0]], [-a[1], a[0], 0.0]])


def rotation_exp(theta: np.ndarray) -> np.ndarray:
    t = float(np.linalg.norm(theta))
    K = _skew(theta)
    if t < 1e-8:
        return np.eye(3) + K + 0.5 * K @ K
    return np.eye(3) + np.sin(t) / t * K + (1.0 - np.cos(t)) / t ** 2 * K @ K


def tangent_operator(theta: np.ndarray) -> np.ndarray:
    """``J(theta)`` with ``delta(exp(theta)) exp(theta)^T = [J delta_theta x]``."""
    t = float(np.linalg.norm(theta))
    K = _skew(theta)
    if t < 1e-6:
        return np.eye(3) + 0.5 * K + K @ K / 6.0
    return np.eye(3) + (1.0 - np.cos(t)) / t ** 2 * K + (t - np.sin(t)) / t ** 3 * K @ K


class FollowerLoads:
    """Follower loads on nodes; ``forces(U)`` and ``stiffness(U)`` on the global DOFs."""

    def __init__(self, rows: npt.ArrayLike, coordinates: npt.ArrayLike, dof_per_node: int, total_dof: int):
        rows = np.asarray(rows, dtype=np.float64).reshape(-1, 9)
        xyz = np.asarray(coordinates, dtype=np.float64)
        scale = max(float(np.ptp(xyz, axis=0).max()) if len(xyz) > 1 else 1.0, 1.0)
        self.total_dof = int(total_dof)
        self.dof_per_node = int(dof_per_node)
        self.loads: List[Tuple[int, np.ndarray, np.ndarray]] = []
        for r in rows:
            d = np.linalg.norm(xyz - r[:3], axis=1)
            n = int(np.argmin(d))
            if d[n] > 1e-6 * scale:
                raise ValueError(f"follower load at {r[:3].tolist()} is not at a node")
            self.loads.append((n, r[3:6].copy(), r[6:9].copy()))

    def __bool__(self) -> bool:
        return bool(self.loads)

    def _dofs(self, n: int) -> Tuple[np.ndarray, np.ndarray]:
        b = n * self.dof_per_node
        return np.arange(b, b + 3), np.arange(b + 3, b + 6)

    @staticmethod
    def _moment(theta: np.ndarray, M0: np.ndarray) -> np.ndarray:
        return tangent_operator(theta).T @ rotation_exp(theta) @ M0

    def forces(self, U: npt.ArrayLike) -> npt.NDArray[np.float64]:
        """Generalised follower forces (total_dof,) at displacement ``U``."""
        U = np.asarray(U, dtype=np.float64)
        g = np.zeros(self.total_dof)
        for n, F0, M0 in self.loads:
            t, r = self._dofs(n)
            th = U[r]
            g[t] += rotation_exp(th) @ F0
            if np.any(M0):
                g[r] += self._moment(th, M0)
        return g

    def stiffness(self, U: npt.ArrayLike, h: float = 1e-7) -> sp.csr_matrix:
        """Load stiffness ``K_L = -d forces / dU`` (sparse, unsymmetric)."""
        U = np.asarray(U, dtype=np.float64)
        rows, cols, vals = [], [], []

        def put(ri, ci, block):
            r, c = np.meshgrid(ri, ci, indexing="ij")
            rows.append(r.ravel())
            cols.append(c.ravel())
            vals.append(block.ravel())

        for n, F0, M0 in self.loads:
            t, r = self._dofs(n)
            th = U[r]
            f = rotation_exp(th) @ F0
            put(t, r, _skew(f) @ tangent_operator(th))          # -d(R F0)/dtheta = [R F0 x] J
            if np.any(M0):
                D = np.column_stack([(self._moment(th + e, M0) - self._moment(th - e, M0)) / (2 * h)
                                     for e in np.eye(3) * h])
                put(r, r, -D)
        if not rows:
            return sp.csr_matrix((self.total_dof, self.total_dof))
        return sp.coo_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
                             shape=(self.total_dof, self.total_dof)).tocsr()


def follower_loads_for_job(job_dir: str, coordinates, dof_per_node: int, total_dof: int) -> Optional[FollowerLoads]:
    """``FollowerLoads`` from ``job_dir/follower_load.txt``, or None."""
    path = os.path.join(job_dir, FOLLOWER_LOAD_FILE)
    if not os.path.exists(path):
        return None
    rows = parse_follower_load(path)
    return FollowerLoads(rows, coordinates, dof_per_node, total_dof) if rows.size else None

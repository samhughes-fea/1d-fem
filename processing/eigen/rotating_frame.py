# processing/eigen/rotating_frame.py
"""Free vibration of a beam structure spinning at a constant rate, in the rotating frame.

Steps (``rotating_frame_vibration``):

1. Element terms ``f_c``, ``K_soft``, ``K_rot`` and ``G`` from each element's Gauss loop
   (``pre_processing.element_library.rotating_frame``), assembled on the global DOFs.
2. Centrifugal state: ``(K + K_soft + K_rot) U_c = f_c`` on the free DOFs.
3. Centrifugal stiffening: ``K_sigma(U_c)`` from every element's
   ``linear_geometric_stiffness_matrix`` (as for linear buckling).
4. ``M q'' + G q' + K_t q = 0`` with ``K_t = K + K_sigma + K_soft + K_rot``. With Coriolis
   terms (``gyroscopic = true``) this is solved as the first-order system
   ``[[0, I], [-M^-1 K_t, -M^-1 G]]``; its eigenvalues ``s = i omega`` give the rotating-frame
   frequencies and complex modes. Without them it is the symmetric ``K_t x = omega^2 M x``.

Constrained DOFs (``bc_dofs``) are removed, not penalised: rotation terms on a support
would otherwise appear as spurious modes.

Public API
----------
RotatingFrameResult, assemble_rotating_frame, rotating_frame_vibration, campbell.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import List, Optional, Sequence

import numpy as np
import numpy.typing as npt
import scipy.linalg as sla
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from pre_processing.element_library.rotating_frame import (
    Rotation,
    RotatingFrameObject,
    element_rotating_frame,
)

__all__ = ["RotatingFrameResult", "assemble_rotating_frame", "rotating_frame_vibration", "campbell"]


@dataclass
class RotatingFrameResult:
    """Frequencies (Hz, rotating frame), modes (total_dof, n) and the centrifugal state."""

    rotation: Rotation
    frequencies_hz: npt.NDArray[np.float64]
    modes: npt.NDArray[np.complex128]
    U_centrifugal: npt.NDArray[np.float64]
    f_centrifugal: npt.NDArray[np.float64]
    element_objects: List[RotatingFrameObject] = field(default_factory=list)


def _dof_maps(elements) -> List[np.ndarray]:
    return [np.asarray(e.assemble_global_dof_indices(), dtype=np.int64).ravel() for e in elements]


def assemble_rotating_frame(elements: Sequence, coordinates: npt.ArrayLike, connectivity: Sequence,
                            rotation: Rotation, total_dof: int):
    """Global ``f_c`` (n,), ``K_soft``, ``K_rot``, ``G`` (sparse) and the element objects."""
    xyz = np.asarray(coordinates, dtype=np.float64)
    f = np.zeros(total_dof)
    rows, cols, ks, kr, gs = [], [], [], [], []
    objs = []
    for elem, dofs, (n0, n1) in zip(elements, _dof_maps(elements), connectivity):
        o = element_rotating_frame(elem, xyz[int(n0)], xyz[int(n1)], rotation, n_dof=dofs.size)
        objs.append(o)
        f[dofs] += o.f_c
        r, c = np.meshgrid(dofs, dofs, indexing="ij")
        rows.append(r.ravel())
        cols.append(c.ravel())
        ks.append(o.K_soft.ravel())
        kr.append(o.K_rot.ravel())
        gs.append(o.G.ravel())
    rows, cols = np.concatenate(rows), np.concatenate(cols)

    def build(data):
        return sp.coo_matrix((np.concatenate(data), (rows, cols)), shape=(total_dof, total_dof)).tocsr()

    return f, build(ks), build(kr), build(gs), objs


def _free(total_dof: int, bc_dofs) -> np.ndarray:
    fixed = np.zeros(total_dof, dtype=bool)
    if bc_dofs is not None and np.size(bc_dofs):
        fixed[np.asarray(bc_dofs, dtype=np.int64)] = True
    return np.flatnonzero(~fixed)


def rotating_frame_vibration(elements: Sequence, coordinates, connectivity, K: sp.spmatrix, M: sp.spmatrix,
                             bc_dofs, rotation: Rotation, num_modes: int, gyroscopic: bool = True,
                             prestress: bool = True) -> RotatingFrameResult:
    """Rotating-frame frequencies and modes (see the module docstring)."""
    from processing.buckling.linear_buckling import assemble_global_geometric_stiffness

    n = K.shape[0]
    K = sp.csr_matrix(K)
    M = sp.csr_matrix(M)
    f_c, K_soft, K_rot, G, objs = assemble_rotating_frame(elements, coordinates, connectivity, rotation, n)
    free = _free(n, bc_dofs)
    ix = np.ix_(free, free)
    U_c = np.zeros(n)
    K_t = K + K_soft + K_rot
    if prestress and np.any(f_c[free]):
        U_c[free] = spla.spsolve(sp.csc_matrix(K_t[ix]), f_c[free])
        K_t = K_t + assemble_global_geometric_stiffness(list(elements), U_c, n)
    Kf = K_t[ix].toarray()
    Mf = M[ix].toarray()
    Kf = 0.5 * (Kf + Kf.T)
    Mf = 0.5 * (Mf + Mf.T)
    m = len(free)
    if gyroscopic and rotation.omega != 0.0:
        Gf = G[ix].toarray()
        Gf = 0.5 * (Gf - Gf.T)
        # first-order form A z' = B z with z = [q, q'], s = i omega
        A = np.block([[np.eye(m), np.zeros((m, m))], [np.zeros((m, m)), Mf]])
        B = np.block([[np.zeros((m, m)), np.eye(m)], [-Kf, -Gf]])
        s, Z = sla.eig(B, A)
        om = s.imag
        keep = np.flatnonzero((om > 0.0) & np.isfinite(om))
        keep = keep[np.argsort(om[keep])][:num_modes]
        freqs = om[keep] / (2.0 * np.pi)
        X = Z[:m, keep]
    else:
        lam, X = sla.eigh(Kf, Mf)
        order = np.argsort(lam)[:num_modes]
        lam, X = lam[order], X[:, order]
        freqs = np.sign(lam) * np.sqrt(np.abs(lam)) / (2.0 * np.pi)   # negative: unstable at this speed
    modes = np.zeros((n, X.shape[1]), dtype=np.complex128)
    for j in range(X.shape[1]):
        x = X[:, j]
        k = int(np.argmax(np.abs(x)))
        modes[free, j] = x / x[k]                  # largest component real and 1
    return RotatingFrameResult(rotation, np.asarray(freqs, dtype=np.float64), modes, U_c, f_c, objs)


def campbell(elements, coordinates, connectivity, K, M, bc_dofs, rotation: Rotation, omegas: Sequence[float],
             num_modes: int, gyroscopic: bool = True, prestress: bool = True) -> npt.NDArray[np.float64]:
    """Rotating-frame frequencies (Hz) at each speed: rows ``[omega, f_1, ..., f_n]``."""
    rows = []
    for w in omegas:
        if float(w) == 0.0:
            rows.append([0.0, *_still(K, M, bc_dofs, num_modes)])
            continue
        r = rotating_frame_vibration(elements, coordinates, connectivity, K, M, bc_dofs, rotation.with_omega(w),
                                     num_modes, gyroscopic, prestress)
        f = np.full(num_modes, np.nan)
        f[:r.frequencies_hz.size] = r.frequencies_hz
        rows.append([float(w), *f])
    return np.array(rows)


def _still(K, M, bc_dofs, num_modes: int) -> npt.NDArray[np.float64]:
    n = K.shape[0]
    free = _free(n, bc_dofs)
    ix = np.ix_(free, free)
    lam = sla.eigh(sp.csr_matrix(K)[ix].toarray(), sp.csr_matrix(M)[ix].toarray(), eigvals_only=True)
    lam = np.sort(lam)[:num_modes]
    out = np.full(num_modes, np.nan)
    out[:lam.size] = np.sign(lam) * np.sqrt(np.abs(lam)) / (2.0 * np.pi)
    return out

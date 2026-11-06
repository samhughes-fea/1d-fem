# pre_processing/element_library/linear/beam/fosdt/gbt/utilities/B_matrix.py
"""Strain-displacement **B** for GBT (variable rows/cols) — infinitesimal strain.

Rows 0–5 reuse Timoshenko **B** (6×12). Extra rows are ``b'``, ``b``, ``a''``, ``a'``, ``a``.
Assumed-strain (B-bar) projection on reduced shear/warping rows is
``reduced_projection``.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Sequence

import numpy as np

from onedfem.pre_processing.element_library.linear.beam.fosdt.gbt.utilities.shape_functions import hermite_x

_N_STD = 12


@dataclass(frozen=True)
class StrainDisplacementOperator:
    """GBT **B** at a station, with optional B-bar on reduced rows."""

    element_length: float
    n_dof: int
    n_modes: int
    n_warping: int
    parent: object
    warping_cols: Callable[[int], tuple]
    mode_cols: Callable[[int], tuple]
    reduced_rows: np.ndarray
    natural_shape: Callable[[np.ndarray], tuple]

    def raw_at(self, xi_g: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        N, dN, d2N = self.natural_shape(np.array([xi_g]))
        B6 = self.parent.physical_coordinate_form(dN, d2N, N)[0]
        m, p = self.n_modes, self.n_warping
        B = np.zeros((6 + 2 * p + 3 * m, self.n_dof))
        B[:6, :_N_STD] = B6
        L1, L2 = 0.5 * (1.0 - xi_g), 0.5 * (1.0 + xi_g)
        L = float(self.element_length)
        for j in range(p):
            c = list(self.warping_cols(j))
            B[6 + j, c] = [-1.0 / L, 1.0 / L]
            B[6 + p + j, c] = [L1, L2]
        H, dH, d2H = hermite_x(float(xi_g), L)
        oa = 6 + 2 * p
        for k in range(m):
            c = list(self.mode_cols(k))
            B[oa + k, c] = d2H
            B[oa + m + k, c] = dH
            B[oa + 2 * m + k, c] = H
        return B, N, dN

    def reduced_projection(self, xi: float, shear_block: int) -> np.ndarray:
        q = int(shear_block)
        r = self.reduced_rows
        xs, ws = np.polynomial.legendre.leggauss(q)
        Br = np.array([self.raw_at(float(x))[0][r] for x in xs])
        out = np.zeros_like(Br[0])
        for k in range(q):
            Pk = np.polynomial.legendre.Legendre.basis(k)
            ck = 0.5 * (2 * k + 1) * np.einsum("g,gij->ij", ws * Pk(xs), Br)
            out += ck * Pk(xi)
        return out

    def at(self, xi_g: float, shear_block: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        B, N, dN = self.raw_at(xi_g)
        B[self.reduced_rows] = self.reduced_projection(float(xi_g), shear_block)
        return B, N, dN

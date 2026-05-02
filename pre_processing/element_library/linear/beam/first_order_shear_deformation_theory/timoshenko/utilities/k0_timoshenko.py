# pre_processing/element_library/linear/beam/first_order_shear_deformation_theory/timoshenko/utilities/k0_timoshenko.py
"""Shared Timoshenko material stiffness assembly for linear ``K_e`` and nonlinear TL ``K_0``."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np


@dataclass(frozen=True)
class TimoshenkoQuadratureOrders:
    """Gauss–Legendre orders for Timoshenko ``Bᵀ D B`` stiffness assembly."""

    axial: int
    bending_y: int
    bending_z: int
    shear_y: int
    shear_z: int
    shear_block: int
    torsion: int

    @property
    def max_full_order(self) -> int:
        return max(
            self.axial,
            self.bending_y,
            self.bending_z,
            self.shear_y,
            self.shear_z,
            self.torsion,
        )

    @property
    def bending_rule_order(self) -> int:
        return max(self.bending_y, self.bending_z)

    @property
    def loop_order(self) -> int:
        """Single Gauss order for TL ``F_int``, ``K_sigma``, and distributed-load loops."""
        return max(self.max_full_order, 2)


def timoshenko_quadrature_orders_from_element_array(
    element_array: np.ndarray,
    *,
    shear_block_order: int | None = None,
) -> TimoshenkoQuadratureOrders:
    """
    Build orders from the per-element integration slice (indices 3–8).

    Parameters
    ----------
    element_array
        Row layout: ``[id, n1, n2, axial, bending_y, bending_z, shear_y, shear_z, torsion, load]``.
    shear_block_order
        Quadrature order for the shear-stiffness block (strain rows 3–4). If ``None``,
        uses ``1`` (legacy reduced integration). Mesh ``shear_y`` / ``shear_z`` columns
        still enter ``max_full_order`` for the baseline full rule.
    """
    axial_order = max(int(element_array[3]), 1)
    bending_y_order = max(int(element_array[4]), 2)
    bending_z_order = max(int(element_array[5]), 2)
    shear_y_order = max(int(element_array[6]), 2) if element_array[6] > 0 else 2
    shear_z_order = max(int(element_array[7]), 2) if element_array[7] > 0 else 2
    torsion_order = max(int(element_array[8]), 1)
    shear_block = 1 if shear_block_order is None else max(int(shear_block_order), 1)
    return TimoshenkoQuadratureOrders(
        axial=axial_order,
        bending_y=bending_y_order,
        bending_z=bending_z_order,
        shear_y=shear_y_order,
        shear_z=shear_z_order,
        shear_block=shear_block,
        torsion=torsion_order,
    )


def assemble_timoshenko_K0(
    shape_function_operator: Any,
    strain_displacement_operator: Any,
    D: np.ndarray,
    detJ: float,
    orders: TimoshenkoQuadratureOrders,
) -> np.ndarray:
    """
    Assemble the 12×12 material stiffness ∫ ``Bᵀ D B`` dx with block-wise quadrature.

    A full rule at ``max_full_order`` forms the baseline (captures coupling); bending
    rows 1–2 and shear rows 3–4 are then substituted by separate Gauss sums at
    ``bending_rule_order`` and ``orders.shear_block``.
    """
    max_order = orders.max_full_order
    bending_order = orders.bending_rule_order
    shear_order = orders.shear_block

    xi_full, w_full = np.polynomial.legendre.leggauss(max_order)
    Ke_full = np.zeros((12, 12), dtype=np.float64)
    for xi_g, w_g in zip(xi_full, w_full):
        N, dN_dξ, d2N_dξ2 = shape_function_operator.natural_coordinate_form(np.array([xi_g]))
        B = strain_displacement_operator.physical_coordinate_form(dN_dξ, d2N_dξ2, N)[0]
        Ke_full += B.T @ D @ B * w_g * detJ

    xi_bending, w_bending = np.polynomial.legendre.leggauss(bending_order)
    Ke_bending_block = np.zeros((12, 12), dtype=np.float64)
    for xi_g, w_g in zip(xi_bending, w_bending):
        N, dN_dξ, d2N_dξ2 = shape_function_operator.natural_coordinate_form(np.array([xi_g]))
        B = strain_displacement_operator.physical_coordinate_form(dN_dξ, d2N_dξ2, N)[0]
        B_bending = B[[1, 2], :]
        D_bending_diag = np.array([D[1, 1], D[2, 2]])
        Ke_bending_block += B_bending.T @ np.diag(D_bending_diag) @ B_bending * w_g * detJ

    xi_shear, w_shear = np.polynomial.legendre.leggauss(shear_order)
    Ke_shear_block = np.zeros((12, 12), dtype=np.float64)
    for xi_g, w_g in zip(xi_shear, w_shear):
        N, dN_dξ, d2N_dξ2 = shape_function_operator.natural_coordinate_form(np.array([xi_g]))
        B = strain_displacement_operator.physical_coordinate_form(dN_dξ, d2N_dξ2, N)[0]
        B_shear = B[[3, 4], :]
        D_shear_diag = np.array([D[3, 3], D[4, 4]])
        Ke_shear_block += B_shear.T @ np.diag(D_shear_diag) @ B_shear * w_g * detJ

    Ke_full_bending = np.zeros((12, 12), dtype=np.float64)
    for xi_g, w_g in zip(xi_full, w_full):
        N, dN_dξ, d2N_dξ2 = shape_function_operator.natural_coordinate_form(np.array([xi_g]))
        B = strain_displacement_operator.physical_coordinate_form(dN_dξ, d2N_dξ2, N)[0]
        B_bending = B[[1, 2], :]
        D_bending_diag = np.array([D[1, 1], D[2, 2]])
        Ke_full_bending += B_bending.T @ np.diag(D_bending_diag) @ B_bending * w_g * detJ

    Ke_full_shear = np.zeros((12, 12), dtype=np.float64)
    for xi_g, w_g in zip(xi_full, w_full):
        N, dN_dξ, d2N_dξ2 = shape_function_operator.natural_coordinate_form(np.array([xi_g]))
        B = strain_displacement_operator.physical_coordinate_form(dN_dξ, d2N_dξ2, N)[0]
        B_shear = B[[3, 4], :]
        D_shear_diag = np.array([D[3, 3], D[4, 4]])
        Ke_full_shear += B_shear.T @ np.diag(D_shear_diag) @ B_shear * w_g * detJ

    return Ke_full - Ke_full_bending - Ke_full_shear + Ke_bending_block + Ke_shear_block

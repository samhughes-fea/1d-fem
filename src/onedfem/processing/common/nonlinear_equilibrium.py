from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Optional

import numpy as np


@dataclass
class NonlinearEquilibriumIterationRecord:
    load_increment_index: int
    load_factor: float
    newton_iter: int
    norm_R_full: float
    norm_F_cond: float
    threshold: float
    norm_delta_u: float
    alpha: float
    residual_ok: bool


@dataclass
class NonlinearEquilibriumResult:
    U_global: np.ndarray
    delta_U_cond: Optional[np.ndarray]
    converged: bool
    iterations_used: int
    last_norm_F_cond: Optional[float]


@dataclass
class ArcLengthCorrectorResult:
    U_global: np.ndarray
    load_factor: float
    converged: bool
    iterations_used: int
    last_norm_F_cond: Optional[float]
    predictor_norm: float


def newton_condensed_residual_converged(
    norm_r_cond: float,
    atol: float,
    rtol: float | None,
    ref_scale: float,
) -> bool:
    rt = 0.0 if rtol is None else float(rtol)
    ref = max(float(ref_scale), np.finfo(float).eps)
    return bool(float(norm_r_cond) <= float(atol) + rt * ref)


def solve_newton_equilibrium(
    *,
    U_global: np.ndarray,
    F_ext_global: np.ndarray,
    load_factor: float,
    load_increment_index: int,
    newton_tol: float,
    newton_max_iter: int,
    newton_tol_delta_u: float,
    newton_relative_tol: float | None,
    newton_relative_reference: str,
    build_system: Callable[[np.ndarray], tuple[float, float, np.ndarray]],
    solve_condensed_step: Callable[[int, float], np.ndarray],
    reconstruct_delta: Callable[[np.ndarray], np.ndarray],
    condensed_residual_norm: Callable[[np.ndarray, np.ndarray], float],
    iteration_callback: Callable[[NonlinearEquilibriumIterationRecord], None],
    line_search_enabled: bool = False,
    line_search_shrink: float = 0.5,
    line_search_max_backtracks: int = 6,
) -> NonlinearEquilibriumResult:
    delta_U_cond = None
    newton_ref_residual = None
    newton_ref_external = None
    converged = False
    last_norm_F_cond: Optional[float] = None
    iterations_used = 0

    for iteration in range(newton_max_iter):
        norm_R_full, norm_R_conv, _ = build_system(U_global)
        last_norm_F_cond = norm_R_conv
        if iteration == 0:
            newton_ref_residual = norm_R_conv
            newton_ref_external = float(np.linalg.norm(np.asarray(F_ext_global, dtype=np.float64).ravel()))

        ref_scale = (
            newton_ref_external if newton_relative_reference == "external_force" else newton_ref_residual
        )
        residual_ok = newton_condensed_residual_converged(
            norm_R_conv,
            newton_tol,
            newton_relative_tol,
            ref_scale,
        )

        delta_U_cond = solve_condensed_step(iteration + 1, float(load_factor))
        delta_U_global = reconstruct_delta(delta_U_cond)

        alpha = 1.0
        if line_search_enabled:
            best_norm = float("inf")
            best_alpha = 1.0
            for k in range(line_search_max_backtracks + 1):
                a = float(line_search_shrink ** k)
                if a < 1e-14:
                    break
                U_try = U_global + a * delta_U_global
                n_try = condensed_residual_norm(U_try, F_ext_global)
                if n_try < best_norm:
                    best_norm = n_try
                    best_alpha = a
            alpha = best_alpha

        U_global = U_global + alpha * delta_U_global
        norm_du = float(np.linalg.norm(alpha * delta_U_global))
        threshold = float(newton_tol)
        if newton_relative_tol is not None:
            rs = max(float(ref_scale), np.finfo(float).eps)
            threshold = float(newton_tol) + float(newton_relative_tol) * rs

        iterations_used = iteration + 1
        iteration_callback(
            NonlinearEquilibriumIterationRecord(
                load_increment_index=load_increment_index,
                load_factor=float(load_factor),
                newton_iter=iteration + 1,
                norm_R_full=float(norm_R_full),
                norm_F_cond=float(norm_R_conv),
                threshold=float(threshold),
                norm_delta_u=float(norm_du),
                alpha=float(alpha),
                residual_ok=bool(residual_ok),
            )
        )
        if residual_ok and norm_du < float(newton_tol_delta_u):
            converged = True
            break

    return NonlinearEquilibriumResult(
        U_global=U_global,
        delta_U_cond=delta_U_cond,
        converged=converged,
        iterations_used=iterations_used,
        last_norm_F_cond=last_norm_F_cond,
    )


def arc_length_predictor_corrector_step(
    *,
    U_prev: np.ndarray,
    load_factor_prev: float,
    predictor_displacement: np.ndarray,
    reference_load_vector: np.ndarray,
    arc_length_radius: float,
    alpha_scale: float = 1.0,
) -> tuple[np.ndarray, float]:
    """Return a simple arc-length predictor using a quadratic constraint.

    The predictor direction is provided by ``predictor_displacement`` and the load
    increment is chosen so that

    ``||ΔU||² + (alpha_scale * Δλ)² ||F_ref||² = arc_length_radius²``.
    """
    du = np.asarray(predictor_displacement, dtype=np.float64).ravel()
    fref = np.asarray(reference_load_vector, dtype=np.float64).ravel()
    du_norm_sq = float(np.dot(du, du))
    f_norm_sq = float(np.dot(fref, fref))
    scaled_f_sq = max(float(alpha_scale) ** 2 * f_norm_sq, np.finfo(float).eps)
    remaining = float(arc_length_radius) ** 2 - du_norm_sq
    delta_lambda = 0.0 if remaining <= 0.0 else float(np.sqrt(remaining / scaled_f_sq))
    return U_prev + du, float(load_factor_prev + delta_lambda)


def solve_arc_length_corrector(
    *,
    U_prev: np.ndarray,
    load_factor_prev: float,
    predictor_displacement: np.ndarray,
    reference_load_vector: np.ndarray,
    arc_length_radius: float,
    alpha_scale: float,
    newton_tol: float,
    newton_max_iter: int,
    build_system_from_state: Callable[[np.ndarray, float], tuple[float, float, np.ndarray]],
    solve_condensed_step_from_state: Callable[[int, float], np.ndarray],
    reconstruct_delta: Callable[[np.ndarray], np.ndarray],
    iteration_callback: Callable[[NonlinearEquilibriumIterationRecord], None],
    load_increment_index: int,
) -> ArcLengthCorrectorResult:
    U_trial, load_factor = arc_length_predictor_corrector_step(
        U_prev=U_prev,
        load_factor_prev=load_factor_prev,
        predictor_displacement=predictor_displacement,
        reference_load_vector=reference_load_vector,
        arc_length_radius=arc_length_radius,
        alpha_scale=alpha_scale,
    )
    predictor_norm = float(np.linalg.norm(np.asarray(U_trial - U_prev, dtype=np.float64).ravel()))
    converged = False
    last_norm: Optional[float] = None
    iterations_used = 0
    fref = np.asarray(reference_load_vector, dtype=np.float64).ravel()
    scaled_f_sq = max((float(alpha_scale) ** 2) * float(np.dot(fref, fref)), np.finfo(float).eps)

    for iteration in range(newton_max_iter):
        norm_R_full, norm_R_conv, _ = build_system_from_state(U_trial, float(load_factor))
        last_norm = norm_R_conv
        delta_u_cond = solve_condensed_step_from_state(iteration + 1, float(load_factor))
        delta_u_global = reconstruct_delta(delta_u_cond)

        du_from_prev = np.asarray(U_trial - U_prev, dtype=np.float64).ravel()
        correction = np.asarray(delta_u_global, dtype=np.float64).ravel()
        g = float(np.dot(du_from_prev, du_from_prev) + scaled_f_sq * ((float(load_factor) - float(load_factor_prev)) ** 2) - float(arc_length_radius) ** 2)
        denom = 2.0 * scaled_f_sq * max(abs(float(load_factor) - float(load_factor_prev)), np.sqrt(np.finfo(float).eps))
        delta_lambda = -g / denom
        U_trial = U_trial + correction
        load_factor = float(load_factor + delta_lambda)
        norm_du = float(np.linalg.norm(correction))
        threshold = float(newton_tol)
        iterations_used = iteration + 1
        iteration_callback(
            NonlinearEquilibriumIterationRecord(
                load_increment_index=load_increment_index,
                load_factor=float(load_factor),
                newton_iter=iteration + 1,
                norm_R_full=float(norm_R_full),
                norm_F_cond=float(norm_R_conv),
                threshold=threshold,
                norm_delta_u=float(norm_du),
                alpha=1.0,
                residual_ok=bool(norm_R_conv <= threshold and abs(g) <= max(threshold, 1e-10)),
            )
        )
        if norm_R_conv <= threshold and abs(g) <= max(threshold, 1e-10):
            converged = True
            break

    return ArcLengthCorrectorResult(
        U_global=np.asarray(U_trial, dtype=np.float64),
        load_factor=float(load_factor),
        converged=converged,
        iterations_used=iterations_used,
        last_norm_F_cond=last_norm,
        predictor_norm=predictor_norm,
    )

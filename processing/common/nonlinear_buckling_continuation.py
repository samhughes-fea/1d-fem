from __future__ import annotations

from dataclasses import dataclass
import csv
import json
import os
from typing import Optional

import numpy as np

from processing.common.nonlinear_equilibrium import arc_length_predictor_corrector_step


@dataclass
class ContinuationConfig:
    continuation_method: str
    load_factors: np.ndarray
    arc_length_radius: float = 1.0
    arc_length_alpha_scale: float = 1.0


@dataclass
class ArcLengthStepResult:
    load_factor: float
    predictor_load_factor: float
    converged: bool
    iterations_used: int
    residual_norm: float | None
    tip_displacement: float


@dataclass
class ImperfectionConfig:
    source: Optional[str] = None
    mode_index: int = 0
    scale: float = 0.0


def seed_initial_imperfection(
    *,
    job_name: str,
    primary_results_dir: str,
    job_results_dir: str,
    U_global: np.ndarray,
    config: ImperfectionConfig,
) -> tuple[np.ndarray, dict | None]:
    if config.source is None or abs(float(config.scale)) <= 0.0:
        return U_global, None
    source = str(config.source).strip().lower()
    if source != "linear_buckling":
        raise ValueError("buckling.imperfection_source must be 'linear_buckling' when set")

    results_dir = os.path.join(primary_results_dir, "modal_results")
    modes_path = os.path.join(results_dir, f"{job_name}_buckling_mode_shapes.txt")
    if not os.path.isfile(modes_path):
        raise FileNotFoundError(
            f"Imperfection seeding requested but linear buckling mode file not found: {modes_path}"
        )
    mode_shapes = np.loadtxt(modes_path, dtype=np.float64)
    if mode_shapes.ndim == 1:
        mode_shapes = mode_shapes.reshape(-1, 1)
    if config.mode_index < 0 or config.mode_index >= mode_shapes.shape[1]:
        raise ValueError(
            f"buckling.imperfection_mode_index={config.mode_index} out of range 0..{mode_shapes.shape[1]-1}"
        )
    phi = np.asarray(mode_shapes[:, config.mode_index], dtype=np.float64).ravel()
    if phi.shape[0] != U_global.shape[0]:
        raise ValueError(f"Imperfection mode length {phi.shape[0]} does not match total DOF {U_global.shape[0]}")
    seeded = np.asarray(U_global, dtype=np.float64).copy() + float(config.scale) * phi
    meta = {
        "source": source,
        "mode_index": int(config.mode_index),
        "scale": float(config.scale),
        "mode_shapes_path": os.path.relpath(modes_path, job_results_dir).replace("\\", "/"),
    }
    return seeded, meta


def predictor_load_factor_for_increment(
    *,
    config: ContinuationConfig,
    increment_index: int,
    current_U: np.ndarray,
    tip_dof: int,
    reference_load_vector: np.ndarray,
) -> float:
    lam = float(config.load_factors[increment_index])
    if config.continuation_method != "arc_length":
        return lam
    predictor_du = np.zeros_like(current_U)
    if current_U.size:
        predictor_du[tip_dof] = float(config.arc_length_radius)
    _, predictor_load_factor = arc_length_predictor_corrector_step(
        U_prev=current_U,
        load_factor_prev=0.0 if increment_index == 0 else float(config.load_factors[increment_index - 1]),
        predictor_displacement=predictor_du,
        reference_load_vector=reference_load_vector,
        arc_length_radius=float(config.arc_length_radius),
        alpha_scale=float(config.arc_length_alpha_scale),
    )
    return float(predictor_load_factor)


def write_continuation_history_header(path: str) -> None:
    with open(path, "w", newline="", encoding="utf-8") as f:
        csv.writer(f).writerow(
            [
                "increment_index",
                "load_factor",
                "converged",
                "iterations_used",
                "residual_norm",
                "tip_dof",
                "tip_displacement",
                "continuation_method",
                "predictor_load_factor",
            ]
        )


def append_continuation_history_row(
    *,
    path: str,
    increment_index: int,
    load_factor: float,
    converged: bool,
    iterations_used: int,
    residual_norm: float | None,
    tip_dof: int,
    tip_displacement: float,
    continuation_method: str,
    predictor_load_factor: float,
) -> None:
    with open(path, "a", newline="", encoding="utf-8") as f:
        csv.writer(f).writerow(
            [
                increment_index,
                float(load_factor),
                int(bool(converged)),
                int(iterations_used),
                "" if residual_norm is None else float(residual_norm),
                int(tip_dof),
                float(tip_displacement),
                continuation_method,
                float(predictor_load_factor),
            ]
        )


def write_continuation_summary(path: str, payload: dict) -> None:
    with open(path, "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2)

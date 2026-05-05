from __future__ import annotations

import csv
import os

import numpy as np


def write_transient_primary_summary(
    *,
    primary_results_dir: str,
    job_name: str,
    t_grid: np.ndarray,
    U: np.ndarray,
    V: np.ndarray,
    A: np.ndarray,
    damping_source: str,
    n_bc_dofs: int,
) -> str:
    path = os.path.join(primary_results_dir, "dynamic_results", f"{job_name}_primary_summary.csv")
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["metric", "value"])
        w.writerow(["num_time_steps", int(max(0, np.asarray(t_grid).shape[0] - 1))])
        w.writerow(["num_time_samples", int(np.asarray(t_grid).shape[0])])
        w.writerow(["total_dof", int(np.asarray(U).shape[1] if np.asarray(U).ndim == 2 else 0)])
        w.writerow(["max_abs_displacement", float(np.max(np.abs(U)))])
        w.writerow(["max_abs_velocity", float(np.max(np.abs(V)))])
        w.writerow(["max_abs_acceleration", float(np.max(np.abs(A)))])
        w.writerow(["damping_source", damping_source])
        w.writerow(["num_bc_dofs", int(n_bc_dofs)])
    return path

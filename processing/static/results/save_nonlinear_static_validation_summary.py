from __future__ import annotations

import csv
import os

import numpy as np


def write_nonlinear_static_tip_history_summary(
    *,
    primary_results_dir: str,
    job_name: str,
    load_factors: np.ndarray,
    tip_displacements: np.ndarray,
) -> str:
    out_dir = os.path.join(primary_results_dir, "nonlinear_static_validation")
    os.makedirs(out_dir, exist_ok=True)
    path = os.path.join(out_dir, f"{job_name}_tip_load_history.csv")
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["load_factor", "tip_displacement"])
        for lam, u in zip(np.asarray(load_factors, dtype=float).ravel(), np.asarray(tip_displacements, dtype=float).ravel()):
            w.writerow([float(lam), float(u)])
    return path

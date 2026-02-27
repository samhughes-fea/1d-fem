#!/usr/bin/env python3
"""
Create mesh variants for Timoshenko jobs (job_0006–job_0011).
Same geometry, loads, and mesh densities (n4, n8, n16, n32, n64, n128) as jobs 0–5,
but element_type = TimoshenkoBeamElement3D.
  - job_0006, 0007, 0008: point load at end, midspan, quarter (same as 0, 1, 2).
  - job_0009, 0010, 0011: UDL, triangular, parabolic (same as 3, 4, 5).
Run from repo root: python pre_processing/mesh_library/create_timoshenko_mesh_variants.py
"""
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from pre_processing.mesh_library.schemes.mesh_generator import (
    generate_mesh,
    save_grid_file,
    save_element_file,
    save_material_file,
    save_section_file,
)

JOBS_DIR = REPO_ROOT / "jobs"
L = 2.0  # m (must match mesh_generator default)
P = -500.0  # N (Fy at load point)
w = 500.0  # N/m

VARIANT_NS = [4, 8, 16, 32, 64, 128, 500]  # n500 = Abaqus validation reference (converged)
# Point load: (base_id, label, formula, x_position)
POINT_JOBS = [
    (6, "End load", "P(x=L)", L),
    (7, "Midspan load", "P(x=L/2)", L / 2),
    (8, "Quarter-point load", "P(x=L/4)", L / 4),
]
# Distributed: (base_id, load_type, fy_fn)
DIST_JOBS = [
    (9, "UDL", lambda x: -w),
    (10, "TRIANGULAR", lambda x: -w * (x / L) if L else 0),
    (11, "PARABOLIC", lambda x: -w * (x / L) ** 2 if L else 0),
]
LOAD_FORMULAS = {"UDL": "q(x) = w", "TRIANGULAR": "q(x) = w * (x/L)", "PARABOLIC": "q(x) = w * (x/L)^2"}


def write_point_load(dir_path: Path, label: str, formula: str, x_load: float) -> None:
    with open(dir_path / "point_load.txt", "w", encoding="utf-8") as f:
        f.write(f"# Load Type: {label}\n")
        f.write(f"# Load Formula: {formula}\n")
        f.write("[Point load]\n")
        f.write("         [x]          [y]          [z]        [F_x]        [F_y]        [F_z]        [M_x]        [M_y]        [M_z]\n")
        f.write(f"    {x_load:.6f}     0.000000     0.000000     0.000000  {P:.6f}     0.000000     0.000000     0.000000     0.000000\n")


def write_distributed_load(dir_path: Path, n: int, load_type: str, fy_at_x) -> None:
    with open(dir_path / "distributed_load.txt", "w", encoding="utf-8") as f:
        f.write(f"# Load Type: {load_type}\n# Load Formula: {LOAD_FORMULAS[load_type]}\n[Distributed load]\n")
        f.write("         [x]          [y]          [z]        [F_x]        [F_y]        [F_z]        [M_x]        [M_y]        [M_z]\n")
        for i in range(n + 1):
            x = L * i / n
            fy = fy_at_x(x)
            f.write(f"    {x:.6f}     0.000000     0.000000     0.000000  {fy:.6f}     0.000000     0.000000     0.000000     0.000000\n")


def write_prescribed_displacement(dir_path: Path) -> None:
    with open(dir_path / "prescribed_displacement.txt", "w", encoding="utf-8") as f:
        f.write("[Prescribed Displacement]\n[id]     [node_id]  [dof]   [value]     [type]          [comment]\n")
        for i, dof in enumerate(["UX", "UY", "UZ", "RX", "RY", "RZ"]):
            f.write(f"{i}        0          {dof}      0.0         displacement     # Fixed support\n")


def write_simulation_settings(dir_path: Path) -> None:
    with open(dir_path / "simulation_settings.txt", "w", encoding="utf-8") as f:
        f.write("[Simulation]\n[Type]      \nStatic             \n#Dynamic          \n#Modal\n")


def main() -> None:
    # Point load (6, 7, 8)
    for base_id, label, formula, x_load in POINT_JOBS:
        for n in VARIANT_NS:
            name = f"job_{base_id:04d}_n{n}"
            dir_path = JOBS_DIR / name
            dir_path.mkdir(exist_ok=True)
            save_dir = str(dir_path)
            node_positions, elements = generate_mesh(
                growth=0, max_nodes=n + 1, uniform_nodes=n + 1, length=L
            )
            save_grid_file(node_positions, save_dir)
            save_element_file(
                elements, save_dir, element_type="TimoshenkoBeamElement3D"
            )
            save_material_file(elements, save_dir)
            save_section_file(elements, save_dir)
            write_point_load(dir_path, label, formula, x_load)
            write_prescribed_displacement(dir_path)
            write_simulation_settings(dir_path)
            print(f"Created {name}")

    # Distributed (9, 10, 11)
    for base_id, load_type, fy_fn in DIST_JOBS:
        for n in VARIANT_NS:
            name = f"job_{base_id:04d}_n{n}"
            dir_path = JOBS_DIR / name
            dir_path.mkdir(exist_ok=True)
            save_dir = str(dir_path)
            node_positions, elements = generate_mesh(
                growth=0, max_nodes=n + 1, uniform_nodes=n + 1, length=L
            )
            save_grid_file(node_positions, save_dir)
            save_element_file(
                elements, save_dir, element_type="TimoshenkoBeamElement3D"
            )
            save_material_file(elements, save_dir)
            save_section_file(elements, save_dir)
            write_distributed_load(dir_path, n, load_type, fy_fn)
            write_prescribed_displacement(dir_path)
            write_simulation_settings(dir_path)
            print(f"Created {name}")


if __name__ == "__main__":
    main()

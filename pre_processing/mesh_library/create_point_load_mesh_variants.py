#!/usr/bin/env python3
"""
Create mesh variants for point-load jobs (job_0000, job_0001, job_0002).
Generates job_XXXX_n4, n8, n16, n32, n64, n128 using the mesh library
(pre_processing/mesh_library/schemes/mesh_generator.py) for geometry and properties,
then adds job-specific point load, BCs, and simulation files.
Run from repo root: python pre_processing/mesh_library/create_point_load_mesh_variants.py
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

VARIANT_NS = [4, 8, 16, 32, 64, 128, 500]  # n500 = Abaqus validation reference (converged)
# (base_id, load_label, load_formula, x_position)
BASE_JOBS = [
    (0, "End load", "P(x=L)", L),
    (1, "Midspan load", "P(x=L/2)", L / 2),
    (2, "Quarter-point load", "P(x=L/4)", L / 4),
]


def write_point_load(dir_path: Path, base_id: int, x: float) -> None:
    labels = ["End load", "Midspan load", "Quarter-point load"]
    formulas = ["P(x=L)", "P(x=L/2)", "P(x=L/4)"]
    label = labels[base_id]
    formula = formulas[base_id]
    with open(dir_path / "point_load.txt", "w", encoding="utf-8") as f:
        f.write(f"# Load Type: {label}\n")
        f.write(f"# Load Formula: {formula}\n")
        f.write("[Point load]\n")
        f.write("         [x]          [y]          [z]        [F_x]        [F_y]        [F_z]        [M_x]        [M_y]        [M_z]\n")
        f.write(f"    {x:.6f}     0.000000     0.000000     0.000000  {P:.6f}     0.000000     0.000000     0.000000     0.000000\n")


def write_prescribed_displacement(dir_path: Path) -> None:
    with open(dir_path / "prescribed_displacement.txt", "w", encoding="utf-8") as f:
        f.write("[Prescribed Displacement]\n")
        f.write("[id]     [node_id]  [dof]   [value]     [type]          [comment]\n")
        for i, dof in enumerate(["UX", "UY", "UZ", "RX", "RY", "RZ"]):
            f.write(f"{i}        0          {dof}      0.0         displacement     # Fixed support\n")


def write_simulation_settings(dir_path: Path) -> None:
    with open(dir_path / "simulation_settings.txt", "w", encoding="utf-8") as f:
        f.write("[Simulation]\n[Type]      \nStatic             \n#Dynamic          \n#Modal\n")


def main() -> None:
    for base_id, load_label, _formula, x_load in BASE_JOBS:
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
                elements, save_dir, element_type="EulerBernoulliBeamElement3D"
            )
            save_material_file(elements, save_dir)
            save_section_file(elements, save_dir)
            write_point_load(dir_path, base_id, x_load)
            write_prescribed_displacement(dir_path)
            write_simulation_settings(dir_path)
            print(f"Created {name}")


if __name__ == "__main__":
    main()

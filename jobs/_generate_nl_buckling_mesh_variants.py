"""
Generate mesh-refinement job variants for the nonlinear buckling imperfect-column benchmark.
Creates job_benchmark_nl_buckling_imperfect_column_n<N> directories for N in N_VALUES.
"""
from __future__ import annotations

import sys
from pathlib import Path

JOBS_DIR = Path(__file__).resolve().parent
N_VALUES = [4, 8, 16, 32, 64, 128, 500]
L = 2.0  # beam length [m]

BASE_JOB = "job_benchmark_nl_buckling_imperfect_column"

ELEMENT_TYPE = "NonlinearEulerBernoulliBeamElement3D"
ELEMENT_ROW = (
    "  {eid:<5} {n1:<5} {n2:<5} "
    "{etype:<40} 3                 3                 3                 "
    "0                 0                 3                 2"
)
MATERIAL_ROW = "{eid:<5} 2.1e+11     8.1000e+10   0.3    7850"
SECTION_ROW = "{eid}  0.00131  0.0  3.23400e-07  2.08769e-06  2.60673e-08"

PRESCRIBED_DISP = """\
[Prescribed Displacement]
[id]     [node_id]  [dof]   [value]     [type]          [comment]
0        0          UX      0.0         displacement     # Fixed support
1        0          UY      0.0         displacement     # Fixed support
2        0          UZ      0.0         displacement     # Fixed support
3        0          RX      0.0         displacement     # Fixed support
4        0          RY      0.0         displacement     # Fixed support
5        0          RZ      0.0         displacement     # Fixed support
"""

POINT_LOAD = """\
# Axial compression benchmark for nonlinear imperfect-column continuation.
[Point load]
         [x]          [y]          [z]        [F_x]        [F_y]        [F_z]        [M_x]        [M_y]        [M_z]
    2.000000     0.000000     0.000000 -1000.000000     0.000000     0.000000     0.000000     0.000000     0.000000
"""

SIM_SETTINGS = """\
[Simulation]
[Type]
buckling

[Buckling]
enabled = true
nonlinear_buckling = true
continuation_method = arc_length
arc_length_radius = 0.01
arc_length_alpha_scale = 1.0
num_increments = 3
imperfection_source = linear_buckling
imperfection_mode_index = 0
imperfection_scale = 1e-4

[Newton]
tolerance = 1e-4
max_iterations = 12
tolerance_delta_u = 1e-9

[Solver]
tolerance = 1e-10
max_iterations = 2000
"""


def _grid_txt(n: int) -> str:
    dx = L / n
    lines = ["[Grid]", "[node_id]   [x]         [y]       [z]"]
    for i in range(n + 1):
        x = round(i * dx, 10)
        lines.append(f"{i:<12}{x:<12.6f}0.0       0.0")
    return "\n".join(lines) + "\n"


def _element_txt(n: int) -> str:
    header = (
        "[Element]\n"
        "[element_id]  [node1]  [node2]  [element_type]                "
        "[axial_order]     [bending_y_order] [bending_z_order] [shear_y_order]   "
        "[shear_z_order]   [torsion_order]   [load_order]"
    )
    rows = [header]
    for i in range(n):
        rows.append(ELEMENT_ROW.format(eid=i, n1=i, n2=i + 1, etype=ELEMENT_TYPE))
    return "\n".join(rows) + "\n"


def _material_txt(n: int) -> str:
    lines = ["[Material]", "[element_id]  [E]         [G]          [nu]   [rho]"]
    for i in range(n):
        lines.append(MATERIAL_ROW.format(eid=i))
    return "\n".join(lines) + "\n"


def _section_txt(n: int) -> str:
    lines = ["[Section]", "[element_id]  [A]  [I_x]  [I_y]  [I_z]  [J_t]"]
    for i in range(n):
        lines.append(SECTION_ROW.format(eid=i))
    return "\n".join(lines) + "\n"


def generate(dry_run: bool = False) -> None:
    created = 0
    for n in N_VALUES:
        dir_name = f"{BASE_JOB}_n{n}"
        job_dir = JOBS_DIR / dir_name
        if not dry_run:
            job_dir.mkdir(exist_ok=True)

        files = {
            "grid.txt": _grid_txt(n),
            "element.txt": _element_txt(n),
            "material.txt": _material_txt(n),
            "section.txt": _section_txt(n),
            "prescribed_displacement.txt": PRESCRIBED_DISP,
            "point_load.txt": POINT_LOAD,
            "simulation_settings.txt": SIM_SETTINGS,
        }
        for fname, content in files.items():
            fpath = job_dir / fname
            if dry_run:
                print(f"  [dry] {fpath}")
            else:
                fpath.write_text(content, encoding="utf-8")
        print(f"{'[dry] ' if dry_run else ''}Created: {dir_name}")
        created += 1

    print(f"\n{'Would create' if dry_run else 'Created'} {created} job directories.")


if __name__ == "__main__":
    dry_run = "--dry-run" in sys.argv
    generate(dry_run=dry_run)

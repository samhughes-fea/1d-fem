"""
Generate mesh-refinement job variants for spectral benchmark families.
Creates job_benchmark_<family>_n<N> directories for N in N_VALUES,
mirroring the static linear job_0000_n<N> convention.
"""
from __future__ import annotations

import os
from pathlib import Path

JOBS_DIR = Path(__file__).resolve().parent
N_VALUES = [4, 8, 16, 32, 64, 128, 500]
L = 2.0  # beam length [m]

ELEMENT_TYPE = "LinearEulerBernoulliBeamElement3D"
ELEMENT_ROW = (
    "  {eid:<5} {n1:<5} {n2:<5} "
    "{etype:<40} 3                 3                 3                 "
    "0                 0                 3                 2"
)

MATERIAL_ROW = "{eid:<5} 2.1e+11     8.1000e+10   0.3    7850"
SECTION_ROW_EIGEN    = "{eid}  0.00131  0.0  3.23400e-07  2.08769e-06  2.60673e-08"
SECTION_ROW_BUCKLING = "{eid}  0.00131  0.0  3.23400e-07  2.08769e-06  2.60673e-08"
SECTION_ROW_HARMONIC = "{eid}  0.00131  0.0  3.23400e-07  2.08769e-06  2.60673e-08"

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

BENCHMARKS: dict[str, dict] = {
    "job_benchmark_eigen_cantilever": {
        "section_row": SECTION_ROW_EIGEN,
        "point_load": (
            "# Tip transverse load (eigen uses K and M only; load kept for completeness).\n"
            "[Point load]\n"
            "         [x]          [y]          [z]        [F_x]        [F_y]        [F_z]        [M_x]        [M_y]        [M_z]\n"
            "    2.000000     0.000000     0.000000     0.000000  -500.000000     0.000000     0.000000     0.000000     0.000000\n"
        ),
        "sim_settings": (
            "[Simulation]\n[Type]\neigen\n\n"
            "[Eigen]\nenabled = true\nnum_modes = 4\n"
        ),
    },
    "job_benchmark_linear_buckling_column": {
        "section_row": SECTION_ROW_BUCKLING,
        "point_load": (
            "# Axial compression at tip for linear-static prestress (Euler-type buckling).\n"
            "[Point load]\n"
            "         [x]          [y]          [z]        [F_x]        [F_y]        [F_z]        [M_x]        [M_y]        [M_z]\n"
            "    2.000000     0.000000     0.000000 -1000.000000     0.000000     0.000000     0.000000     0.000000     0.000000\n"
        ),
        "sim_settings": (
            "[Simulation]\n[Type]\nbuckling\n\n"
            "[Buckling]\nenabled = true\nnum_modes = 4\n"
            "buckling_prestress = linear_static\nbuckling_load_factor = 1.0\n"
        ),
    },
    "job_benchmark_harmonic_sdof": {
        "section_row": SECTION_ROW_HARMONIC,
        "point_load": (
            "# Tip transverse harmonic excitation.\n"
            "[Point load]\n"
            "         [x]          [y]          [z]        [F_x]        [F_y]        [F_z]        [M_x]        [M_y]        [M_z]\n"
            "    2.000000     0.000000     0.000000     0.000000  -500.000000     0.000000     0.000000     0.000000     0.000000\n"
        ),
        "sim_settings": (
            "[Simulation]\n[Type]\nharmonic\n\n"
            "[Harmonic]\nenabled = true\n"
            "frequency_min_hz = 8.0\nfrequency_max_hz = 12.0\n"
            "num_frequency_points = 5\nmodal_damping_ratio = 0.02\n"
        ),
    },
}


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


def _section_txt(n: int, template: str) -> str:
    lines = ["[Section]", "[element_id]  [A]  [I_x]  [I_y]  [I_z]  [J_t]"]
    for i in range(n):
        lines.append(template.format(eid=i))
    return "\n".join(lines) + "\n"


def _precurvature_txt(n: int) -> str:
    lines = [
        "# Reference curvature / twist (1/m) in local element frame.",
        "[Precurvature]",
        "[element_id] [k_x0] [k_y0] [k_z0]",
    ]
    for i in range(n):
        lines.append(f"{i} 0.0 0.0 0.0")
    return "\n".join(lines) + "\n"


def generate(dry_run: bool = False) -> None:
    created = 0
    for base_name, cfg in BENCHMARKS.items():
        for n in N_VALUES:
            dir_name = f"{base_name}_n{n}"
            job_dir = JOBS_DIR / dir_name
            if not dry_run:
                job_dir.mkdir(exist_ok=True)

            files = {
                "grid.txt": _grid_txt(n),
                "element.txt": _element_txt(n),
                "material.txt": _material_txt(n),
                "section.txt": _section_txt(n, cfg["section_row"]),
                "precurvature.txt": _precurvature_txt(n),
                "prescribed_displacement.txt": PRESCRIBED_DISP,
                "point_load.txt": cfg["point_load"],
                "simulation_settings.txt": cfg["sim_settings"],
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
    import sys
    dry_run = "--dry-run" in sys.argv
    generate(dry_run=dry_run)

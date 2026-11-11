"""
Generate missing mesh-refinement job variants for nonlinear static benchmark families.

For each load case, reads the base job directory and produces n-variant directories
by regenerating grid.txt and element.txt for each target n while copying all other
input files verbatim from the base job.

Usage:
    python jobs/_generate_nl_static_mesh_variants.py           # create missing dirs
    python jobs/_generate_nl_static_mesh_variants.py --dry-run # print what would be created
"""
from __future__ import annotations

import sys
from pathlib import Path

JOBS_DIR = Path(__file__).resolve().parent
L = 2.0  # beam length [m]

ELEMENT_TYPE = "NonlinearEulerBernoulliBeamElement3D"
_ELEM_ROW = (
    "{eid:<5} {n1:<5} {n2:<5} "
    "{etype:<42} 3                 3                 3                 "
    "0                 0                 3                 2"
)

# Files copied verbatim from base job (content is independent of n)
_STATIC_FILES = [
    "simulation_settings.txt",
    "prescribed_displacement.txt",
    "section.txt",
    "material.txt",
]

# Per load-case config: load_file is copied verbatim; missing_n lists dirs to create
BENCHMARKS: dict[str, dict] = {
    "job_benchmark_nl_static_cantilever_tip": {
        "load_file": "point_load.txt",
        "missing_n": [4, 8, 32, 128],
    },
    "job_benchmark_nl_static_midspan_point": {
        "load_file": "point_load.txt",
        "missing_n": [4, 8, 32, 128],
    },
    "job_benchmark_nl_static_quarter_point": {
        "load_file": "point_load.txt",
        "missing_n": [4, 8, 32, 128, 500],
    },
    "job_benchmark_nl_static_udl": {
        "load_file": "distributed_load.txt",
        "missing_n": [8, 32, 128],
    },
    "job_benchmark_nl_static_triangular": {
        "load_file": "distributed_load.txt",
        "missing_n": [8, 32, 128],
    },
    "job_benchmark_nl_static_parabolic": {
        "load_file": "distributed_load.txt",
        "missing_n": [8, 32, 128],
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
        rows.append(_ELEM_ROW.format(eid=i, n1=i, n2=i + 1, etype=ELEMENT_TYPE))
    return "\n".join(rows) + "\n"


def _scale_section(base_section: str, n: int) -> str:
    """Re-generate section.txt with n rows by repeating the per-element property line."""
    lines = base_section.splitlines()
    header_lines = [l for l in lines if l.startswith("[") or l.startswith("#")]
    # Extract the property columns from first data row (skip header lines)
    data_lines = [l for l in lines if l and not l.startswith("[") and not l.startswith("#")]
    if not data_lines:
        return base_section
    # Take everything after the element_id column
    first_cols = data_lines[0].split(None, 1)
    if len(first_cols) < 2:
        return base_section
    props = first_cols[1]
    out = ["[Section]", "[element_id]  [A]  [I_x]  [I_y]  [I_z]  [J_t]"]
    for i in range(n):
        out.append(f"{i}  {props}")
    return "\n".join(out) + "\n"


def _scale_material(base_material: str, n: int) -> str:
    """Re-generate material.txt with n rows."""
    lines = base_material.splitlines()
    data_lines = [l for l in lines if l and not l.startswith("[") and not l.startswith("#")]
    if not data_lines:
        return base_material
    first_cols = data_lines[0].split(None, 1)
    if len(first_cols) < 2:
        return base_material
    props = first_cols[1]
    out = ["[Material]", "[element_id]  [E]         [G]          [nu]   [rho]"]
    for i in range(n):
        out.append(f"{i:<5} {props}")
    return "\n".join(out) + "\n"


def generate(dry_run: bool = False) -> None:
    created = 0
    skipped = 0
    for base_name, cfg in BENCHMARKS.items():
        base_dir = JOBS_DIR / base_name
        if not base_dir.is_dir():
            print(f"WARNING: base job not found: {base_dir}")
            continue

        # Read static files from base job
        static_content: dict[str, str] = {}
        for fname in _STATIC_FILES:
            fpath = base_dir / fname
            if not fpath.is_file():
                print(f"WARNING: {base_name}/{fname} not found; skipping base.")
                break
            static_content[fname] = fpath.read_text(encoding="utf-8")
        else:
            pass

        load_file = cfg["load_file"]
        load_path = base_dir / load_file
        if not load_path.is_file():
            print(f"WARNING: {base_name}/{load_file} not found; skipping.")
            continue
        load_content = load_path.read_text(encoding="utf-8")

        for n in cfg["missing_n"]:
            dir_name = f"{base_name}_n{n}"
            job_dir = JOBS_DIR / dir_name
            if job_dir.is_dir():
                skipped += 1
                continue

            files: dict[str, str] = {
                "grid.txt": _grid_txt(n),
                "element.txt": _element_txt(n),
                "section.txt": _scale_section(static_content.get("section.txt", ""), n),
                "material.txt": _scale_material(static_content.get("material.txt", ""), n),
                "simulation_settings.txt": static_content.get("simulation_settings.txt", ""),
                "prescribed_displacement.txt": static_content.get("prescribed_displacement.txt", ""),
                load_file: load_content,
            }

            if dry_run:
                print(f"[dry] Would create: {dir_name}")
                for fname in files:
                    print(f"  {fname}")
            else:
                job_dir.mkdir(exist_ok=True)
                for fname, content in files.items():
                    (job_dir / fname).write_text(content, encoding="utf-8")
                print(f"Created: {dir_name}")
            created += 1

    label = "Would create" if dry_run else "Created"
    print(f"\n{label} {created} directories ({skipped} already existed, skipped).")


if __name__ == "__main__":
    dry_run = "--dry-run" in sys.argv
    generate(dry_run=dry_run)

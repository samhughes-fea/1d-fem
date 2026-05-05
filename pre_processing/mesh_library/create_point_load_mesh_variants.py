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

from pre_processing.mesh_library.benchmark_generation import (
    linear_static_family_variant,
    nonlinear_static_point_benchmark,
)
from pre_processing.mesh_library.job_generation import write_job_from_spec

JOBS_DIR = REPO_ROOT / "jobs"
L = 2.0  # m (must match mesh_generator default)

VARIANT_NS = [4, 8, 16, 32, 64, 128, 500]  # n500 = Abaqus validation reference (converged)
# (base_id, load_label, load_formula, x_position)
BASE_JOBS = [
    (0, "End load", "P(x=L)", L),
    (1, "Midspan load", "P(x=L/2)", L / 2),
    (2, "Quarter-point load", "P(x=L/4)", L / 4),
]


def main() -> None:
    for base_id, load_label, _formula, x_load in BASE_JOBS:
        for n in VARIANT_NS:
            name = f"job_{base_id:04d}_n{n}"
            spec = linear_static_family_variant(
                name,
                num_elements=n,
                x_load=x_load,
                load_label=load_label,
                load_formula=_formula,
                element_type="EulerBernoulliBeamElement3D",
            )
            write_job_from_spec(JOBS_DIR, spec)
            print(f"Created {name}")

    benchmark_jobs = [
        ("job_benchmark_nl_static_cantilever_tip", 4, "End load", "P(x=L)", L),
        ("job_benchmark_nl_static_midspan_point", 4, "Midspan load", "P(x=L/2)", L / 2),
        ("job_benchmark_nl_static_quarter_point", 4, "Quarter-point load", "P(x=L/4)", L / 4),
        ("job_benchmark_nl_static_cantilever_tip_n16", 16, "End load", "P(x=L)", L),
        ("job_benchmark_nl_static_midspan_point_n16", 16, "Midspan load", "P(x=L/2)", L / 2),
        ("job_benchmark_nl_static_quarter_point_n16", 16, "Quarter-point load", "P(x=L/4)", L / 4),
        ("job_benchmark_nl_static_cantilever_tip_n64", 64, "End load", "P(x=L)", L),
        ("job_benchmark_nl_static_midspan_point_n64", 64, "Midspan load", "P(x=L/2)", L / 2),
        ("job_benchmark_nl_static_quarter_point_n64", 64, "Quarter-point load", "P(x=L/4)", L / 4),
        ("job_benchmark_nl_static_cantilever_tip_n500", 500, "End load", "P(x=L)", L),
        ("job_benchmark_nl_static_midspan_point_n500", 500, "Midspan load", "P(x=L/2)", L / 2),
        ("job_benchmark_nl_static_quarter_point_n500", 500, "Quarter-point load", "P(x=L/4)", L / 4),
    ]
    for name, n, label, formula, x_load in benchmark_jobs:
        write_job_from_spec(
            JOBS_DIR,
            nonlinear_static_point_benchmark(
                name,
                num_elements=n,
                load_label=label,
                load_formula=formula,
                x_load=x_load,
            ),
        )
        print(f"Created {name}")


if __name__ == "__main__":
    main()

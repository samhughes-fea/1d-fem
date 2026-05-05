#!/usr/bin/env python3
"""
Create mesh variants for distributed-load jobs (job_0003, job_0004, job_0005).
Generates job_0003_n4, n8, n16, n32, n64, n128 (UDL), job_0004_n*, job_0005_n* (triangular, parabolic)
using the mesh library (pre_processing/mesh_library/schemes/mesh_generator.py) for geometry and
properties, then adds job-specific load/BC/simulation files.
Run from repo root: python pre_processing/mesh_library/create_distributed_mesh_variants.py
"""
import sys

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from pre_processing.mesh_library.benchmark_generation import (
    linear_static_family_variant,
    nonlinear_static_distributed_benchmark,
)
from pre_processing.mesh_library.job_generation import write_job_from_spec

JOBS_DIR = REPO_ROOT / "jobs"
L = 2.0  # m (must match mesh_generator default)
w = 500.0  # N/m


def _udl(x: float) -> float:
    return -w


def _triangular(x: float) -> float:
    return -w * (x / L) if L else 0


def _parabolic(x: float) -> float:
    return -w * (x / L) ** 2 if L else 0

VARIANT_NS = [4, 8, 16, 32, 64, 128, 500]  # n500 = Abaqus validation reference (converged)
BASE_JOBS = [
    (3, "UDL", _udl),
    (4, "TRIANGULAR", _triangular),
    (5, "PARABOLIC", _parabolic),
]
LOAD_FORMULAS = {"UDL": "q(x) = w", "TRIANGULAR": "q(x) = w * (x/L)", "PARABOLIC": "q(x) = w * (x/L)^2"}

JOB_NAME_OVERRIDES = {
    ("UDL", 4): "job_benchmark_nl_static_udl_n4",
    ("TRIANGULAR", 4): "job_benchmark_nl_static_triangular_n4",
    ("PARABOLIC", 4): "job_benchmark_nl_static_parabolic_n4",
    ("UDL", 16): "job_benchmark_nl_static_udl_n16",
    ("TRIANGULAR", 16): "job_benchmark_nl_static_triangular_n16",
    ("PARABOLIC", 16): "job_benchmark_nl_static_parabolic_n16",
    ("UDL", 64): "job_benchmark_nl_static_udl_n64",
    ("TRIANGULAR", 64): "job_benchmark_nl_static_triangular_n64",
    ("PARABOLIC", 64): "job_benchmark_nl_static_parabolic_n64",
    ("TRIANGULAR", 500): "job_benchmark_nl_static_triangular_n500",
    ("PARABOLIC", 500): "job_benchmark_nl_static_parabolic_n500",
}

def resolve_job_name(base_id: int, load_type: str, n: int) -> str:
    return JOB_NAME_OVERRIDES.get((load_type, n), f"job_{base_id:04d}_n{n}")


def main() -> None:
    for base_id, load_type, fy_fn in BASE_JOBS:
        for n in VARIANT_NS:
            name = resolve_job_name(base_id, load_type, n)
            spec = linear_static_family_variant(
                name,
                num_elements=n,
                distributed=(load_type, LOAD_FORMULAS[load_type], fy_fn),
                element_type="LinearEulerBernoulliBeamElement3D",
            )
            write_job_from_spec(JOBS_DIR, spec)
            print(f"Created {name}")

    benchmark_jobs = [
        ("job_benchmark_nl_static_udl", 4, "UDL", LOAD_FORMULAS["UDL"], _udl),
        ("job_benchmark_nl_static_triangular", 4, "TRIANGULAR", LOAD_FORMULAS["TRIANGULAR"], _triangular),
        ("job_benchmark_nl_static_parabolic", 4, "PARABOLIC", LOAD_FORMULAS["PARABOLIC"], _parabolic),
        ("job_benchmark_nl_static_udl_n4", 4, "UDL", LOAD_FORMULAS["UDL"], _udl),
        ("job_benchmark_nl_static_triangular_n4", 4, "TRIANGULAR", LOAD_FORMULAS["TRIANGULAR"], _triangular),
        ("job_benchmark_nl_static_parabolic_n4", 4, "PARABOLIC", LOAD_FORMULAS["PARABOLIC"], _parabolic),
        ("job_benchmark_nl_static_udl_n16", 16, "UDL", LOAD_FORMULAS["UDL"], _udl),
        ("job_benchmark_nl_static_triangular_n16", 16, "TRIANGULAR", LOAD_FORMULAS["TRIANGULAR"], _triangular),
        ("job_benchmark_nl_static_parabolic_n16", 16, "PARABOLIC", LOAD_FORMULAS["PARABOLIC"], _parabolic),
        ("job_benchmark_nl_static_udl_n64", 64, "UDL", LOAD_FORMULAS["UDL"], _udl),
        ("job_benchmark_nl_static_triangular_n64", 64, "TRIANGULAR", LOAD_FORMULAS["TRIANGULAR"], _triangular),
        ("job_benchmark_nl_static_parabolic_n64", 64, "PARABOLIC", LOAD_FORMULAS["PARABOLIC"], _parabolic),
        ("job_benchmark_nl_static_udl_n500", 500, "UDL", LOAD_FORMULAS["UDL"], _udl),
        ("job_benchmark_nl_static_triangular_n500", 500, "TRIANGULAR", LOAD_FORMULAS["TRIANGULAR"], _triangular),
        ("job_benchmark_nl_static_parabolic_n500", 500, "PARABOLIC", LOAD_FORMULAS["PARABOLIC"], _parabolic),
    ]
    for name, n, load_type, formula, fy_fn in benchmark_jobs:
        write_job_from_spec(
            JOBS_DIR,
            nonlinear_static_distributed_benchmark(
                name,
                num_elements=n,
                load_type=load_type,
                load_formula=formula,
                fy_fn=fy_fn,
            ),
        )
        print(f"Created {name}")


if __name__ == "__main__":
    main()

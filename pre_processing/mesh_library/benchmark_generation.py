from __future__ import annotations

from dataclasses import replace

from pre_processing.mesh_library.job_generation import (
    NONLINEAR_STATIC_SETTINGS,
    DistributedLoadSpec,
    ElementSpec,
    JobSpec,
    MeshSpec,
    PointLoadSpec,
    distributed_samples,
    fixed_cantilever_support,
)


LENGTH = 2.0
POINT_LOAD = -500.0
LINEAR_ELEMENT = "LinearEulerBernoulliBeamElement3D"
NONLINEAR_ELEMENT = "NonlinearEulerBernoulliBeamElement3D"


def nonlinear_static_point_benchmark(job_name: str, *, num_elements: int, load_label: str, load_formula: str, x_load: float) -> JobSpec:
    return JobSpec(
        job_name=job_name,
        mesh=MeshSpec(length=LENGTH, num_elements=num_elements),
        element=ElementSpec(element_type=NONLINEAR_ELEMENT),
        prescribed_displacements=fixed_cantilever_support(),
        point_loads=(PointLoadSpec(x=x_load, F_y=POINT_LOAD),),
        simulation_settings_text=NONLINEAR_STATIC_SETTINGS,
        metadata_comments={"Load Type": load_label, "Load Formula": load_formula},
        readme_reference_text=_point_reference_readme(job_name) if num_elements == 500 else None,
    )


def nonlinear_static_distributed_benchmark(job_name: str, *, num_elements: int, load_type: str, load_formula: str, fy_fn) -> JobSpec:
    return JobSpec(
        job_name=job_name,
        mesh=MeshSpec(length=LENGTH, num_elements=num_elements),
        element=ElementSpec(element_type=NONLINEAR_ELEMENT),
        prescribed_displacements=fixed_cantilever_support(),
        distributed_load=DistributedLoadSpec(samples=distributed_samples(num_elements, LENGTH, fy_fn)),
        simulation_settings_text=NONLINEAR_STATIC_SETTINGS,
        metadata_comments={"Load Type": load_type, "Load Formula": load_formula},
        readme_reference_text=_distributed_reference_readme(load_type) if num_elements == 500 else None,
    )


def linear_static_family_variant(job_name: str, *, num_elements: int, x_load: float | None = None, load_label: str | None = None, load_formula: str | None = None, distributed: tuple[str, str, object] | None = None, element_type: str = LINEAR_ELEMENT) -> JobSpec:
    base = JobSpec(
        job_name=job_name,
        mesh=MeshSpec(length=LENGTH, num_elements=num_elements),
        element=ElementSpec(element_type=element_type),
        prescribed_displacements=fixed_cantilever_support(),
    )
    if distributed is not None:
        load_type, formula, fy_fn = distributed
        return replace(
            base,
            distributed_load=DistributedLoadSpec(samples=distributed_samples(num_elements, LENGTH, fy_fn)),
            metadata_comments={"Load Type": load_type, "Load Formula": formula},
        )
    return replace(
        base,
        point_loads=(PointLoadSpec(x=float(x_load), F_y=POINT_LOAD),),
        metadata_comments={"Load Type": str(load_label), "Load Formula": str(load_formula)},
    )


def _point_reference_readme(job_name: str) -> str:
    n64_job = job_name.replace("_n500", "_n64")
    slug = job_name.removeprefix("job_benchmark_nl_static_").removesuffix("_n500").replace("_", " ")
    return (
        f"# Fine Abaqus reference job: nonlinear-static cantilever {slug}\n\n"
        "This job is the first planned fine-reference Abaqus anchor for the nonlinear-static validation suite.\n\n"
        "## Purpose\n\n"
        "Provide the high-fidelity external reference counterpart for:\n\n"
        f"- [`{n64_job}/`](../{n64_job}/)\n\n"
        f"using the same canonical {slug} nonlinear-static cantilever definition at a much finer mesh level.\n\n"
        "## Phase-3 scope\n\n"
        "This phase defines the job root and contract only. It does not yet claim a completed fine-reference comparison.\n"
    )


def _distributed_reference_readme(load_type: str) -> str:
    root_slug = load_type.lower()
    title = load_type.lower().capitalize()
    n64_job = f"job_benchmark_nl_static_{root_slug}_n64"
    return (
        f"# Fine Abaqus reference job: nonlinear-static cantilever {title} distributed load\n\n"
        "This job is the planned fine-reference Abaqus anchor for the nonlinear-static distributed-load validation suite.\n\n"
        "## Purpose\n\n"
        "Provide the high-fidelity external reference counterpart for:\n\n"
        f"- [`{n64_job}/`](../{n64_job}/)\n\n"
        f"using the same canonical {title.lower()} distributed nonlinear-static cantilever definition at a much finer mesh level.\n\n"
        "## Phase-3 scope\n\n"
        "This phase defines the job root and contract only. It does not yet claim a completed fine-reference comparison.\n"
    )

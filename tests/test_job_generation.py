from __future__ import annotations

from pathlib import Path

from pre_processing.mesh_library.benchmark_generation import (
    nonlinear_static_distributed_benchmark,
    nonlinear_static_point_benchmark,
)
from pre_processing.mesh_library.job_generation import (
    ElementSpec,
    JobSpec,
    MeshSpec,
    PointLoadSpec,
    fixed_cantilever_support,
    write_job_from_spec,
)


def test_generic_job_writer_preserves_parser_contract(tmp_path: Path) -> None:
    spec = JobSpec(
        job_name="job_contract_case",
        mesh=MeshSpec(length=2.0, num_elements=4),
        element=ElementSpec(element_type="NonlinearEulerBernoulliBeamElement3D"),
        prescribed_displacements=fixed_cantilever_support(),
        point_loads=(PointLoadSpec(x=1.0, F_y=-500.0),),
        simulation_settings_text="[Simulation]\n[Type]\nstatic_nonlinear\n",
        metadata_comments={"Load Type": "Midspan load", "Load Formula": "P(x=L/2)"},
    )

    job_dir = write_job_from_spec(tmp_path, spec)

    assert (job_dir / "grid.txt").read_text(encoding="utf-8").startswith("[Grid]\n[node_id]")
    element_text = (job_dir / "element.txt").read_text(encoding="utf-8")
    assert "[Element]" in element_text
    assert "NonlinearEulerBernoulliBeamElement3D" in element_text
    point_load_text = (job_dir / "point_load.txt").read_text(encoding="utf-8")
    assert "# Load Type: Midspan load" in point_load_text
    assert "[Point load]" in point_load_text
    assert (job_dir / "prescribed_displacement.txt").read_text(encoding="utf-8").count("# Fixed support") == 6


def test_point_benchmark_profile_writes_reference_readme_for_n500(tmp_path: Path) -> None:
    spec = nonlinear_static_point_benchmark(
        "job_benchmark_nl_static_midspan_point_n500",
        num_elements=500,
        load_label="Midspan load",
        load_formula="P(x=L/2)",
        x_load=1.0,
    )

    job_dir = write_job_from_spec(tmp_path, spec)

    readme = (job_dir / "README_REFERENCE.md").read_text(encoding="utf-8")
    assert "job_benchmark_nl_static_midspan_point_n64" in readme
    assert "fine-reference Abaqus anchor" in readme


def test_distributed_benchmark_profile_writes_load_samples_and_readme(tmp_path: Path) -> None:
    spec = nonlinear_static_distributed_benchmark(
        "job_benchmark_nl_static_triangular_n500",
        num_elements=4,
        load_type="TRIANGULAR",
        load_formula="q(x) = w * (x/L)",
        fy_fn=lambda x: -500.0 * (x / 2.0),
    )

    job_dir = write_job_from_spec(tmp_path, spec)
    text = (job_dir / "distributed_load.txt").read_text(encoding="utf-8")
    assert "# Load Type: TRIANGULAR" in text
    assert "-500.000000" in text
    assert "[Distributed load]" in text

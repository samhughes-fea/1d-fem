from __future__ import annotations

from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parent.parent


def test_validation_directory_standard_files_exist() -> None:
    assert (PROJECT_ROOT / "docs" / "conventions" / "ABAQUS_VALIDATION_DIRECTORY_STANDARD.md").is_file()
    for family in ("static", "eigen", "buckling", "transient", "harmonic"):
        family_dir = PROJECT_ROOT / "post_processing" / "validation_visualisers" / family
        assert family_dir.is_dir()
        assert any(p.name.startswith("README_VALIDATION_") for p in family_dir.iterdir())
        for sub in ("reference_cases", "plots", "tables", "output"):
            assert (family_dir / sub).is_dir()


def test_nonlinear_static_validation_phase1_assets_exist() -> None:
    root = PROJECT_ROOT
    assert (root / "post_processing" / "validation_visualisers" / "static" / "README_VALIDATION_NONLINEAR_STATIC.md").is_file()
    assert (root / "jobs" / "job_benchmark_nl_static_cantilever_tip" / "simulation_settings.txt").is_file()
    assert (root / "docs" / "conventions" / "NONLINEAR_STATIC_ABAQUS_VALIDATION_PLAN.md").is_file()


def test_nonlinear_static_tip_displacement_contract_assets_exist() -> None:
    root = PROJECT_ROOT
    assert (root / "docs" / "conventions" / "NONLINEAR_STATIC_TIP_DISPLACEMENT_CONTRACT.md").is_file()
    text = (root / "post_processing" / "validation_visualisers" / "static" / "README_VALIDATION_NONLINEAR_STATIC.md").read_text(encoding="utf-8")
    assert "tip displacement vs load level" in text


def test_nonlinear_static_tip_history_writer_outputs_csv(tmp_path) -> None:
    from processing.static.results.save_nonlinear_static_validation_summary import write_nonlinear_static_tip_history_summary

    out = write_nonlinear_static_tip_history_summary(
        primary_results_dir=str(tmp_path),
        job_name="nl_static_case",
        load_factors=[0.5, 1.0],
        tip_displacements=[-1.0e-3, -2.0e-3],
    )
    p = Path(out)
    assert p.is_file()
    text = p.read_text(encoding="utf-8")
    assert "load_factor,tip_displacement" in text


def test_nonlinear_static_abaqus_reference_contract_assets_exist() -> None:
    root = PROJECT_ROOT
    doc = root / "docs" / "conventions" / "NONLINEAR_STATIC_ABAQUS_REFERENCE_CONTRACT.md"
    assert doc.is_file()
    text = doc.read_text(encoding="utf-8")
    assert "tip_load_history.csv" in text
    readme = (root / "post_processing" / "validation_visualisers" / "static" / "README_VALIDATION_NONLINEAR_STATIC.md").read_text(encoding="utf-8")
    assert "tip_load_history.csv" in readme


def test_nonlinear_static_point_family_roots_exist() -> None:
    root = PROJECT_ROOT / "jobs"
    for name in (
        "job_benchmark_nl_static_cantilever_tip",
        "job_benchmark_nl_static_midspan_point",
        "job_benchmark_nl_static_quarter_point",
    ):
        d = root / name
        assert d.is_dir()
        assert (d / "simulation_settings.txt").is_file()


def test_nonlinear_static_distributed_family_roots_exist() -> None:
    root = PROJECT_ROOT / "jobs"
    for name in (
        "job_benchmark_nl_static_udl",
        "job_benchmark_nl_static_triangular",
        "job_benchmark_nl_static_parabolic",
    ):
        d = root / name
        assert d.is_dir()
        assert (d / "simulation_settings.txt").is_file()


def test_nonlinear_static_n16_ladder_seed_exists() -> None:
    root = PROJECT_ROOT / "jobs"
    for name in (
        "job_benchmark_nl_static_cantilever_tip_n16",
        "job_benchmark_nl_static_midspan_point_n16",
        "job_benchmark_nl_static_quarter_point_n16",
        "job_benchmark_nl_static_udl_n16",
        "job_benchmark_nl_static_triangular_n16",
        "job_benchmark_nl_static_parabolic_n16",
    ):
        d = root / name
        assert d.is_dir()
        assert (d / "simulation_settings.txt").is_file()


def test_nonlinear_static_n4_ladder_seed_exists() -> None:
    root = PROJECT_ROOT / "jobs"
    for name in (
        "job_benchmark_nl_static_cantilever_tip",
        "job_benchmark_nl_static_midspan_point",
        "job_benchmark_nl_static_quarter_point",
        "job_benchmark_nl_static_udl_n4",
        "job_benchmark_nl_static_triangular_n4",
        "job_benchmark_nl_static_parabolic_n4",
    ):
        d = root / name
        assert d.is_dir()
        assert (d / "simulation_settings.txt").is_file()


def test_nonlinear_static_n64_point_seed_exists() -> None:
    root = PROJECT_ROOT / "jobs"
    for name in (
        "job_benchmark_nl_static_cantilever_tip_n64",
        "job_benchmark_nl_static_midspan_point_n64",
        "job_benchmark_nl_static_quarter_point_n64",
    ):
        d = root / name
        assert d.is_dir()
        assert (d / "simulation_settings.txt").is_file()


def test_nonlinear_static_n64_distributed_seed_exists() -> None:
    root = PROJECT_ROOT / "jobs"
    for name in (
        "job_benchmark_nl_static_udl_n64",
        "job_benchmark_nl_static_triangular_n64",
        "job_benchmark_nl_static_parabolic_n64",
    ):
        d = root / name
        assert d.is_dir()
        assert (d / "simulation_settings.txt").is_file()


def test_nonlinear_static_tip_n500_reference_job_exists() -> None:
    root = PROJECT_ROOT / "jobs" / "job_benchmark_nl_static_cantilever_tip_n500"
    assert root.is_dir()
    assert (root / "simulation_settings.txt").is_file()
    assert (root / "README_REFERENCE.md").is_file()


def test_nonlinear_static_suite_docs_record_common_fem_artifact_contract() -> None:
    root = PROJECT_ROOT
    suite_doc = (root / "docs" / "conventions" / "NONLINEAR_STATIC_ABAQUS_VALIDATION_SUITE.md").read_text(encoding="utf-8")
    assert "nonlinear_static_validation" in suite_doc
    assert "tip_load_history.csv" in suite_doc
    readme = (root / "post_processing" / "validation_visualisers" / "static" / "README_VALIDATION_NONLINEAR_STATIC.md").read_text(encoding="utf-8")
    assert "tip_load_history.csv" in readme


def test_nonlinear_static_suite_docs_record_common_abaqus_artifact_contract() -> None:
    root = PROJECT_ROOT
    suite_doc = (root / "docs" / "conventions" / "NONLINEAR_STATIC_ABAQUS_VALIDATION_SUITE.md").read_text(encoding="utf-8")
    assert "U_global.csv" in suite_doc
    assert "rotation_source.txt" in suite_doc
    assert "frame_index" in suite_doc
    readme = (root / "post_processing" / "validation_visualisers" / "static" / "README_VALIDATION_NONLINEAR_STATIC.md").read_text(encoding="utf-8")
    assert "U_global.csv" in readme
    assert "rotation_source.txt" in readme

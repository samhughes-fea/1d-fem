from __future__ import annotations

from pathlib import Path

from post_processing.validation_visualisers.static.compare_nonlinear_static_tip_history import (
    benchmark_paths,
    compare_tip_histories,
    infer_abaqus_reference_job_name,
    run_nonlinear_static_tip_history_validation,
    run_nonlinear_static_tip_history_validation_from_benchmark_paths,
    write_tip_history_comparison_csv,
)
from post_processing.validation_visualisers.static.run_nonlinear_static_validation import (
    discover_suite_job_names,
    evaluate_benchmark_pass_fail,
    run_benchmark,
    run_suite,
)


def test_compare_tip_histories_and_write_report(tmp_path):
    fem = tmp_path / "fem.csv"
    aba = tmp_path / "abaqus.csv"
    fem.write_text("load_factor,tip_displacement\n0.5,-0.001\n1.0,-0.002\n", encoding="utf-8")
    aba.write_text("load_factor,tip_displacement\n0.5,-0.0011\n1.0,-0.0019\n", encoding="utf-8")

    rows = compare_tip_histories(fem, aba)
    assert len(rows) == 2
    assert rows[0]["abs_error"] > 0.0
    assert rows[0]["load_factor_alignment_error"] == 0.0

    out = tmp_path / "output" / "comparison.csv"
    write_tip_history_comparison_csv(rows, out)
    assert out.is_file()
    text = out.read_text(encoding="utf-8")
    assert "abs_error" in text
    assert "load_factor_alignment_error" in text


def test_compare_tip_histories_interpolates_abaqus_to_fem_load_factor(tmp_path):
    fem = tmp_path / "fem_interp.csv"
    aba = tmp_path / "abaqus_interp.csv"
    fem.write_text("load_factor,tip_displacement\n0.5,-0.001\n1.0,-0.002\n", encoding="utf-8")
    aba.write_text("load_factor,tip_displacement\n0.0,0.0\n1.0,-0.002\n", encoding="utf-8")

    rows = compare_tip_histories(fem, aba)

    assert len(rows) == 2
    assert rows[0]["load_factor_abaqus"] == 0.5
    assert rows[0]["tip_displacement_abaqus"] == -0.001
    assert rows[0]["load_factor_alignment_error"] == 0.0
    assert rows[0]["abs_error"] == 0.0


def test_read_tip_history_csv_is_header_aware_for_abaqus_layout(tmp_path):
    from post_processing.validation_visualisers.static.compare_nonlinear_static_tip_history import read_tip_history_csv

    aba = tmp_path / "abaqus_header_aware.csv"
    aba.write_text(
        "frame_index,load_factor,tip_displacement\n0,0.0,0.0\n1,1.0,-0.003\n",
        encoding="utf-8",
    )

    rows = read_tip_history_csv(aba)

    assert rows == [(0.0, 0.0), (1.0, -0.003)]


def test_compare_tip_histories_handles_mixed_fem_and_abaqus_csv_layouts(tmp_path):
    fem = tmp_path / "fem_mixed.csv"
    aba = tmp_path / "abaqus_mixed.csv"
    fem.write_text("load_factor,tip_displacement\n0.5,-0.001\n1.0,-0.002\n", encoding="utf-8")
    aba.write_text(
        "frame_index,load_factor,tip_displacement\n0,0.0,0.0\n1,1.0,-0.002\n",
        encoding="utf-8",
    )

    rows = compare_tip_histories(fem, aba)

    assert rows[0]["load_factor_abaqus"] == 0.5
    assert rows[0]["tip_displacement_abaqus"] == -0.001
    assert rows[1]["tip_displacement_abaqus"] == -0.002


def test_run_nonlinear_static_tip_history_validation_writes_outputs(tmp_path):
    fem = tmp_path / "fem.csv"
    aba = tmp_path / "abaqus.csv"
    fem.write_text("load_factor,tip_displacement\n0.5,-0.001\n1.0,-0.002\n", encoding="utf-8")
    aba.write_text("load_factor,tip_displacement\n0.0,0.0\n1.0,-0.002\n", encoding="utf-8")
    out = run_nonlinear_static_tip_history_validation(
        fem_csv=fem,
        abaqus_csv=aba,
        output_dir=tmp_path / "validation_out",
    )
    assert Path(out["comparison_csv"]).is_file()
    assert Path(out["summary_csv"]).is_file()
    assert Path(out["plot_png"]).is_file()
    summary_text = Path(out["summary_csv"]).read_text(encoding="utf-8")
    assert "max_load_factor_alignment_error" in summary_text


def test_benchmark_paths_and_missing_reference_tolerant_flow(tmp_path):
    paths = benchmark_paths(tmp_path / "results", tmp_path / "abaqus", "job_benchmark_nl_static_cantilever_tip")
    assert paths["fem_csv"].endswith("job_benchmark_nl_static_cantilever_tip_tip_load_history.csv")
    out = run_nonlinear_static_tip_history_validation_from_benchmark_paths(
        results_root=tmp_path / "results",
        abaqus_root=tmp_path / "abaqus",
        job_name="job_benchmark_nl_static_cantilever_tip",
    )
    assert out["ready"] is False
    assert Path(out["status_file"]).is_file()


def test_infer_abaqus_reference_job_name_for_n64_pairing() -> None:
    assert infer_abaqus_reference_job_name("job_benchmark_nl_static_cantilever_tip_n64") == "job_benchmark_nl_static_cantilever_tip_n500"


def test_explicit_abaqus_reference_override_in_benchmark_runner(tmp_path):
    out = run_benchmark(
        tmp_path / "results",
        job_name="job_benchmark_nl_static_cantilever_tip_n64",
        abaqus_reference_job_name="job_benchmark_nl_static_cantilever_tip_n500",
    )
    assert out["ready"] is False
    assert out["abaqus_reference_job_name"] == "job_benchmark_nl_static_cantilever_tip_n500"


def test_run_benchmark_runner_returns_missing_reference_status(tmp_path):
    out = run_benchmark(tmp_path / "results")
    assert out["ready"] is False
    assert Path(out["status_file"]).is_file()


def test_discover_suite_job_names_finds_nonlinear_static_benchmarks():
    names = discover_suite_job_names()
    assert "job_benchmark_nl_static_cantilever_tip" in names


def test_run_suite_writes_summary(tmp_path):
    out = run_suite(tmp_path / "results")
    summary = Path(out["summary_csv"])
    assert summary.is_file()
    text = summary.read_text(encoding="utf-8")
    assert "job_name,ready,passed,status_file,max_load_factor_alignment_error,max_abs_error,mean_abs_error" in text


def test_evaluate_benchmark_pass_fail_accepts_first_reference_scale() -> None:
    assert evaluate_benchmark_pass_fail(
        {
            "max_load_factor_alignment_error": 0.0,
            "max_abs_error": 2.8889808153266788e-06,
        }
    ) is True


def test_evaluate_benchmark_pass_fail_rejects_large_error() -> None:
    assert evaluate_benchmark_pass_fail(
        {
            "max_load_factor_alignment_error": 0.0,
            "max_abs_error": 1.0e-3,
        }
    ) is False

"""Beam stages report on the shared logger."""

from __future__ import annotations

import pytest

from onedfem.stage_log import LOGGER_NAME, StageCollector, StageReport, collecting, fold_reports, stage


def test_a_stage_records_seconds_and_an_error_keeps_its_time() -> None:
    with collecting() as sink:
        with stage("onedfem", "kept"):
            pass
        with pytest.raises(RuntimeError):
            with stage("onedfem", "failed"):
                raise RuntimeError("stage failed")
    folded = fold_reports(sink.take())
    assert folded["kept"] >= 0.0
    failed = folded["failed"]
    assert isinstance(failed, dict)
    assert failed["seconds"] >= 0.0
    assert failed["error"] == "stage failed"


def test_public_contract() -> None:
    import dataclasses
    import logging

    assert LOGGER_NAME == "blade.stage"
    assert [item.name for item in dataclasses.fields(StageReport)] == [
        "package", "stage", "seconds", "error", "depth", "started"]
    assert callable(stage) and callable(fold_reports)
    assert issubclass(StageCollector, logging.Handler)


def test_a_processing_stage_keeps_its_class_log_and_reports(tmp_path) -> None:
    import logging

    from onedfem.processing.common.stage_logging import init_stage_logger, timed_operation

    class ExampleStage:
        def __init__(self) -> None:
            self._log = init_stage_logger("ExampleStage", tmp_path / "primary_results")

        @timed_operation
        def run(self) -> int:
            self._log.info("hello")
            return 3

    with collecting() as sink:
        assert ExampleStage().run() == 3
    found = sink.take()
    assert found[0].package == "onedfem"
    assert found[0].stage == "ExampleStage"
    assert found[0].error is None
    assert found[0].seconds >= 0.0
    text = (tmp_path / "logs" / "ExampleStage.log").read_text(encoding="utf-8")
    assert "hello" in text
    assert logging.getLogger("onedfem.processing.stage.ExampleStage").propagate is False

# processing/common/stage_logging.py
"""Per-stage file logging for processing operation classes (spectral, harmonic, transient, etc.)."""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Optional, Union


def init_stage_logger(
    class_name: str, job_results_dir: Optional[Union[str, Path]]
) -> logging.Logger:
    """
    Return a logger with optional ``logs/<ClassName>.log`` under the job tree.

    When ``job_results_dir`` is set, the log directory is ``parent(logs)`` of that path
    (same layout as static ``PrepareLocalSystem`` when given a path under ``primary_results``).

    For **element-level** file logging (per stiffness/force matrix under element subfolders),
    use ``pre_processing.element_library.base_logger_operator.BaseLoggerOperator`` instead;
    this helper is for **job-stage** processing operations (spectral, harmonic, transient, etc.).
    """
    log = logging.getLogger(f"processing.stage.{class_name}")
    log.handlers.clear()
    log.setLevel(logging.DEBUG)
    log.propagate = False
    if job_results_dir:
        p = Path(job_results_dir)
        log_dir = p.parent / "logs"
        log_dir.mkdir(parents=True, exist_ok=True)
        log_path = log_dir / f"{class_name}.log"
        try:
            fh = logging.FileHandler(log_path, mode="w", encoding="utf-8")
            fh.setFormatter(
                logging.Formatter("%(asctime)s [%(levelname)s] %(message)s")
            )
            log.addHandler(fh)
        except OSError:
            pass
    sh = logging.StreamHandler()
    sh.setLevel(logging.INFO)
    sh.setFormatter(logging.Formatter("[%(levelname)s] %(message)s"))
    log.addHandler(sh)
    return log

# processing/dynamic/operations/_logging.py
"""Re-export shared stage logger for transient operations."""

from processing.common.stage_logging import init_stage_logger

__all__ = ["init_stage_logger"]

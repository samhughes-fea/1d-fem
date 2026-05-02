# tests/test_nonlinear_newton_residual_metric.py
"""Tests for Newton condensed-residual convergence metric and settings parsing."""

import os
import tempfile

from pre_processing.parsing.simulation_settings_parser import parse_simulation_settings
from simulation_runner.static.nonlinear_static_simulation import (
    newton_condensed_residual_converged,
)


class TestNewtonCondensedResidualConverged:
    def test_absolute_only(self):
        assert newton_condensed_residual_converged(1e-7, 1e-6, None, 1.0)
        assert not newton_condensed_residual_converged(1e-5, 1e-6, None, 1.0)

    def test_relative_adds_to_atol(self):
        # 4 <= 1 + 0.1 * 30
        assert newton_condensed_residual_converged(4.0, 1.0, 0.1, 30.0)
        # 5 > 1 + 0.1 * 30
        assert not newton_condensed_residual_converged(5.0, 1.0, 0.1, 30.0)

    def test_tiny_ref_scale(self):
        assert newton_condensed_residual_converged(1e-12, 1e-10, 0.0, 0.0)


def test_newton_parser_relative_keys():
    with tempfile.NamedTemporaryFile(mode="w", suffix=".txt", delete=False) as f:
        f.write(
            """[Simulation]
[Type]
static_nonlinear

[Newton]
tolerance = 0.1
relative_tolerance = 0.0001
relative_reference = external_force
"""
        )
        temp_path = f.name
    try:
        settings = parse_simulation_settings(temp_path)
        assert settings["newton"]["tolerance"] == 0.1
        assert settings["newton"]["relative_tolerance"] == 0.0001
        assert settings["newton"]["relative_reference"] == "external_force"
    finally:
        os.unlink(temp_path)


def test_newton_defaults_include_relative():
    with tempfile.NamedTemporaryFile(mode="w", suffix=".txt", delete=False) as f:
        f.write(
            """[Simulation]
[Type]
static_nonlinear
"""
        )
        temp_path = f.name
    try:
        settings = parse_simulation_settings(temp_path)
        assert settings["newton"].get("relative_tolerance") is None
        assert settings["newton"].get("relative_reference") == "first_residual"
    finally:
        os.unlink(temp_path)

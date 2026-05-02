"""Parser [Nonlinear] section, load increment schedule, and runner flags (no full job run)."""

import os
import tempfile

import numpy as np

from pre_processing.parsing.simulation_settings_parser import parse_simulation_settings
from simulation_runner.static.nonlinear_static_simulation import NonlinearStaticSimulationRunner


def test_parser_nonlinear_section_defaults_merge():
    with tempfile.NamedTemporaryFile(mode="w", suffix=".txt", delete=False) as f:
        f.write(
            """[Simulation]
[Type]
static_nonlinear
"""
        )
        path = f.name
    try:
        s = parse_simulation_settings(path)
        assert s["nonlinear"]["num_increments"] == 1
        assert s["nonlinear"]["line_search"] is False
        assert s["nonlinear"]["load_factors"] is None
        assert s["nonlinear"]["line_search_max_backtracks"] == 6
        assert s["nonlinear"]["line_search_shrink"] == 0.5
    finally:
        os.unlink(path)


def test_parser_nonlinear_overrides():
    with tempfile.NamedTemporaryFile(mode="w", suffix=".txt", delete=False) as f:
        f.write(
            """[Simulation]
[Type]
static_nonlinear

[Nonlinear]
num_increments = 3
line_search = true
line_search_max_backtracks = 4
line_search_shrink = 0.25
load_factors = 0.3, 0.7, 1.0
"""
        )
        path = f.name
    try:
        s = parse_simulation_settings(path)
        assert s["nonlinear"]["num_increments"] == 3
        assert s["nonlinear"]["line_search"] is True
        assert s["nonlinear"]["line_search_max_backtracks"] == 4
        assert s["nonlinear"]["line_search_shrink"] == 0.25
        assert s["nonlinear"]["load_factors"] == [0.3, 0.7, 1.0]
    finally:
        os.unlink(path)


def test_compute_load_factors_uniform_schedule():
    """Linspace(1/n, 1, n) when load_factors omitted."""
    r = object.__new__(NonlinearStaticSimulationRunner)
    r.nonlinear_load_factors = None
    r.nonlinear_num_increments = 4
    out = NonlinearStaticSimulationRunner._compute_load_factors(r)
    np.testing.assert_allclose(out, [0.25, 0.5, 0.75, 1.0])


def test_compute_load_factors_explicit_list():
    r = object.__new__(NonlinearStaticSimulationRunner)
    r.nonlinear_load_factors = [0.1, 0.5, 1.0]
    r.nonlinear_num_increments = 99
    out = NonlinearStaticSimulationRunner._compute_load_factors(r)
    np.testing.assert_allclose(out, [0.1, 0.5, 1.0])


def test_newton_condensed_helper_matches_runner_threshold_semantics():
    """Sanity: atol + rtol * ref matches documented NR residual gate."""
    from simulation_runner.static.nonlinear_static_simulation import newton_condensed_residual_converged

    assert newton_condensed_residual_converged(0.05, 0.1, 1e-3, 100.0)
    # threshold = 0.1 + 1e-3 * 100 = 0.2 — strictly above boundary should fail
    assert not newton_condensed_residual_converged(0.21, 0.1, 1e-3, 100.0)

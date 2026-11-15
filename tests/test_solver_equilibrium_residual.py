"""Equilibrium residual of a scaled condensed solve.

Jacobi scaling replaces ``K`` and ``F`` and the solve returns the unscaled
displacement. The residual must be ``||K U - F||`` on that unscaled system,
not ``||K_scaled U - F_scaled||``.
"""

import csv
from pathlib import Path

import numpy as np
import pytest
import scipy.sparse as sp

from onedfem.processing.static.operations.solver import (
    SolveCondensedSystem,
    _RELATIVE_RESIDUAL_WARN,
    _relative_residual_is_high,
    _row_col_scale,
)


def test_solve_condensed_system_scaled_system_reports_equilibrium_residual(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """A scaled solve reports the unscaled residual, not the mixed one.

    The mixed residual is the still-scaled stiffness and force applied to the
    unscaled displacement. It sits above the warning tolerance. The equilibrium
    residual does not, and it is the value written to the solve history.
    """
    # Given a condensed system whose diagonal makes Jacobi scaling nontrivial
    monkeypatch.delenv("FEM_DISABLE_SCALING", raising=False)
    stiffness = np.array(
        [
            [1.0e4, 0.2, 0.0],
            [0.2, 1.0, -0.01],
            [0.0, -0.01, 1.0e-2],
        ],
        dtype=np.float64,
    )
    displacement = np.array([0.01, -0.2, 3.0], dtype=np.float64)
    force = stiffness @ displacement
    condensed = sp.csr_matrix(stiffness)
    job_results = tmp_path / "primary_results"
    job_results.mkdir()
    solver = SolveCondensedSystem(
        K_cond=condensed,
        F_cond=force,
        solver_name="direct",
        job_results_dir=str(job_results),
        preconditioner=None,
    )

    # When the condensed system is solved
    solved = solver.solve()

    # Then the returned displacement is the physical solution
    assert solved is not None
    assert solved == pytest.approx(displacement, rel=1e-9, abs=1e-12)

    scaled_stiffness, scaled_force, _scale = _row_col_scale(condensed.copy(), force.copy())
    mixed = float(np.linalg.norm(scaled_stiffness @ solved - scaled_force))
    equilibrium = float(np.linalg.norm(stiffness @ solved - force))
    relative = equilibrium / float(np.linalg.norm(force))

    # The mixed residual trips the warning tolerance; the equilibrium residual does not.
    assert mixed > _RELATIVE_RESIDUAL_WARN
    assert mixed > 1.0
    assert equilibrium < _RELATIVE_RESIDUAL_WARN
    assert relative < _RELATIVE_RESIDUAL_WARN

    history = tmp_path / "logs" / "inner_solve_history.csv"
    with open(history, newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    assert len(rows) == 1
    reported = float(rows[0]["condensed_residual_norm"])
    reported_relative = float(rows[0]["final_relative_residual"])

    # Reporting the mixed residual fails; the equilibrium residual passes.
    # Dense and sparse products of this stiffness differ by a few ulps of the
    # 1e4 diagonal entry, so the match is at 1e-12 rather than exact.
    assert reported == pytest.approx(equilibrium)
    assert reported < _RELATIVE_RESIDUAL_WARN
    assert reported != pytest.approx(mixed)
    assert reported_relative == pytest.approx(relative)
    assert solver.diagnostics["equilibrium_residual"] == pytest.approx(reported)
    assert solver.diagnostics["relative_residual"] == pytest.approx(reported_relative)

    log_text = (tmp_path / "logs" / "SolveCondensedSystem.log").read_text(encoding="utf-8")
    assert "Residual is high" not in log_text
    assert "relative" in log_text


def test_relative_residual_is_high_mixed_residual_fails_equilibrium_passes() -> None:
    """The warning gate rejects the mixed residual and accepts the equilibrium one.

    On the RM1 1.2 at U_out beam the mixed residual is 0.216. The equilibrium
    residual is 5.96e-6 N, relative 1.42e-9 against a 4.19e3 N force, so an
    absolute check at 1e-6 still warns and the relative check does not. A
    textbook cantilever trips the same absolute check at 3.12; its equilibrium
    residual is 1.3e-12.
    """
    # Given the residuals of a solved beam and of the mixed-scale check
    mixed_rm1 = 0.216
    mixed_cantilever = 3.12
    equilibrium_rm1 = 5.96e-6
    force_rm1 = 4.19e3
    relative_rm1 = equilibrium_rm1 / force_rm1
    equilibrium_cantilever = 1.3e-12

    # When the warning gate classifies them
    mixed_rm1_is_high = _relative_residual_is_high(mixed_rm1)
    mixed_cantilever_is_high = _relative_residual_is_high(mixed_cantilever)
    equilibrium_rm1_is_high = _relative_residual_is_high(relative_rm1)
    equilibrium_cantilever_is_high = _relative_residual_is_high(equilibrium_cantilever)

    # Then the mixed residuals fail and the equilibrium residuals pass
    assert equilibrium_rm1 > _RELATIVE_RESIDUAL_WARN
    assert relative_rm1 == pytest.approx(1.42e-9, rel=1e-2)
    assert mixed_rm1_is_high
    assert mixed_cantilever_is_high
    assert not equilibrium_rm1_is_high
    assert not equilibrium_cantilever_is_high

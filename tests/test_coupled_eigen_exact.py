"""1d-fem's eigen solver on a fully coupled anisotropic section, against the exact prismatic beam
(``verification/eigen/coupled_prismatic_exact.py``)."""

from __future__ import annotations

import os
import sys

import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))
sys.path.insert(0, HERE)

from verification.eigen.coupled_prismatic_exact import prismatic_frequencies  # noqa: E402

pytest.importorskip("anisotropic_section_model.section_analysis")
import test_section_matrix_job as tm  # noqa: E402


def test_exact_reference_reduces_to_the_uncoupled_cantilever():
    E, G, b, h, L, rho = 70e9, 26e9, 0.05, 0.1, 8.0, 2700.0
    A, Iy, Iz = b * h, b * h ** 3 / 12, h * b ** 3 / 12
    K = np.diag([E * A, E * Iy, E * Iz, 1e3 * G * A, 1e3 * G * A, G * 1e-6])
    M = np.diag([rho * A] * 3 + [rho * 1e-9, 1e-9 * rho * Iy, 1e-9 * rho * Iz])
    f = prismatic_frequencies(K, M, L, 2.0)
    beta = 1.8751040687
    exact = sorted(beta ** 2 / (2 * np.pi) * np.sqrt(E * I / (rho * A * L ** 4)) for I in (Iy, Iz))
    np.testing.assert_allclose(f[:2], exact, rtol=1e-6)


def test_coupled_section_converges_to_the_exact_frequencies():
    """Six modes of an unbalanced, off-centre composite box: bending, torsion and their couplings."""
    from anisotropic_section_model.section_analysis import SectionAnalysis

    an = tm.common.build_analysis(tm.SECTION)
    M = SectionAnalysis.from_thin_walled(an).mass_matrix()
    errors = []
    for n in (6, 12, 24):
        fe = tm._run(f"eigen_exact_{n}", "eigen", table="engine", n_el=n)
        exact = prismatic_frequencies(an.stiffness, M, tm.LENGTH, 1.05 * fe.max())
        assert exact.size >= fe.size
        errors.append(np.abs(fe / exact[:fe.size] - 1.0))
    errors = np.array(errors)
    assert np.all(errors[1:] < errors[:-1] + 1e-6), errors      # converges with the mesh
    assert errors[-1].max() < 0.01, errors[-1]
    assert errors[-1][:2].max() < 5e-4, errors[-1]

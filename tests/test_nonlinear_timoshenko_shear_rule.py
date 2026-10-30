"""NonlinearTimoshenkoBeamElement3D integrates transverse shear on the same reduced rule as K_0.

Before the fix ``F_int`` and ``K_delta`` used the full Gauss rule for every strain while
``K_0`` used one point for shear. The linear two-node element then locks in ``F_int`` but
not in its tangent, Newton cannot reach a consistent equilibrium, and a slender cantilever
came out at -16.6 times its exact tip deflection.
"""

from __future__ import annotations

import logging
import os
import sys
import tempfile

import numpy as np
import pytest

sys.path.insert(0, os.path.dirname(__file__))

from test_point_moments import _exact_tip, _line, _run  # noqa: E402
from test_warping_nl_linear_material_agreement import _job_root, _section_gamma  # noqa: E402


def _elements(warping: bool, L: float = 0.5):
    from pre_processing.element_library.linear.beam.first_order_shear_deformation_theory.timoshenko.linear_timoshenko_3D import (
        LinearTimoshenkoBeamElement3D,
    )
    from pre_processing.element_library.nonlinear.timoshenko.nonlinear_timoshenko_3D import (
        NonlinearTimoshenkoBeamElement3D,
    )

    orders = {k: np.array([v]) for k, v in dict(axial=3, bending_y=3, bending_z=3, shear_y=2, shear_z=2,
                                                 torsion=3, load=3).items()}
    ed = {"ids": np.array([0]), "connectivity": np.array([[0, 1]]),
          "types": np.array(["LinearTimoshenkoBeamElement3D"]), "integration_orders": orders}
    if warping:
        ed["warping"] = np.array([1], dtype=np.int8)
    kw = dict(element_id=0, grid_dictionary={"coordinates": np.array([[0.0, 0, 0], [L, 0, 0]])},
              material_dictionary={"E": np.array([2.1e11]), "G": np.array([8.1e10]), "nu": np.array([0.3]),
                                   "rho": np.array([0.0])},
              section_dictionary=_section_gamma(1e-7), point_load_array=np.empty((0, 9)),
              distributed_load_array=np.empty((0, 9)), job_results_dir=_job_root(tempfile.mkdtemp()))
    lin = LinearTimoshenkoBeamElement3D(element_dictionary=ed, **kw)
    nl = NonlinearTimoshenkoBeamElement3D(
        element_dictionary=dict(ed, types=np.array(["NonlinearTimoshenkoBeamElement3D"])), **kw)
    return lin, nl


@pytest.mark.parametrize("warping", [False, True])
def test_internal_force_is_consistent_with_the_linear_stiffness(warping):
    logging.disable(logging.WARNING)
    lin, nl = _elements(warping)
    Ke = lin.element_stiffness_matrix().K_e
    n = Ke.shape[0]
    h = 1e-7
    J = np.column_stack([(nl.internal_force_vector(e) - nl.internal_force_vector(-e)) / (2 * h)
                         for e in np.eye(n) * h])
    scale = np.abs(Ke).max()
    assert np.abs(J - Ke).max() < 1e-10 * scale
    assert np.abs(nl.tangent_stiffness_matrix(np.zeros(n)) - Ke).max() < 1e-12 * scale


def test_tangent_is_the_derivative_of_the_internal_force():
    logging.disable(logging.WARNING)
    _, nl = _elements(False)
    rng = np.random.default_rng(3)
    U = rng.normal(size=12) * 1e-3
    h = 1e-8
    J = np.column_stack([(nl.internal_force_vector(U + e) - nl.internal_force_vector(U - e)) / (2 * h)
                         for e in np.eye(12) * h])
    K = nl.tangent_stiffness_matrix(U)
    assert np.abs(J - K).max() < 1e-5 * np.abs(K).max()


def test_slender_cantilever_matches_the_linear_element(tmp_path):
    L = 2.0
    row = [L, 0, 0, 0, 0, 1.0, 0, 0, 0]
    U_lin = _run(tmp_path, "LinearTimoshenkoBeamElement3D", _line(8, L), [row])
    U_nl = _run(tmp_path, "NonlinearTimoshenkoBeamElement3D", _line(8, L), [row], nonlinear=True)
    assert U_nl[-1, 2] == pytest.approx(U_lin[-1, 2], rel=1e-5)
    assert U_nl[-1, 2] == pytest.approx(_exact_tip(row[3:], L)[2], rel=0.01)

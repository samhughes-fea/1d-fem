"""Native unshearable GEBT: linear limit, vs TL EB, elastica, factory increment."""

from __future__ import annotations

import shutil

import numpy as np
import pytest

from onedfem.pre_processing.element_library.element_factory import ElementFactory
from tests.test_element_architecture_layout import _make_element, _newton_increments, _ts_dicts, _job_dirs


def test_gebt_unshearable_zero_state_matches_linear_eb():
    root, job = _job_dirs()
    try:
        gebt = _make_element("GEBTUnshearableBeamElement3D", job)
        lin = _make_element("LinearEulerBernoulliBeamElement3D", job)
        U0 = np.zeros(12)
        K = gebt.tangent_stiffness_matrix(U0)
        np.testing.assert_allclose(K, lin.element_stiffness_matrix().K_e, rtol=1e-8, atol=1e-6)
        np.testing.assert_allclose(gebt.internal_force_vector(U0), 0.0, atol=1e-8)
    finally:
        shutil.rmtree(root, ignore_errors=True)


def test_gebt_unshearable_moderate_load_vs_tl_eb():
    root, job = _job_dirs()
    try:
        P = -80.0
        F = np.zeros(12)
        F[7] = P
        gebt = _make_element("GEBTUnshearableBeamElement3D", job)
        tl = _make_element("NonlinearEulerBernoulliBeamElement3D", job)
        U_g, _ = _newton_increments(gebt, F, n_inc=1, n_iter=20, tol=1e-6)
        U_t, _ = _newton_increments(tl, F, n_inc=1, n_iter=20, tol=1e-6)
        np.testing.assert_allclose(U_g[7], U_t[7], rtol=0.05, atol=1e-4)
        assert abs(U_g[7] - U_t[7]) > 1e-16
    finally:
        shutil.rmtree(root, ignore_errors=True)


def test_gebt_unshearable_bisshopp_drucker_cantilever():
    """P L^2 / EI = 1: tip drop 0.30172 L, inbound 0.05643 L (many elements)."""
    from tests.test_corotational_frames import E, RECT, _line, _newton, _tip, Mesh

    L, n = 1.0, 16
    pts = _line([0, 0, 0], [1, 0, 0], L, n)
    root, job = _job_dirs()
    try:
        mesh = Mesh("GEBTUnshearableBeamElement3D", pts, props=RECT)
        for el in mesh.elements:
            el.job_results_dir = job
        EI = E * RECT["I_y"]
        U = _newton(mesh, _tip(mesh, n, [0, 0, -EI / L**2, 0, 0, 0]), list(range(6)), 8)
        tip = pts[-1] + U.reshape(-1, 6)[-1, :3]
        assert tip[2] == pytest.approx(-0.30172 * L, rel=5e-2)
        assert L - tip[0] == pytest.approx(0.05643 * L, rel=8e-2)
    finally:
        shutil.rmtree(root, ignore_errors=True)


def test_gebt_unshearable_factory_static_increment():
    root, job = _job_dirs()
    factory = ElementFactory(job_results_dir=job)
    ed, gd, mat, sec = _ts_dicts("GEBTUnshearableBeamElement3D")
    ed["warping"] = np.zeros(1, dtype=np.int8)
    try:
        els = factory.create_elements_batch(
            element_ids=np.array([0]),
            element_dictionary=ed,
            grid_dictionary=gd,
            material_dictionary=mat,
            section_dictionary=sec,
            point_load_array=np.empty((0, 9)),
            distributed_load_array=np.empty((0, 9)),
            enable_parallel=False,
        )
        el = els[0]
        el.logger_operator = None
        U = np.zeros(12)
        U[7] = 1e-5
        K = np.asarray(el.tangent_stiffness_matrix(U))
        F = np.asarray(el.internal_force_vector(U))
        assert np.all(np.isfinite(K)) and np.all(np.isfinite(F))
        assert K.shape == (12, 12)
    finally:
        shutil.rmtree(root, ignore_errors=True)

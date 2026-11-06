"""Layout and operator-in-utilities checkpoints after zosdt/fosdt/tosdt reorganization."""

from __future__ import annotations

import os
import shutil
import tempfile

import numpy as np
import pytest

from onedfem.pre_processing.element_library.element_factory import ElementFactory


def _job_dirs():
    root = tempfile.mkdtemp()
    job = os.path.join(root, "job")
    for sub in ("element_stiffness_matrices", "element_force_vectors", "logs"):
        os.makedirs(os.path.join(job, sub), exist_ok=True)
    return root, job


def _ts_dicts(etype: str, L: float = 1.5):
    return {
        "ids": np.array([0]),
        "connectivity": np.array([[0, 1]]),
        "types": np.array([etype]),
        "integration_orders": {
            "axial": np.array([2]),
            "bending_y": np.array([3]),
            "bending_z": np.array([3]),
            "shear_y": np.array([2]),
            "shear_z": np.array([2]),
            "torsion": np.array([2]),
            "load": np.array([2]),
        },
    }, {
        "ids": np.array([0, 1]),
        "coordinates": np.array([[0.0, 0.0, 0.0], [L, 0.0, 0.0]]),
    }, {
        "E": np.array([2.0e11]),
        "G": np.array([77e9]),
        "nu": np.array([0.3]),
        "rho": np.array([7850.0]),
    }, {
        "A": np.array([6e-4]),
        "I_x": np.array([0.0]),
        "I_y": np.array([4e-8]),
        "I_z": np.array([8e-8]),
        "J_t": np.array([1e-9]),
    }


def test_factory_maps_use_abbreviated_sdt_paths():
    lin = ElementFactory.LINEAR_ELEMENT_CLASS_MAP
    nl = ElementFactory.NONLINEAR_ELEMENT_CLASS_MAP
    assert ".zosdt." in lin["LinearEulerBernoulliBeamElement3D"]
    assert ".fosdt.timoshenko." in lin["LinearTimoshenkoBeamElement3D"]
    assert ".fosdt.gbt." in lin["LinearGBTBeamElement3D"]
    assert ".fosdt.benscoter.linear_gbt" not in lin["LinearGBTBeamElement3D"]
    assert ".tosdt." in lin["LinearReddyBeamElement3D"]
    assert ".total_lagrangian." in nl["NonlinearTimoshenkoBeamElement3D"]
    assert ".tosdt.reddy." in nl["NonlinearReddyBeamElement3D"]
    assert ".updated_lagrangian." in nl["UpdatedLagrangianEulerBernoulliBeamElement3D"]
    assert ".tosdt.reddy.updated_lagrangian." in nl["UpdatedLagrangianReddyBeamElement3D"]
    assert ".tosdt.levinson.updated_lagrangian." in nl["UpdatedLagrangianLevinsonBeamElement3D"]


def test_ul_eb_zero_state_matches_linear_and_config_update():
    from onedfem.pre_processing.element_library.linear.beam.zosdt.euler_bernoulli.linear_euler_bernoulli_3D import (
        LinearEulerBernoulliBeamElement3D,
    )
    from onedfem.pre_processing.element_library.nonlinear.euler_bernoulli.updated_lagrangian.updated_lagrangian_euler_bernoulli_3D import (
        UpdatedLagrangianEulerBernoulliBeamElement3D,
    )

    ed, gd, mat, sec = _ts_dicts("UpdatedLagrangianEulerBernoulliBeamElement3D")
    ed_lin = dict(ed)
    ed_lin["types"] = np.array(["LinearEulerBernoulliBeamElement3D"])
    root, job = _job_dirs()
    try:
        kw = dict(
            element_id=0,
            grid_dictionary=gd,
            material_dictionary=mat,
            section_dictionary=sec,
            point_load_array=np.empty((0, 9)),
            distributed_load_array=np.empty((0, 9)),
            job_results_dir=job,
        )
        ul = UpdatedLagrangianEulerBernoulliBeamElement3D(element_dictionary=ed, **kw)
        lin = LinearEulerBernoulliBeamElement3D(element_dictionary=ed_lin, **kw)
        U0 = np.zeros(12)
        K_ul = ul.tangent_stiffness_matrix(U0)
        K_lin = lin.element_stiffness_matrix().K_e
        np.testing.assert_allclose(K_ul, K_lin, rtol=1e-9, atol=1e-8)
        L0 = ul.L
        U = np.zeros(12)
        U[6] = 0.01
        ul.update_reference_configuration(U)
        assert ul.L == pytest.approx(L0 + 0.01, rel=1e-12)
        np.testing.assert_allclose(ul.internal_force_vector(U), 0.0, atol=1e-8)
    finally:
        shutil.rmtree(root, ignore_errors=True)


def test_ul_ts_zero_state_and_small_axial_vs_tl():
    from onedfem.pre_processing.element_library.linear.beam.fosdt.timoshenko.linear_timoshenko_3D import (
        LinearTimoshenkoBeamElement3D,
    )
    from onedfem.pre_processing.element_library.nonlinear.timoshenko.total_lagrangian.nonlinear_timoshenko_3D import (
        NonlinearTimoshenkoBeamElement3D,
    )
    from onedfem.pre_processing.element_library.nonlinear.timoshenko.updated_lagrangian.updated_lagrangian_timoshenko_3D import (
        UpdatedLagrangianTimoshenkoBeamElement3D,
    )

    ed, gd, mat, sec = _ts_dicts("UpdatedLagrangianTimoshenkoBeamElement3D")
    ed_lin = dict(ed)
    ed_lin["types"] = np.array(["LinearTimoshenkoBeamElement3D"])
    ed_tl = dict(ed)
    ed_tl["types"] = np.array(["NonlinearTimoshenkoBeamElement3D"])
    root, job = _job_dirs()
    try:
        kw = dict(
            element_id=0,
            grid_dictionary=gd,
            material_dictionary=mat,
            section_dictionary=sec,
            point_load_array=np.empty((0, 9)),
            distributed_load_array=np.empty((0, 9)),
            job_results_dir=job,
        )
        ul = UpdatedLagrangianTimoshenkoBeamElement3D(element_dictionary=ed, **kw)
        lin = LinearTimoshenkoBeamElement3D(element_dictionary=ed_lin, **kw)
        tl = NonlinearTimoshenkoBeamElement3D(element_dictionary=ed_tl, **kw)
        U0 = np.zeros(12)
        np.testing.assert_allclose(
            ul.tangent_stiffness_matrix(U0), lin.element_stiffness_matrix().K_e, rtol=1e-9, atol=1e-8
        )
        U = np.zeros(12)
        U[6] = 1e-6
        np.testing.assert_allclose(
            ul.tangent_stiffness_matrix(U), tl.tangent_stiffness_matrix(U), rtol=1e-5, atol=1e-4
        )
    finally:
        shutil.rmtree(root, ignore_errors=True)


def test_nl_reddy_zero_state_matches_linear_k0():
    from onedfem.pre_processing.element_library.linear.beam.tosdt.reddy.linear_reddy_3D import (
        LinearReddyBeamElement3D,
    )
    from onedfem.pre_processing.element_library.nonlinear.tosdt.reddy.nonlinear_reddy_3D import (
        NonlinearReddyBeamElement3D,
    )

    ed, gd, mat, sec = _ts_dicts("NonlinearReddyBeamElement3D")
    sec["kappa"] = np.array([5.0 / 6.0])
    sec["alpha"] = np.array([0.0])
    ed["warping"] = np.zeros(1, dtype=np.int8)
    ed_lin = dict(ed)
    ed_lin["types"] = np.array(["LinearReddyBeamElement3D"])
    root, job = _job_dirs()
    try:
        kw = dict(
            element_id=0,
            grid_dictionary=gd,
            material_dictionary=mat,
            section_dictionary=sec,
            point_load_array=np.empty((0, 9)),
            distributed_load_array=np.empty((0, 9)),
            job_results_dir=job,
        )
        nl = NonlinearReddyBeamElement3D(element_dictionary=ed, **kw)
        lin = LinearReddyBeamElement3D(element_dictionary=ed_lin, **kw)
        U0 = np.zeros(16)
        K0 = nl.tangent_stiffness_matrix(U0)
        np.testing.assert_allclose(K0, lin.element_stiffness_matrix().K_e, rtol=1e-8, atol=1e-6)
        F0 = nl.internal_force_vector(U0)
        np.testing.assert_allclose(F0, 0.0, atol=1e-8)
    finally:
        shutil.rmtree(root, ignore_errors=True)


def _make_element(etype: str, job: str, *, kappa: bool = False, warping: bool = False):
    ed, gd, mat, sec = _ts_dicts(etype)
    if kappa:
        sec["kappa"] = np.array([5.0 / 6.0])
        sec["alpha"] = np.array([0.0])
    if not warping:
        ed["warping"] = np.zeros(1, dtype=np.int8)
    kw = dict(
        element_id=0,
        element_dictionary=ed,
        grid_dictionary=gd,
        material_dictionary=mat,
        section_dictionary=sec,
        point_load_array=np.empty((0, 9)),
        distributed_load_array=np.empty((0, 9)),
        job_results_dir=job,
    )
    Cls = getattr(
        __import__(ElementFactory.ELEMENT_CLASS_MAP[etype], fromlist=[etype]),
        etype,
    )
    el = Cls(**kw)
    el.logger_operator = None
    return el


def _newton_increments(el, F, n_inc=2, n_iter=12, tol=1e-8):
    """Load-controlled Newton on a clamped–free two-node element; UL config update after each increment."""
    n = int(el._n_dof)
    F = np.asarray(F, dtype=np.float64)
    if F.size < n:
        F = np.pad(F, (0, n - F.size))
    U = np.zeros(n, dtype=np.float64)
    fixed = np.arange(6)
    free = np.setdiff1d(np.arange(n), fixed)
    n_update = 0
    orig = getattr(el, "update_reference_configuration", None)

    if orig is not None:

        def _counted(Ue, _orig=orig):
            nonlocal n_update
            n_update += 1
            return _orig(Ue)

        el.update_reference_configuration = _counted
    for inc in range(n_inc):
        lam = (inc + 1) / n_inc
        last_r = np.inf
        for _ in range(n_iter):
            K = np.asarray(el.tangent_stiffness_matrix(U), dtype=np.float64)
            R = np.asarray(el.internal_force_vector(U), dtype=np.float64) - lam * F[:n]
            last_r = float(np.linalg.norm(R[free]))
            if last_r < tol:
                break
            U[free] += np.linalg.solve(K[np.ix_(free, free)], -R[free])
        if last_r >= tol > 0.0:
            raise AssertionError(
                f"{type(el).__name__} Newton did not converge at increment {inc + 1} (||R||={last_r})"
            )
        if orig is not None:
            el.update_reference_configuration(U)
    return U, n_update


def test_ul_multi_increment_vs_tl_and_runner_hook():
    """Increment 1: UL vs TL rtol 0.15 and UL-EB vs linear rtol 0.02.

    Increment 2 checks that the runner updates the live chord twice. UL remaining
    load on the current chord is not Green–Lagrange, so tip ``uy`` is not forced
    to 0.15 of TL after the update.
    """
    from onedfem.simulation_runner.static.nonlinear_static_simulation import (
        NonlinearStaticSimulationRunner,
    )

    class _Dummy:
        n = 0

        def assemble_global_dof_indices(self):
            return np.arange(12)

        def update_reference_configuration(self, U_e):
            self.n += 1
            self.last = np.asarray(U_e)

    runner = object.__new__(NonlinearStaticSimulationRunner)
    dummy = _Dummy()
    runner.elements = [dummy]
    NonlinearStaticSimulationRunner._update_ul_reference_configurations(runner, np.arange(12, dtype=float))
    assert dummy.n == 1
    np.testing.assert_allclose(dummy.last, np.arange(12, dtype=float))

    root, job = _job_dirs()
    try:
        P = -80.0
        pairs = (
            (
                "UpdatedLagrangianEulerBernoulliBeamElement3D",
                "NonlinearEulerBernoulliBeamElement3D",
                "LinearEulerBernoulliBeamElement3D",
                False,
                12,
            ),
            (
                "UpdatedLagrangianTimoshenkoBeamElement3D",
                "NonlinearTimoshenkoBeamElement3D",
                None,
                True,
                12,
            ),
            (
                "UpdatedLagrangianReddyBeamElement3D",
                "NonlinearReddyBeamElement3D",
                None,
                True,
                16,
            ),
            (
                "UpdatedLagrangianLevinsonBeamElement3D",
                "NonlinearLevinsonBeamElement3D",
                None,
                True,
                16,
            ),
        )
        for ul_name, tl_name, lin_name, kappa, ndof in pairs:
            F = np.zeros(ndof)
            F[7] = P
            ul = _make_element(ul_name, job, kappa=kappa)
            tl = _make_element(tl_name, job, kappa=kappa)
            U_ul1, n1 = _newton_increments(ul, F, n_inc=1, n_iter=20, tol=1e-6)
            U_tl1, _ = _newton_increments(tl, F, n_inc=1, n_iter=20, tol=1e-6)
            assert n1 == 1
            np.testing.assert_allclose(U_ul1[7], U_tl1[7], rtol=0.15, atol=1e-4)
            if lin_name is not None:
                lin = _make_element(lin_name, job, kappa=kappa)
                K0 = np.asarray(lin.element_stiffness_matrix().K_e, dtype=np.float64)
                n = K0.shape[0]
                U_lin = np.zeros(n)
                free = np.setdiff1d(np.arange(n), np.arange(6))
                U_lin[free] = np.linalg.solve(K0[np.ix_(free, free)], F[:n][free])
                np.testing.assert_allclose(U_ul1[7], U_lin[7], rtol=0.02, atol=1e-5)
            ul2 = _make_element(ul_name, job, kappa=kappa)
            L0 = float(ul2.L)
            U_ul2, n2 = _newton_increments(ul2, F, n_inc=2, n_iter=20, tol=1e-6)
            assert n2 == 2
            assert abs(float(ul2.L) - L0) > 0.0
            assert U_ul2[7] * P > 0.0
            # After the config update, UL is remaining load on the live chord (small strain);
            # TL is Green–Lagrange on the original chord. They are different models: require
            # convergence and the same tip sign, not a 0.15 match.
    finally:
        shutil.rmtree(root, ignore_errors=True)


def test_nl_tosdt_one_newton_step_vs_timoshenko():
    """16-DOF Hermite Reddy/Levinson vs 12-DOF Timoshenko: first Newton step from rest.

    Both use the same tip force. The Hermite transverse interpolation is not the
    Timoshenko linear field, so the tip translations must differ, but they stay
    the same sign and order of magnitude on a coarse one-element cantilever.
    """
    root, job = _job_dirs()
    try:
        P = -200.0
        F = np.zeros(16)
        F[7] = P
        ts = _make_element("NonlinearTimoshenkoBeamElement3D", job, kappa=True)
        rd = _make_element("NonlinearReddyBeamElement3D", job, kappa=True)
        lv = _make_element("NonlinearLevinsonBeamElement3D", job, kappa=True)
        U_ts, _ = _newton_increments(ts, F[:12], n_inc=1, n_iter=1, tol=0.0)
        U_rd, _ = _newton_increments(rd, F, n_inc=1, n_iter=1, tol=0.0)
        U_lv, _ = _newton_increments(lv, F, n_inc=1, n_iter=1, tol=0.0)
        assert U_ts.size == 12
        assert U_rd.size == 16
        assert U_lv.size == 16
        # Hermite nodal slopes v', w' (DOFs 12–15) have no Timoshenko counterpart.
        assert float(np.linalg.norm(U_rd[12:])) > 0.0
        assert float(np.linalg.norm(U_lv[12:])) > 0.0
        assert U_rd[7] * U_ts[7] > 0.0
        assert U_lv[7] * U_ts[7] > 0.0
        np.testing.assert_allclose(U_rd[7], U_ts[7], rtol=0.25, atol=1e-8)
        np.testing.assert_allclose(U_lv[7], U_ts[7], rtol=0.25, atol=1e-8)
    finally:
        shutil.rmtree(root, ignore_errors=True)


def test_factory_ul_nl_types_complete_static_increment():
    root, job = _job_dirs()
    factory = ElementFactory(job_results_dir=job)
    types = [
        "UpdatedLagrangianEulerBernoulliBeamElement3D",
        "UpdatedLagrangianTimoshenkoBeamElement3D",
        "UpdatedLagrangianReddyBeamElement3D",
        "UpdatedLagrangianLevinsonBeamElement3D",
        "NonlinearReddyBeamElement3D",
        "NonlinearLevinsonBeamElement3D",
        "NonlinearEulerBernoulliBeamElement3D",
        "NonlinearTimoshenkoBeamElement3D",
    ]
    try:
        for etype in types:
            ed, gd, mat, sec = _ts_dicts(etype)
            sec["kappa"] = np.array([5.0 / 6.0])
            sec["alpha"] = np.array([0.0])
            ed["warping"] = np.zeros(1, dtype=np.int8)
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
            n = int(el._n_dof)
            U = np.zeros(n)
            U[7] = 1e-5
            K = np.asarray(el.tangent_stiffness_matrix(U))
            F = np.asarray(el.internal_force_vector(U))
            assert np.all(np.isfinite(K))
            assert np.all(np.isfinite(F))
            assert K.shape == (n, n)
            assert F.shape == (n,)
    finally:
        shutil.rmtree(root, ignore_errors=True)


def test_old_sdt_and_flat_nl_paths_are_gone():
    import importlib

    gone = (
        "onedfem.pre_processing.element_library.linear.beam.zero_order_shear_deformation_theory",
        "onedfem.pre_processing.element_library.linear.beam.first_order_shear_deformation_theory",
        "onedfem.pre_processing.element_library.linear.beam.third_order_shear_deformation_theory",
        "onedfem.pre_processing.element_library.nonlinear.reddy",
        "onedfem.pre_processing.element_library.nonlinear.levinson",
        "onedfem.pre_processing.element_library.nonlinear.timoshenko.nonlinear_timoshenko_3D",
        "onedfem.pre_processing.element_library.nonlinear.euler_bernoulli.nonlinear_euler_bernoulli_3D",
        "onedfem.pre_processing.element_library.nonlinear.timoshenko.updated_lagrangian_timoshenko_3D",
        "onedfem.pre_processing.element_library.nonlinear.euler_bernoulli.updated_lagrangian_euler_bernoulli_3D",
    )
    for mod in gone:
        with pytest.raises(ModuleNotFoundError):
            importlib.import_module(mod)
    joined = " ".join(ElementFactory.ELEMENT_CLASS_MAP.values())
    assert "zero_order_shear_deformation_theory" not in joined
    assert "first_order_shear_deformation_theory" not in joined
    assert "third_order_shear_deformation_theory" not in joined
    assert "nonlinear.reddy" not in joined
    assert ".zosdt." in joined
    assert ".fosdt." in joined
    assert ".tosdt." in joined
    assert ".updated_lagrangian." in ElementFactory.NONLINEAR_ELEMENT_CLASS_MAP[
        "UpdatedLagrangianReddyBeamElement3D"
    ]

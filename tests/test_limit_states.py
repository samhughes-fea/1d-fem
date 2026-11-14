"""Deflection, extreme-fibre stress, buckling and blade modes a caller reads back.

The static profiles are checked against the Euler-Bernoulli cantilever, which a
Hermite beam represents exactly. Buckling factors and blade frequencies are the
eigenpairs the existing solvers already computed: reading them does not solve
again and does not assemble a rotor. The mesh below is fixed; it is not refined
to chase a frequency or compliance gate.
"""

from __future__ import annotations

import os
import shutil
import tempfile

import numpy as np
import pytest

from onedfem.post_processing.limit_states import (
    blade_vibration_modes,
    buckling_result,
    deflection_profile,
    extreme_fibre_profile,
)

ELEMENT = "LinearEulerBernoulliBeamElement3D"


def _cantilever(n_elem: int, L: float, E: float, I: float, rho: float, tip_force):
    """Straight cantilever along x, root at the origin. ``tip_force`` is ``(Fx, Fy, Fz)``."""
    from onedfem.pre_processing.element_library.linear.beam.zosdt.euler_bernoulli.linear_euler_bernoulli_3D import (
        LinearEulerBernoulliBeamElement3D,
    )

    n_nodes = n_elem + 1
    dx = L / n_elem
    coordinates = np.array([[i * dx, 0.0, 0.0] for i in range(n_nodes)], dtype=np.float64)
    ids = np.arange(n_nodes, dtype=np.int32)
    section = {
        "A": np.full(n_elem, 0.01),
        "I_x": np.full(n_elem, 1.0e-9),
        "I_y": np.full(n_elem, I),
        "I_z": np.full(n_elem, I),
        "J_t": np.full(n_elem, 1.0e-9),
    }
    material = {
        "E": np.full(n_elem, E),
        "G": np.full(n_elem, E / 2.6),
        "nu": np.full(n_elem, 0.3),
        "rho": np.full(n_elem, rho),
    }
    element_dictionary = {
        "ids": np.arange(n_elem, dtype=np.int32),
        "connectivity": np.array([[i, i + 1] for i in range(n_elem)], dtype=np.int32),
        "types": np.array([ELEMENT] * n_elem),
        "integration_orders": {name: np.full(n_elem, 4, dtype=np.int64)
                               for name in ("axial", "bending_y", "bending_z", "shear_y", "shear_z", "torsion", "load")},
    }
    grid = {"coordinates": coordinates, "ids": ids}
    Fx, Fy, Fz = tip_force
    point_load = np.array([[L, 0.0, 0.0, Fx, Fy, Fz, 0.0, 0.0, 0.0]], dtype=np.float64)
    distributed = np.empty((0, 9))
    temp = tempfile.mkdtemp()
    job = os.path.join(temp, "job")
    os.makedirs(os.path.join(job, "element_stiffness_matrices"), exist_ok=True)
    os.makedirs(os.path.join(job, "element_force_vectors"), exist_ok=True)
    elements, stiffness, mass, element_objects, force_objects = [], [], [], [], []
    for eid in range(n_elem):
        element = LinearEulerBernoulliBeamElement3D(
            element_id=eid, element_dictionary=element_dictionary, grid_dictionary=grid,
            material_dictionary=material, section_dictionary=section, point_load_array=point_load,
            distributed_load_array=distributed, job_results_dir=job,
        )
        elements.append(element)
        element_objects.append(element.element_stiffness_matrix())
        force_objects.append(element.element_force_vector())
        stiffness.append(np.asarray(element_objects[-1].K_e, dtype=np.float64))
        mass.append(np.asarray(element.element_mass_matrix().M_e, dtype=np.float64))
    return {
        "elements": elements,
        "coordinates": coordinates,
        "stiffness": stiffness,
        "mass": mass,
        "element_objects": element_objects,
        "force_objects": force_objects,
        "temp": temp,
        "job": job,
        "grid": grid,
        "mesh": {"node_ids": ids, "coordinates": coordinates},
        "element_dictionary": element_dictionary,
        "material": material,
        "section": section,
        "point_load": point_load,
        "distributed": distributed,
    }


def _static_displacement(built) -> np.ndarray:
    elements = built["elements"]
    ndof = built["coordinates"].shape[0] * 6
    K = np.zeros((ndof, ndof))
    F = np.zeros(ndof)
    for element, Ke in zip(elements, built["stiffness"]):
        dof = np.asarray(element.assemble_global_dof_indices(), dtype=int)
        K[np.ix_(dof, dof)] += Ke
        F[dof] += np.asarray(element.element_force_vector().F_e, dtype=np.float64).ravel()
    free = np.arange(6, ndof)
    U = np.zeros(ndof)
    U[free] = np.linalg.solve(K[np.ix_(free, free)], F[free])
    return U


def _runner_settings(built, simulation_settings):
    return {
        "elements": np.asarray(built["elements"], dtype=object),
        "mesh_dictionary": built["mesh"],
        "grid_dictionary": built["grid"],
        "element_dictionary": built["element_dictionary"],
        "material_dictionary": built["material"],
        "section_dictionary": built["section"],
        "point_load_array": built["point_load"],
        "distributed_load_array": built["distributed"],
        "job_results_dir": built["job"],
        "simulation_settings": simulation_settings,
        "element_stiffness_matrices": np.stack(built["stiffness"]),
        "element_mass_matrices": np.stack(built["mass"]),
        "element_objects": np.asarray(built["element_objects"], dtype=object),
        "force_objects": np.asarray(built["force_objects"], dtype=object),
        "prescribed_displacement_dict": None,
    }


def test_deflection_profile_of_a_tip_loaded_cantilever():
    """Hermite beam, tip force in +z: the sampled profile is the cubic cantilever."""
    L, E, I, P = 2.0, 210.0e9, 1.0e-6, 1000.0
    built = _cantilever(4, L, E, I, 7850.0, (0.0, 0.0, P))
    try:
        U = _static_displacement(built)
        profile = deflection_profile(built["coordinates"], U, ELEMENT, n_per_element=5)
        analytic = P * profile.s**2 * (3.0 * L - profile.s) / (6.0 * E * I)
        np.testing.assert_allclose(profile.displacement[:, 2], analytic, rtol=1e-8, atol=1e-12)
        np.testing.assert_allclose(profile.displacement[:, 0], 0.0, atol=1e-12)
        np.testing.assert_allclose(profile.displacement[:, 1], 0.0, atol=1e-12)
        assert profile.s[0] == 0.0
        assert profile.s[-1] == pytest.approx(L)
        assert profile.position.shape == (profile.s.size, 3)
    finally:
        shutil.rmtree(built["temp"], ignore_errors=True)


def test_extreme_fibre_stress_and_strain_along_the_outer_profile():
    """At z = +c the fibre is in compression under an upward tip load; stress is E times strain."""
    L, E, I, P, c = 2.0, 210.0e9, 1.0e-6, 1000.0, 0.05
    built = _cantilever(4, L, E, I, 7850.0, (0.0, 0.0, P))
    try:
        U = _static_displacement(built)
        outer = np.array([[0.0, c], [0.0, -c], [0.02, 0.0]])
        fibre = extreme_fibre_profile(built["elements"], U, outer, xi=0.0)
        moment = P * (L - fibre.s)
        sigma_top = -moment * c / I
        np.testing.assert_allclose(fibre.stress[:, 0], sigma_top, rtol=1e-6, atol=1e-6)
        np.testing.assert_allclose(fibre.stress[:, 1], -sigma_top, rtol=1e-6, atol=1e-6)
        np.testing.assert_allclose(fibre.stress[:, 2], 0.0, atol=1e-4)
        np.testing.assert_allclose(fibre.stress, E * fibre.strain, rtol=1e-12, atol=0.0)
        np.testing.assert_allclose(fibre.extreme_stress[:, 0], sigma_top, rtol=1e-6, atol=1e-6)
        np.testing.assert_allclose(fibre.extreme_stress[:, 1], -sigma_top, rtol=1e-6, atol=1e-6)
        np.testing.assert_allclose(fibre.extreme_strain[:, 0], sigma_top / E, rtol=1e-6, atol=0.0)
        assert fibre.y.tolist() == [0.0, 0.0, 0.02]
        assert fibre.z.tolist() == [c, -c, 0.0]
    finally:
        shutil.rmtree(built["temp"], ignore_errors=True)


def test_extreme_fibre_profile_rejects_a_bad_contour():
    with pytest.raises(ValueError, match="outer_profile"):
        extreme_fibre_profile([], np.zeros(12), np.zeros(3))


def test_buckling_result_reads_the_solved_eigenpairs_and_does_not_solve_again():
    factors = np.array([3.5, 12.0])
    shapes = np.arange(8.0).reshape(4, 2)

    class _Solved:
        buckling_load_factors = factors
        buckling_mode_shapes = shapes

        def solve_buckling_linearized_problem(self, *args, **kwargs):
            raise AssertionError("buckling was solved again")

    result = buckling_result(_Solved())
    np.testing.assert_array_equal(result.load_factors, factors)
    np.testing.assert_array_equal(result.mode_shapes, shapes)
    again = buckling_result((factors, shapes))
    np.testing.assert_array_equal(again.load_factors, factors)
    with pytest.raises(RuntimeError, match="buckling result"):
        buckling_result(object())


def test_blade_vibration_modes_read_the_solved_modes_and_do_not_build_a_rotor():
    frequencies = np.array([1.25, 7.5])
    shapes = np.ones((6, 2))

    class _Solved:
        vibration_frequencies_hz = frequencies
        vibration_mode_shapes = shapes

        def solve_vibration_eigenpairs(self, *args, **kwargs):
            raise AssertionError("blade modes were recomputed")

        def _run_rotating_frame_vibration(self, *args, **kwargs):
            raise AssertionError("a rotor mode was solved")

    modes = blade_vibration_modes(_Solved())
    np.testing.assert_array_equal(modes.frequencies_hz, frequencies)
    np.testing.assert_array_equal(modes.mode_shapes, shapes)
    packed = blade_vibration_modes((frequencies, shapes))
    np.testing.assert_array_equal(packed.frequencies_hz, frequencies)
    with pytest.raises(RuntimeError, match="blade vibration modes"):
        blade_vibration_modes(object())


def test_linear_buckling_run_returns_the_cantilever_euler_load():
    """Eight elements, left as they are. The returned factor is the solver's, near the Euler load."""
    from onedfem.simulation_runner.buckling.buckling_simulation import LinearBucklingSimulationRunner

    L, E, I, P_ref = 2.5, 210.0e9, 1.0e-6, 1.0
    built = _cantilever(8, L, E, I, 7850.0, (-P_ref, 0.0, 0.0))
    try:
        settings = _runner_settings(built, {"buckling": {"num_modes": 4, "buckling_prestress": "linear_static",
                                                         "buckling_load_factor": 1.0}})
        runner = LinearBucklingSimulationRunner(settings=settings, job_name="limit_state_buckling")
        result = runner.run()
        np.testing.assert_allclose(result.load_factors, np.asarray(runner.buckling_load_factors).ravel())
        np.testing.assert_allclose(result.load_factors, buckling_result(runner).load_factors)
        euler = (np.pi ** 2) * E * I / (4.0 * L ** 2)
        critical = np.min(np.abs(result.load_factors[:4] * P_ref))
        np.testing.assert_allclose(critical, euler, rtol=0.15)
    finally:
        shutil.rmtree(built["temp"], ignore_errors=True)


def test_eigen_run_returns_blade_modes_near_the_cantilever_frequency():
    """Eight elements, left as they are. The returned modes are the beam's, not a rotor's."""
    from onedfem.simulation_runner.eigen.eigen_simulation import EigenSimulationRunner

    L, E, I, rho, A = 2.0, 70.0e9, 1.0e-5, 2700.0, 0.01
    built = _cantilever(8, L, E, I, rho, (0.0, 0.0, 0.0))
    try:
        settings = _runner_settings(built, {"eigen": {"num_modes": 2}})
        runner = EigenSimulationRunner(settings=settings, job_name="limit_state_modes")
        modes = runner.run()
        assert "rotation" not in (runner.simulation_settings or {})
        np.testing.assert_allclose(modes.frequencies_hz, blade_vibration_modes(runner).frequencies_hz)
        beta = 1.875104068711961
        analytic = (beta ** 2 / L ** 2) * np.sqrt(E * I / (rho * A)) / (2.0 * np.pi)
        np.testing.assert_allclose(np.min(modes.frequencies_hz[:2]), analytic, rtol=0.02)
        assert modes.mode_shapes.shape[1] == 2
    finally:
        shutil.rmtree(built["temp"], ignore_errors=True)

"""Corotational elements on coupled sections: large deflection, large twist and stability.

References:

- Bisshopp and Drucker (1945): a cantilever with tip load ``P L^2 / EI = 1`` deflects by
  ``0.30172 L`` vertically and ``0.05643 L`` horizontally.
- A tip moment ``M = pi EI / L`` bends a cantilever into a semicircle.
- Uniform twist rate ``k`` with free warping and free axial end (Wagner stiffening):
  ``T = GJ k + 1/2 E (I_rr - I_p^2 / A) k^3``, with ``I_rr = int r^4 dA``.
- Restrained warping at the root: the nonlinear Vlasov equation
  ``GJ phi' - E Gamma phi''' + 1/2 E I_n phi'^3 = T``, solved here by collocation.
- Lateral-torsional buckling under uniform moment, including the prebuckling in-plane
  deflection (Trahair 1993, eq. 6.24)::

      M_cr = M_cr0 / sqrt((1 - EI_z / EI_y) (1 - GJ (1 + pi^2 E Gamma / (GJ L^2)) / EI_y))
"""

from __future__ import annotations

import numpy as np
import pytest
from scipy.integrate import solve_bvp

from tests.test_second_order import (
    E,
    G,
    Beam,
    channel,
    end_moments,
    i_beam,
    i_beam_properties,
)

CR_EB, CR_TS, CR_BS = (
    "CorotationalEulerBernoulliBeamElement3D",
    "CorotationalTimoshenkoBeamElement3D",
    "CorotationalBenscoterBeamElement3D",
)


def assemble(beam, U):
    Fi = np.zeros(beam.n)
    Kt = np.zeros((beam.n, beam.n))
    for el, g in zip(beam.elements, beam.maps):
        Fi[g] += el.internal_force_vector(U[g])
        Kt[np.ix_(g, g)] += el.tangent_stiffness_matrix(U[g])
    return Fi, Kt


def newton_path(beam, F_of, loads, fixed, tol=1e-8, max_it=30):
    """Load-controlled Newton; yields ``(load, U, K_t)`` at each converged step."""
    free = np.setdiff1d(np.arange(beam.n), fixed)
    U = np.zeros(beam.n)
    for lam in loads:
        F = F_of(lam)
        for _ in range(max_it):
            Fi, Kt = assemble(beam, U)
            r = (F - Fi)[free]
            if np.linalg.norm(r) <= tol * np.linalg.norm(F):
                break
            U[free] += np.linalg.solve(Kt[np.ix_(free, free)], r)
        else:
            raise AssertionError(f"Newton did not converge at load {lam}")
        yield lam, U.copy(), Kt[np.ix_(free, free)]


def _tip_load(beam, dof, value):
    F = np.zeros(beam.n)
    F[beam.dof(beam.n_el, dof)] = value
    return F


# --------------------------------------------------------------------- #
# Element properties
# --------------------------------------------------------------------- #
@pytest.mark.parametrize("etype", [CR_EB, CR_TS, CR_BS])
def test_rigid_motion_is_stress_free_and_tangent_is_consistent(etype):
    from pre_processing.element_library.nonlinear.large_rotations.corotational.corotational_section import (
        rotation_exp,
    )

    beam = Beam(etype, 0.5, 1, sa=channel())
    el = beam.elements[0]
    n = el._n_dof
    a = np.array([0.4, -0.7, 1.1])
    U = np.zeros(n)
    U[3:6] = U[9:12] = a
    U[6:9] = (rotation_exp(a) - np.eye(3)) @ np.array([0.5, 0.0, 0.0])
    np.testing.assert_allclose(el.local_deformation(U), 0.0, atol=1e-12)
    assert np.abs(el.internal_force_vector(U)).max() < 1e-6 * np.abs(el._K_core).max() * 1e-6
    # small state: the linear element; larger state: tangent = derivative of the force
    np.testing.assert_allclose(el.tangent_stiffness_matrix(np.zeros(n)), el._K_core, atol=1e-9 * np.abs(el._K_core).max())
    rng = np.random.default_rng(0)
    U = rng.normal(size=n) * 2e-2
    Kt = el.tangent_stiffness_matrix(U)
    h = 1e-7
    Kfd = np.column_stack([(el.internal_force_vector(U + h * e) - el.internal_force_vector(U - h * e)) / (2 * h)
                           for e in np.eye(n)])
    np.testing.assert_allclose(Kt, Kfd, atol=1e-6 * np.abs(Kfd).max())


# --------------------------------------------------------------------- #
# Large deflection
# --------------------------------------------------------------------- #
@pytest.mark.parametrize("etype", [CR_EB, CR_TS])
def test_elastica_tip_load(etype):
    L, n = 1.0, 16
    beam = Beam(etype, L, n, rect=(0.02, 0.05))
    EI = E * beam.rect["I_y"]
    *_, (_, U, _) = newton_path(beam, lambda s: _tip_load(beam, 2, -s * EI / L**2), np.linspace(0.25, 1, 4),
                                beam.clamp_root())
    # Timoshenko adds the shear deflection of a beam with L/h = 20
    tol = 1e-3 if etype == CR_EB else 3e-3
    assert -U[beam.dof(n, 2)] / L == pytest.approx(0.30172, rel=tol)
    assert -U[beam.dof(n, 0)] / L == pytest.approx(0.05643, rel=3 * tol)


@pytest.mark.parametrize("etype", [CR_EB, CR_TS])
def test_tip_moment_rolls_a_semicircle(etype):
    L, n = 1.0, 32
    beam = Beam(etype, L, n, rect=(0.02, 0.05))
    EI = E * beam.rect["I_y"]
    *_, (_, U, _) = newton_path(beam, lambda s: _tip_load(beam, 4, -s * np.pi * EI / L), np.linspace(0.125, 1, 8),
                                beam.clamp_root())
    assert U[beam.dof(n, 0)] == pytest.approx(-L, rel=2e-3)          # tip back above the root
    assert U[beam.dof(n, 2)] == pytest.approx(2 * L / np.pi, rel=2e-3)
    assert abs(U[beam.dof(n, 4)]) == pytest.approx(np.pi, rel=1e-3)


# --------------------------------------------------------------------- #
# Large twist
# --------------------------------------------------------------------- #
def _i_beam_twist_constants():
    from numpy.polynomial.legendre import leggauss

    h, b, tf, tw = 0.3, 0.15, 6e-3, 4e-3
    p = i_beam_properties()
    x, w = leggauss(10)
    Irr = sum(tf * 0.5 * b * np.sum(w * ((0.5 * b * x) ** 2 + z0**2) ** 2) for z0 in (h / 2, -h / 2))
    Irr += tw * 0.5 * h * np.sum(w * (0.5 * h * x) ** 4)
    I_n = Irr - (p["I_y"] + p["I_z"]) ** 2 / p["A"]
    return p, I_n


def test_wagner_torsional_stiffening_free_warping():
    p, I_n = _i_beam_twist_constants()
    an = i_beam()
    GJ = 1.0 / an.flexibility[5, 5]
    L, n = 2.0, 8
    beam = Beam(CR_TS, L, n, sa=an)
    for T, U, _ in newton_path(beam, lambda T: _tip_load(beam, 3, T), np.linspace(250, 2000, 8), beam.clamp_root()):
        k = U[beam.dof(n, 3)] / L
        assert T == pytest.approx(GJ * k + 0.5 * E * I_n * k**3, rel=1e-6)
    assert k * L < 0.5 * T / GJ * L  # far from linear: the linear twist would be 1.9 rad
    # the section shortens by the mean of r^2 k^2 / 2
    assert U[beam.dof(n, 0)] == pytest.approx(-0.5 * (p["I_y"] + p["I_z"]) / p["A"] * k**2 * L, rel=1e-6)


def test_benscoter_large_twist_restrained_root():
    p, I_n = _i_beam_twist_constants()
    an = i_beam()
    GJ, EG = 1.0 / an.flexibility[5, 5], E * p["Gamma"]
    L, n, T = 2.0, 16, 6000.0

    def rhs(x, y):
        return np.vstack([y[1], y[2], (GJ * y[1] + 0.5 * E * I_n * y[1] ** 3 - T) / EG])

    x = np.linspace(0, L, 200)
    y0 = np.zeros((3, x.size))
    y0[1] = T / GJ * (1 - np.exp(-5 * x))
    ref = solve_bvp(rhs, lambda a, b: np.array([a[0], a[1], b[2]]), x, y0, tol=1e-10, max_nodes=100000)
    beam = Beam(CR_BS, L, n, sa=an)
    *_, (_, U, _) = newton_path(beam, lambda s: _tip_load(beam, 3, s), np.linspace(500, T, 12), list(range(7)))
    # Benscoter adds the warping shear flexibility that Vlasov leaves out
    assert U[beam.dof(n, 3)] == pytest.approx(ref.sol(L)[0], rel=5e-3)
    linear_twist = np.linalg.solve(beam.K[7:, 7:], _tip_load(beam, 3, T)[7:])[beam.dof(n, 3) - 7]
    assert U[beam.dof(n, 3)] < 0.8 * linear_twist


# --------------------------------------------------------------------- #
# Stability along the path
# --------------------------------------------------------------------- #
@pytest.mark.parametrize("etype", [CR_EB, CR_BS])
def test_lateral_torsional_stability_loss_with_prebuckling_deflection(etype):
    """The tangent loses positive definiteness where the closed form with prebuckling deflection says."""
    p = i_beam_properties()
    an = i_beam()
    L, n = 6.0, 16
    Gam = p["Gamma"] if etype == CR_BS else 0.0
    GJ = 1.0 / an.flexibility[5, 5]
    EIy, EIz = E * p["I_y"], E * p["I_z"]
    M0 = np.pi / L * np.sqrt(EIz * (GJ + np.pi**2 * E * Gam / L**2))
    ref = M0 / np.sqrt((1 - EIz / EIy) * (1 - (GJ + np.pi**2 * E * Gam / L**2) / EIy))
    beam = Beam(etype, L, n, sa=an)
    imperfection = np.zeros(beam.n)
    imperfection[beam.dof(n // 2, 1)] = 1e-3 * M0 / L
    loads = ref * np.array([0.5, 0.8, 0.9, 0.95])
    lam_min = [(m, np.linalg.eigvalsh(0.5 * (K + K.T)).min())
               for m, _, K in newton_path(beam, lambda m: end_moments(beam, m) + imperfection, loads, beam.fork())]
    (m1, e1), (m2, e2) = lam_min[-2:]
    m_cr = m2 - e2 * (m2 - m1) / (e2 - e1)
    assert m_cr == pytest.approx(ref, rel=1e-2)


def test_large_twist_example_through_the_nonlinear_runner(tmp_path):
    import importlib.util
    import os

    path = os.path.join(os.path.dirname(__file__), "..", "examples", "thin_walled", "large_twist", "run.py")
    spec = importlib.util.spec_from_file_location("large_twist_run", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    s = mod.main(n_el=8, plot=False, work_dir=str(tmp_path), torques=(6000.0,))
    row = s["rows"][0]
    assert row["corotational"] == pytest.approx(row["vlasov_nonlinear"], rel=1e-2)
    assert row["corotational"] < 0.8 * row["linear"]

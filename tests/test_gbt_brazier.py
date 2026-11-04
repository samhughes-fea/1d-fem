"""Brazier flattening of corotational GBT beams.

References:

- ``brazier_energy_derivatives`` is the gradient and Hessian of its energy (central
  differences).
- An isotropic tube (a 24-sided polygon) under a tip moment ovalises as Brazier's long
  tube: ``zeta = (1 - nu^2) kappa^2 R^4 / t^2`` (Brazier 1927), with the bending stiffness
  falling by ``1.5 zeta`` against the same beam without the term.
"""

from __future__ import annotations

import csv
import logging
import os
import sys

import numpy as np
import pytest

from pre_processing.element_library.nonlinear.large_rotations.corotational.corotational_section import (
    brazier_energy_derivatives,
)

HERE = os.path.dirname(os.path.abspath(__file__))


def test_brazier_energy_derivatives_match_central_differences():
    rng = np.random.default_rng(0)
    m, n = 3, 7
    T = rng.normal(size=(m, 2, n))
    G = rng.normal(size=(m, 2, 4))
    Q = np.einsum("kci,ldi->klcd", G, G)

    def energy(z):
        a, k, E = z[:m], z[m:m + 2], z[m + 2:]
        return a @ np.einsum("kce,c,e->k", T, k, E) + 0.5 * np.einsum("klcd,k,l,c,d->", Q, a, a, k, k)

    z = rng.normal(size=m + 2 + n)
    g, H = brazier_energy_derivatives(T, Q, z)
    h = 1e-6
    eye = np.eye(z.size)
    gf = np.array([(energy(z + h * e) - energy(z - h * e)) / (2 * h) for e in eye])
    Hf = np.array([(brazier_energy_derivatives(T, Q, z + h * e)[0] - brazier_energy_derivatives(T, Q, z - h * e)[0])
                   / (2 * h) for e in eye])
    np.testing.assert_allclose(g, gf, atol=1e-7 * np.abs(gf).max())
    np.testing.assert_allclose(H, Hf, atol=1e-7 * np.abs(Hf).max())
    np.testing.assert_array_equal(H, H.T)


R, T_WALL, E, NU, N_WALLS = 0.1, 2e-3, 70e9, 0.3, 24
L, N_EL, N_MODES = 4.0, 16, 4
_NL = ("[Simulation]\n[Type]\nstatic_nonlinear\n\n[Newton]\ntolerance = 1e-6\nrelative_tolerance = 1e-10\n"
       "relative_reference = first_residual\nmax_iterations = 40\ntolerance_delta_u = 1e-14\n\n[Solver]\n"
       "type = direct\ntolerance = 1e-12\n\n[Nonlinear]\nnum_increments = 6\nline_search = false\n")


def _tube():
    sys.path.insert(0, os.path.join(HERE, "..", "examples", "thin_walled"))
    import common

    G = E / (2 * (1 + NU))
    pts = {i + 1: (R * np.cos(2 * np.pi * i / N_WALLS), R * np.sin(2 * np.pi * i / N_WALLS)) for i in range(N_WALLS)}
    walls = [(i + 1, i + 1, (i + 1) % N_WALLS + 1, 1) for i in range(N_WALLS)]
    return common, common.SectionSpec(pts, walls, {1: common.Material(E, E, G, NU, G, G, 2700.0)},
                                      {1: [(1, 0.0, T_WALL)]}, n_sub=2)


def _run(tmp_path, moment, name):
    from workflow_orchestrator.run_job import process_job

    common, sec = _tube()
    spec = common.CantileverSpec(L, N_EL, "CorotationalGBTBeamElement3D", sec, tip_load=[0, 0, 0, 0, moment, 0],
                                 modes=N_MODES, warping=1)
    jd = str(tmp_path / "jobs" / name)
    common.write_job(jd, spec)
    common.write_supports(jd, spec)
    with open(os.path.join(jd, "simulation_settings.txt"), "w", encoding="utf-8") as f:
        f.write(_NL)
    res = common.results_dir(str(tmp_path / "res" / name))
    logging.getLogger().setLevel(logging.WARNING)
    process_job(jd, res, {}, {}, force_serial=True, max_processes_per_job=1)
    with open(os.path.join(res, "primary_results", "global", "U_global.csv"), newline="", encoding="utf-8") as f:
        U = np.array([float(r[1]) for r in list(csv.reader(f))[1:]])
    return U.reshape(N_EL + 1, -1)


@pytest.mark.slow
def test_tube_ovalises_as_brazier(tmp_path, monkeypatch):
    from pre_processing.section_integration.section_analysis import SectionAnalysis

    M = 1.0e4                                   # zeta about 1.2 %
    U = _run(tmp_path, M, "brazier")
    k = N_EL // 2
    dx = L / N_EL
    kappa = (U[k + 1, 4] - U[k - 1, 4]) / (2 * dx)
    n_w = U.shape[1] - 6 - 2 * N_MODES
    a = U[k, 6 + n_w::2][:N_MODES]
    # the ovalisation pair (unit peak in-plane displacement) carries it; the next modes do not
    zeta = np.hypot(a[0], a[1]) / R
    theory = (1 - NU ** 2) * kappa ** 2 * R ** 4 / T_WALL ** 2
    assert zeta == pytest.approx(theory, rel=0.02)
    assert np.abs(a[2:]).max() < 1e-3 * np.hypot(a[0], a[1])

    monkeypatch.setattr(SectionAnalysis, "gbt_brazier_terms", lambda self, *args: None)
    U0 = _run(tmp_path, M, "rigid_section")
    kappa0 = (U0[k + 1, 4] - U0[k - 1, 4]) / (2 * dx)
    assert np.abs(U0[:, 6:]).max() < 1e-12
    assert kappa0 / kappa == pytest.approx(1 - 1.5 * theory, abs=0.1 * 1.5 * theory)

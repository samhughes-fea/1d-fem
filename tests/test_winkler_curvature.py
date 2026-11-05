"""Winkler curvature correction of curved elements' sections (``SectionAnalysis.with_curvature``).

On a curved axis the fibres inside the bend are shorter than those outside, so the
section's axial and bending block is ``int E a a^T / g dA`` with ``g = 1 - z k_y - y k_z``.
Checked against the integrals themselves and, on a thick quarter-circle bar under a tip
moment, against the Winkler-Bach curved beam (rotation ``M L / (E A e R)``).
"""

import csv
import os
import shutil

import numpy as np
import pytest
from scipy.integrate import quad

from onedfem.pre_processing.section_integration.section_analysis import SectionAnalysis

E, G = 70e9, 26e9


def _rect(b, h):
    """b along z, h along y."""
    A, I_y, I_z = b * h, h * b**3 / 12, b * h**3 / 12
    return SectionAnalysis.from_isotropic(E=E, G=G, A=A, I_x=I_y + I_z, I_y=I_y, I_z=I_z, J_t=0.2 * b**3 * h,
                                          kappa=5 / 6, rho=2700.0)


def test_section_block_is_the_curved_fibre_integral():
    b, h, k = 0.1, 0.5, 1.0
    sa = _rect(b, h)
    D = sa.with_curvature(0.0, k).stiffness
    g = lambda y: 1.0 - k * y                                   # noqa: E731
    ref = np.array([[quad(lambda y: E * b * f(y) * q(y) / g(y), -h / 2, h / 2)[0]
                     for q in (lambda y: 1.0, lambda y: -y)] for f in (lambda y: 1.0, lambda y: -y)])
    got = D[np.ix_([0, 2], [0, 2])]
    np.testing.assert_allclose(got, ref, rtol=1e-6)
    # first order: extension-bending coupling -k EI (here a_3 = -y, so D[eps, kappa_z] = -k EI_z)
    assert D[0, 2] == pytest.approx(-k * E * b * h**3 / 12, rel=0.05)
    # untouched: shear, torsion, and the other bending plane's coupling to this one
    S = sa.stiffness
    np.testing.assert_allclose(D[3:, 3:], S[3:, 3:])
    assert sa.with_curvature(0.0, 0.0) is sa


def test_centre_of_curvature_inside_the_section_is_refused():
    with pytest.raises(ValueError, match="centre of curvature"):
        _rect(0.1, 0.5).with_curvature(0.0, 10.0)


def _write_job(job, n, b, h, R, M):
    t = np.linspace(0.0, np.pi / 2, n + 1)
    xyz = np.column_stack([R * np.sin(t), R * (1 - np.cos(t)), 0 * t])
    job.mkdir(parents=True)
    (job / "grid.txt").write_text("[Grid]\n[node_id] [x] [y] [z]\n"
                                  + "".join(f"{i} {x:.15g} {y:.15g} {z:.15g}\n" for i, (x, y, z) in enumerate(xyz)))
    (job / "element.txt").write_text(
        "[Element]\n[element_id] [node1] [node2] [element_type] [axial_order] [bending_y_order] [bending_z_order] "
        "[shear_y_order] [shear_z_order] [torsion_order] [load_order]\n"
        + "".join(f"{e} {e} {e + 1} LinearTimoshenkoBeamElement3D 2 2 2 2 2 2 2\n" for e in range(n)))
    (job / "material.txt").write_text("[Material]\n[element_id] [E] [G] [nu] [rho]\n"
                                      + "".join(f"{e} {E} {G} 0.3 2700\n" for e in range(n)))
    A, I_y, I_z = b * h, h * b**3 / 12, b * h**3 / 12
    (job / "section.txt").write_text("[Section]\n[element_id] [A] [I_x] [I_y] [I_z] [J_t] [kappa] [alpha]\n"
                                     + "".join(f"{e} {A} {I_y + I_z} {I_y} {I_z} {0.2 * b**3 * h} 0.8333 0.0\n"
                                               for e in range(n)))
    (job / "precurvature.txt").write_text("[Precurvature]\n[element_id] [k_x0] [k_y0] [k_z0]\n"
                                          + "".join(f"{e} 0 0 {1 / R}\n" for e in range(n)))
    (job / "prescribed_displacement.txt").write_text(
        "[Prescribed Displacement]\n[id] [node_id] [dof] [value] [type] [comment]\n"
        + "".join(f"{i} 0 {d} 0.0 displacement\n" for i, d in enumerate(("UX", "UY", "UZ", "RX", "RY", "RZ"))))
    x, y, z = xyz[-1]
    (job / "point_load.txt").write_text("[Point load]\n[x] [y] [z] [F_x] [F_y] [F_z] [M_x] [M_y] [M_z]\n"
                                        f"{x:.15g} {y:.15g} {z:.15g} 0 0 0 0 0 {M}\n")
    (job / "simulation_settings.txt").write_text(
        "[Simulation]\n[Type]\nstatic\n\n[Solver]\ntype = direct\ntolerance = 1e-12\nmax_iterations = 500\n")
    return A, I_z


def test_thick_curved_bar_follows_winkler_bach(tmp_path):
    from onedfem.workflow_orchestrator.run_job import process_job, setup_job_results_directory

    n, b, h, R, M = 16, 0.1, 0.5, 1.0, 1e5
    A, I_z = _write_job(tmp_path / "job", n, b, h, R, M)
    res = setup_job_results_directory("test_thick_curved_bar_follows_winkler_bach")
    try:
        process_job(str(tmp_path / "job"), res, {}, {}, force_serial=True, max_processes_per_job=1)
        with open(os.path.join(res, "primary_results", "global", "U_global.csv"), newline="") as fh:
            U = np.array([float(r[1]) for r in list(csv.reader(fh))[1:]]).reshape(n + 1, -1)
    finally:
        shutil.rmtree(res, ignore_errors=True)
    L = np.pi * R / 2
    e = R - h / np.log((R + h / 2) / (R - h / 2))              # centroid to neutral axis
    winkler = M * L / (E * A * e * R)
    straight = M * L / (E * I_z)
    assert U[n, 5] == pytest.approx(winkler, rel=2e-3)
    assert abs(winkler / straight - 1.0) > 0.01                 # the correction is what is measured

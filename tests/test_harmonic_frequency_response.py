"""§4 harmonic frequency sweep — kernels + end-to-end ``process_job``."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pytest
from scipy.sparse import csr_matrix, eye

PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from processing.harmonic.damping_table import interpolate_zeta_hz, load_zeta_vs_frequency_hz
from processing.harmonic.frequency_response import (
    _dynamic_matrix,
    dynamic_matrix_csc_signature,
    harmonic_damping_matrix,
    mass_proportional_damping,
    solve_one_frequency,
    solve_one_frequency_partitioned,
    sweep_displacements,
)
from processing.harmonic.modal_superposition import (
    harmonic_displacement_modal_superposition,
    harmonic_truncation_metrics_vs_direct,
    undamped_natural_modes,
)
from processing.common.harmonic_complex_post import build_harmonic_complex_snapshots
from simulation_runner.harmonic.harmonic_simulation import effective_harmonic_config
from workflow_orchestrator.run_job import process_job, setup_job_results_directory


def test_effective_harmonic_config_defaults() -> None:
    cfg = effective_harmonic_config({"type": "harmonic", "harmonic": {}})
    assert "frequency_min_hz" in cfg
    assert cfg["num_frequency_points"] >= 2
    assert cfg["mp_damping_reference"] == "geometric_mean"
    assert cfg["use_modal_superposition"] is False
    assert cfg["harmonic_linear_solver"] == "spsolve"


def test_solve_one_frequency_sdof_matches_dense() -> None:
    n = 4
    omega = 10.0
    K = eye(n, format="csr", dtype=np.float64)
    M = eye(n, format="csr", dtype=np.float64) * 0.25
    zeta = 0.05
    omega_ref = 12.0
    C = mass_proportional_damping(M, zeta, omega_ref)
    F = np.zeros(n, dtype=np.float64)
    F[1] = 1.0
    u_sparse = solve_one_frequency(K, M, C, F, omega)
    A_dense = (
        -omega**2 * M.toarray() + 1j * omega * C.toarray() + K.toarray()
    ).astype(np.complex128)
    u_dense = np.linalg.solve(A_dense, F.astype(np.complex128))
    np.testing.assert_allclose(u_sparse, u_dense, rtol=1e-9, atol=1e-9)


def test_harmonic_sdof_matches_analytical_reference() -> None:
    n = 1
    omega = 8.0
    K = eye(n, format="csr", dtype=np.float64) * 100.0
    M = eye(n, format="csr", dtype=np.float64) * 1.0
    C = eye(n, format="csr", dtype=np.float64) * 0.0
    F = np.array([1.0], dtype=np.float64)
    u = solve_one_frequency(K, M, C, F, omega)
    u_ref = np.array([1.0 / (100.0 - omega**2)], dtype=np.complex128)
    np.testing.assert_allclose(u, u_ref, rtol=1e-9, atol=1e-10)


def test_sweep_two_points_column_layout() -> None:
    n = 3
    K = eye(n, format="csr") * 100.0
    M = eye(n, format="csr") * 1.0
    C = eye(n, format="csr") * 0.1
    F = np.array([0.0, 1.0, 0.0])
    freqs = np.array([2.0, 4.0])
    U = sweep_displacements(K, M, C, F, freqs)
    assert U.shape == (n, 2)


def test_harmonic_damping_matrix_rayleigh_matches_sum() -> None:
    n = 5
    K = eye(n, format="csr") * 3.0
    M = eye(n, format="csr") * 2.0
    zeta = 0.01
    omega_ref = 5.0
    ra, rb = 0.02, 1e-4
    C = harmonic_damping_matrix(M, K, zeta, omega_ref, ra, rb)
    C_exp = mass_proportional_damping(M, zeta, omega_ref) + M * ra + K * rb
    np.testing.assert_allclose(C.toarray(), C_exp.toarray(), rtol=1e-12)


def test_sweep_parallel_matches_serial() -> None:
    n = 4
    K = eye(n, format="csr") * 50.0
    M = eye(n, format="csr") * 0.5
    C = eye(n, format="csr") * 0.05
    F = np.array([1.0, 0.0, 0.5, 0.0])
    freqs = np.linspace(1.0, 10.0, 9)
    Us = sweep_displacements(K, M, C, F, freqs, parallel=False)
    Up = sweep_displacements(K, M, C, F, freqs, parallel=True, max_workers=4)
    np.testing.assert_allclose(Us, Up, rtol=1e-12, atol=1e-12)


def test_partitioned_prescribed_dof_residual() -> None:
    """Prescribed-motion partition satisfies reduced equilibrium."""
    n = 3
    K = eye(n, format="csr") * 10.0
    M = eye(n, format="csr") * 1.0
    C = eye(n, format="csr") * 0.2
    omega = 4.0
    rng = np.random.default_rng(42)
    F = rng.standard_normal(n).astype(np.float64)
    u_p = np.zeros(n, dtype=np.complex128)
    u_p[0] = 0.25 + 0.1j
    unk = np.array([1, 2], dtype=np.int64)
    u_sol = solve_one_frequency_partitioned(K, M, C, F, omega, u_p, unk)
    np.testing.assert_allclose(u_sol[0], u_p[0])
    om2 = omega * omega
    Kc = K.astype(np.complex128)
    Mc = M.astype(np.complex128)
    Cc = C.astype(np.complex128)
    A = Kc + (-om2) * Mc + (1j * omega) * Cc
    res = A @ u_sol - np.asarray(F, dtype=np.complex128)
    assert np.linalg.norm(res[unk]) < 1e-9 * max(1.0, np.linalg.norm(F[unk]))


def test_solve_one_frequency_splu_matches_spsolve() -> None:
    n = 3
    K = eye(n, format="csr") * 4.0
    M = eye(n, format="csr") * 0.5
    C = eye(n, format="csr") * 0.1
    F = np.array([0.0, 1.0, 0.0])
    omega = 3.0
    u0 = solve_one_frequency(K, M, C, F, omega, linear_solver="spsolve")
    u1 = solve_one_frequency(K, M, C, F, omega, linear_solver="splu")
    np.testing.assert_allclose(u0, u1, rtol=1e-10, atol=1e-12)


def test_sweep_splu_csc_buffer_matches_without_buffer(monkeypatch: pytest.MonkeyPatch) -> None:
    n = 3
    K = eye(n, format="csr") * 50.0
    M = eye(n, format="csr") * 0.5
    C = eye(n, format="csr") * 0.05
    F = np.array([1.0, 0.0, 0.0])
    freqs = np.linspace(1.0, 8.0, 7)
    monkeypatch.setenv("FEM_HARMONIC_SPLU_CSC_BUFFER", "1")
    U_buf = sweep_displacements(K, M, C, F, freqs, linear_solver="splu", parallel=False)
    monkeypatch.setenv("FEM_HARMONIC_SPLU_CSC_BUFFER", "0")
    U_plain = sweep_displacements(K, M, C, F, freqs, linear_solver="splu", parallel=False)
    np.testing.assert_allclose(U_buf, U_plain, rtol=1e-12, atol=1e-12)


def test_undamped_natural_modes_sparse_branch_matches_dense() -> None:
    n = 24
    K = eye(n, format="csr") * 60.0
    M = eye(n, format="csr") * 0.35
    om_d, _Phi_d = undamped_natural_modes(K, M, 6, dense_threshold=10_000)
    om_s, _Phi_s = undamped_natural_modes(K, M, 6, dense_threshold=0)
    np.testing.assert_allclose(om_d, om_s, rtol=1e-6, atol=1e-8)


def test_sweep_zeta_per_frequency_table_matches_interpolated_scalars() -> None:
    n = 2
    K = eye(n, format="csr") * 100.0
    M = eye(n, format="csr") * 1.0
    F = np.array([1.0, 0.0])
    freqs = np.array([10.0, 20.0, 30.0])
    z0, z1, z2 = 0.01, 0.02, 0.03
    C0 = harmonic_damping_matrix(M, K, z0, 2.0 * np.pi * 10.0, 0.0, 0.0)
    C1 = harmonic_damping_matrix(M, K, z1, 2.0 * np.pi * 20.0, 0.0, 0.0)
    C2 = harmonic_damping_matrix(M, K, z2, 2.0 * np.pi * 30.0, 0.0, 0.0)
    u0 = solve_one_frequency(K, M, C0, F, 2.0 * np.pi * 10.0)
    u1 = solve_one_frequency(K, M, C1, F, 2.0 * np.pi * 20.0)
    u2 = solve_one_frequency(K, M, C2, F, 2.0 * np.pi * 30.0)
    C_dummy = C0
    U = sweep_displacements(
        K,
        M,
        C_dummy,
        F,
        freqs,
        mp_damping_reference="current_frequency",
        zeta_per_frequency=np.array([z0, z1, z2]),
    )
    np.testing.assert_allclose(U[:, 0], u0, rtol=1e-9, atol=1e-9)
    np.testing.assert_allclose(U[:, 1], u1, rtol=1e-9, atol=1e-9)
    np.testing.assert_allclose(U[:, 2], u2, rtol=1e-9, atol=1e-9)


def test_dynamic_matrix_sparsity_signature_stable() -> None:
    n = 3
    K = eye(n, format="csr") * 2.0
    M = eye(n, format="csr") * 0.5
    C1 = eye(n, format="csr") * 0.1
    C2 = eye(n, format="csr") * 0.2
    A0 = _dynamic_matrix(K, M, C1, 1.0)
    A1 = _dynamic_matrix(K, M, C2, 2.0)
    assert dynamic_matrix_csc_signature(A0) == dynamic_matrix_csc_signature(A1)


def test_modal_zeta_per_mode_scalars() -> None:
    n = 3
    K = eye(n, format="csr") * 10.0
    M = eye(n, format="csr") * 0.5
    om, Phi = undamped_natural_modes(K, M, 3)
    F = np.array([0.0, 1.0, 0.0])
    fhz = np.array([0.5])
    Um = harmonic_displacement_modal_superposition(om, Phi, F, fhz, 0.02, zeta_per_mode=np.full(3, 0.02))
    U0 = harmonic_displacement_modal_superposition(om, Phi, F, fhz, 0.02)
    np.testing.assert_allclose(Um, U0, rtol=1e-10, atol=1e-12)


def test_modal_truncation_error_small_system() -> None:
    n = 5
    K = eye(n, format="csr") * 5.0
    M = eye(n, format="csr") * 0.2
    z = 0.0
    omega_ref = 4.0
    C = harmonic_damping_matrix(M, K, z, omega_ref, 0.0, 0.0)
    rng = np.random.default_rng(0)
    F = rng.standard_normal(n)
    freqs = np.array([1.0])
    w = 2.0 * np.pi * float(freqs[0])
    u_dir = solve_one_frequency(K, M, C, F, w)
    om, Phi = undamped_natural_modes(K, M, 3)
    u_mod = harmonic_displacement_modal_superposition(om, Phi, F, freqs, z)
    mean_r, max_r = harmonic_truncation_metrics_vs_direct(u_mod.reshape(-1, 1), u_dir.reshape(-1, 1))
    assert max_r > 0.01


def test_load_zeta_table_roundtrip(tmp_path) -> None:
    p = tmp_path / "z.txt"
    p.write_text("#hz zeta\n10 0.02\n100 0.05\n", encoding="utf-8")
    fh, zt = load_zeta_vs_frequency_hz(str(p))
    xs = np.array([50.0])
    zi = interpolate_zeta_hz(xs, fh, zt)
    assert zi.shape == (1,)
    assert 0.02 < zi[0] < 0.05


def test_mp_current_frequency_changes_response_vs_geometric() -> None:
    n = 2
    K = eye(n, format="csr") * 100.0
    M = eye(n, format="csr") * 1.0
    zeta = 0.05
    omega_geo = 2.0 * np.pi * np.sqrt(10.0 * 50.0)
    C_geo = harmonic_damping_matrix(M, K, zeta, omega_geo, 0.0, 0.0)
    F = np.array([1.0, 0.0])
    freqs = np.array([10.0])
    U_geo = sweep_displacements(K, M, C_geo, F, freqs)
    U_cur = sweep_displacements(
        K,
        M,
        C_geo,
        F,
        freqs,
        mp_damping_reference="current_frequency",
        zeta_mp=zeta,
        rayleigh_alpha_mp=0.0,
        rayleigh_beta_mp=0.0,
    )
    assert np.linalg.norm(U_geo - U_cur) > 1e-9


def test_modal_superposition_matches_direct_undamped_full_basis() -> None:
    n = 4
    K = eye(n, format="csr") * 100.0
    M = eye(n, format="csr") * 1.0
    zeta = 0.0
    omega_ref = 10.0
    C = harmonic_damping_matrix(M, K, zeta, omega_ref, 0.0, 0.0)
    rng = np.random.default_rng(0)
    F = rng.standard_normal(n)
    freqs = np.array([3.0])
    omega_t = 2.0 * np.pi * float(freqs[0])
    u_direct = solve_one_frequency(K, M, C, F, omega_t)
    om_nat, Phi = undamped_natural_modes(K, M, n)
    u_modal = harmonic_displacement_modal_superposition(om_nat, Phi, F, freqs, zeta)
    np.testing.assert_allclose(u_modal.ravel(), u_direct, rtol=1e-8, atol=1e-8)


def test_load_phase_rotates_response() -> None:
    n = 2
    K = eye(n, format="csr") * 20.0
    M = eye(n, format="csr") * 1.0
    C = eye(n, format="csr") * 0.1
    F = np.array([1.0, 0.0])
    omega = 2.0
    u0 = solve_one_frequency(K, M, C, F, omega)
    phase = 0.7
    u1 = solve_one_frequency(K, M, C, F * np.exp(1j * phase), omega)
    np.testing.assert_allclose(u1, u0 * np.exp(1j * phase), rtol=1e-9, atol=1e-9)


@pytest.mark.integration
def test_job_smoke_harmonic_process_job() -> None:
    job_dir = PROJECT_ROOT / "jobs" / "job_smoke_harmonic"
    assert job_dir.is_dir()
    res_dir = setup_job_results_directory("pytest_harmonic_smoke_e2e")
    jt: dict = {}
    je: dict = {}
    process_job(
        str(job_dir),
        res_dir,
        jt,
        je,
        force_serial=True,
        max_processes_per_job=1,
    )
    root = Path(res_dir)
    assert (root / "logs" / "run_manifest.json").is_file()
    harm_dir = root / "primary_results" / "harmonic_results"
    assert (harm_dir / "job_smoke_harmonic_frequencies_hz.txt").is_file()
    assert (harm_dir / "job_smoke_harmonic_displacement_real.txt").is_file()
    manifest = root / "logs" / "primary_artifacts.json"
    assert manifest.is_file(), "harmonic runner should write primary_artifacts.json like static/spectral"
    payload = json.loads(manifest.read_text(encoding="utf-8"))
    assert payload.get("family") == "harmonic"
    arts = payload.get("artifacts") or {}
    assert "frequencies_hz" in arts
    assert arts["frequencies_hz"].startswith("primary_results/harmonic_results/")


def test_build_harmonic_complex_snapshots_both_components() -> None:
    Ucol = np.array([1.0 + 2.0j, 3.0 + 4.0j], dtype=np.complex128)
    snaps = build_harmonic_complex_snapshots(
        U_column=Ucol,
        frequency_index=2,
        displacement_component="complex_components",
    )
    assert len(snaps) == 2
    assert snaps[0].component_name == "real"
    assert snaps[1].component_name == "imag"
    np.testing.assert_allclose(snaps[0].U_global, np.array([1.0, 3.0]))
    np.testing.assert_allclose(snaps[1].U_global, np.array([2.0, 4.0]))

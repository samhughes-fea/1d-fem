# processing/harmonic/frequency_response.py
"""Frequency-domain kernels for §4 harmonic / steady-state response."""

from __future__ import annotations

import os
from concurrent.futures import ThreadPoolExecutor
from typing import Optional

import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import spsolve, splu

def mass_proportional_damping(M_mod: csr_matrix, zeta: float, omega_ref_rad: float) -> csr_matrix:
    """
    Mass-proportional damping: C = 2 * zeta * omega_ref * M.

    Parameters
    ----------
    M_mod : csr_matrix
        Mass matrix after boundary conditions (same layout as K_mod).
    zeta : float
        Modal damping ratio (dimensionless).
    omega_ref_rad : float
        Reference angular frequency (rad/s) for scaling C.
    """
    scale = 2.0 * float(zeta) * float(omega_ref_rad)
    return M_mod * scale


def rayleigh_damping(M_mod: csr_matrix, K_mod: csr_matrix, alpha: float, beta: float) -> csr_matrix:
    """
    Rayleigh damping contribution: alpha * M + beta * K (alpha [1/s], beta [s]).
    """
    a = float(alpha)
    b = float(beta)
    if a == 0.0 and b == 0.0:
        return M_mod * 0.0
    out = M_mod * a
    if b != 0.0:
        out = out + K_mod * b
    return out


def harmonic_damping_matrix(
    M_mod: csr_matrix,
    K_mod: csr_matrix,
    zeta: float,
    omega_ref_rad: float,
    rayleigh_alpha: float,
    rayleigh_beta: float,
) -> csr_matrix:
    """
    Combined damping: mass-proportional reference scaling plus optional Rayleigh terms.

    C = 2*zeta*omega_ref*M + rayleigh_alpha*M + rayleigh_beta*K
    """
    C = mass_proportional_damping(M_mod, zeta, omega_ref_rad)
    if float(rayleigh_alpha) != 0.0 or float(rayleigh_beta) != 0.0:
        C = C + rayleigh_damping(M_mod, K_mod, rayleigh_alpha, rayleigh_beta)
    return C


def frequency_grid_hz(f_min_hz: float, f_max_hz: float, n_points: int) -> np.ndarray:
    """Uniformly spaced frequencies in Hz (inclusive endpoints)."""
    if n_points < 2:
        raise ValueError("num_frequency_points must be >= 2")
    return np.linspace(float(f_min_hz), float(f_max_hz), int(n_points))


def geometric_mean_reference_omega_rad(frequencies_hz: np.ndarray) -> float:
    """``2π √(f_min · f_max)`` for a sweep vector (matches harmonic runner reference)."""
    f = np.asarray(frequencies_hz, dtype=np.float64).ravel()
    if f.size < 1:
        raise ValueError("frequencies_hz must be non-empty")
    return 2.0 * np.pi * float(np.sqrt(float(f[0]) * float(f[-1])))


def _dynamic_matrix(
    K_mod: csr_matrix,
    M_mod: csr_matrix,
    C_mod: csr_matrix,
    omega_rad: float,
) -> csr_matrix:
    omega = float(omega_rad)
    om2 = omega * omega
    Kc = K_mod.astype(np.complex128)
    Mc = M_mod.astype(np.complex128)
    Cc = C_mod.astype(np.complex128)
    return Kc + (-om2) * Mc + (1j * omega) * Cc


def dynamic_matrix_csc_signature(A: csr_matrix) -> tuple:
    """Hashable sparsity signature for **A** in CSC form (shape, nnz, structure)."""
    Ac = A.tocsc()
    return (Ac.shape, Ac.nnz, Ac.indices.tobytes(), Ac.indptr.tobytes())


def _splu_csc_buffer_enabled() -> bool:
    """Reuse one CSC matrix shell for ``splu`` when sparsity is fixed (``A(ω)`` same pattern)."""
    v = os.environ.get("FEM_HARMONIC_SPLU_CSC_BUFFER", "").strip().lower()
    if v in ("0", "false", "no", "off"):
        return False
    if v in ("1", "true", "yes", "on"):
        return True
    return True


def _solve_complex_sparse_system(
    A: csr_matrix,
    b: np.ndarray,
    *,
    linear_solver: str,
    sparsity_signature_cache: Optional[list] = None,
    splu_csc_workspace: Optional[dict] = None,
) -> np.ndarray:
    """
    Solve A x = b with A complex sparse; ``linear_solver`` is ``spsolve`` or ``splu``.

    When environment ``FEM_HARMONIC_VERIFY_A_PATTERN=1`` is set, or when
    *sparsity_signature_cache* is a list, the CSC structure of *A* is recorded on the
    first solve and later solves must match (SuperLU reuse of symbolic factorization is
    not exposed in SciPy's high-level ``splu``; this hook validates fixed sparsity across ω).

    When *splu_csc_workspace* is a dict and ``FEM_HARMONIC_SPLU_CSC_BUFFER`` is enabled
    (default *on*), *splu* reuses a single allocated CSC matrix: only ``.data`` is
    overwritten each frequency. **Numeric** factorization is still recomputed every step;
    this reduces Python alloc / conversion overhead, not SuperLU symbolics.
    Set ``FEM_HARMONIC_SPLU_CSC_BUFFER=0`` to disable.
    """
    b = np.asarray(b, dtype=np.complex128).ravel()
    if sparsity_signature_cache is not None:
        sig = dynamic_matrix_csc_signature(A)
        if not sparsity_signature_cache:
            sparsity_signature_cache.append(sig)
        elif sparsity_signature_cache[0] != sig:
            raise RuntimeError(
                "harmonic: dynamic matrix CSC sparsity changed between frequency samples "
                "(symbolic factor reuse requires a fixed pattern)"
            )

    if linear_solver == "splu":
        Acsc = A.tocsc()
        use_buf = (
            splu_csc_workspace is not None
            and _splu_csc_buffer_enabled()
        )
        if use_buf:
            sig = dynamic_matrix_csc_signature(A)
            if splu_csc_workspace.get("signature") is None:
                splu_csc_workspace["signature"] = sig
                splu_csc_workspace["csc"] = Acsc.copy()
            elif splu_csc_workspace["signature"] != sig:
                raise RuntimeError(
                    "harmonic: splu CSC buffer sparsity mismatch (disable FEM_HARMONIC_SPLU_CSC_BUFFER)"
                )
            else:
                splu_csc_workspace["csc"].data[:] = Acsc.data
            return splu(splu_csc_workspace["csc"]).solve(b)
        return splu(Acsc).solve(b)
    if linear_solver == "spsolve":
        return spsolve(A, b)
    raise ValueError(f"linear_solver must be 'spsolve' or 'splu', got {linear_solver!r}")


def solve_one_frequency(
    K_mod: csr_matrix,
    M_mod: csr_matrix,
    C_mod: csr_matrix,
    F_mod: np.ndarray,
    omega_rad: float,
    *,
    linear_solver: str = "spsolve",
    sparsity_signature_cache: Optional[list] = None,
    splu_csc_workspace: Optional[dict] = None,
) -> np.ndarray:
    """
    Solve (-omega^2 M + i omega C + K) u = F for complex u.

    All matrices are real-valued; the system matrix is formed in complex arithmetic.
    *F_mod* may be real or complex (e.g. phased harmonic loads).

    *linear_solver*: ``spsolve`` (default) or ``splu`` (explicit SuperLU; useful when benchmarking).
    """
    A = _dynamic_matrix(K_mod, M_mod, C_mod, omega_rad)
    b = np.asarray(F_mod, dtype=np.complex128).ravel()
    return _solve_complex_sparse_system(
        A,
        b,
        linear_solver=linear_solver,
        sparsity_signature_cache=sparsity_signature_cache,
        splu_csc_workspace=splu_csc_workspace,
    )


def solve_one_frequency_partitioned(
    K_mod: csr_matrix,
    M_mod: csr_matrix,
    C_mod: csr_matrix,
    F_mod: np.ndarray,
    omega_rad: float,
    u_prescribed: np.ndarray,
    dof_unknown: np.ndarray,
    *,
    linear_solver: str = "spsolve",
    sparsity_signature_cache: Optional[list] = None,
    splu_csc_workspace: Optional[dict] = None,
) -> np.ndarray:
    """
    Solve A u = F with prescribed complex displacements on DOFs not listed in *dof_unknown*.

    Unknown vector u satisfies u[k] = u_prescribed[k] for k not in dof_unknown
    (typically fixed + motion-prescribed DOFs set in *u_prescribed*).
    """
    A = _dynamic_matrix(K_mod, M_mod, C_mod, omega_rad)
    u_t = np.asarray(u_prescribed, dtype=np.complex128).ravel()
    if u_t.shape[0] != K_mod.shape[0]:
        raise ValueError("u_prescribed length must match global DOF count")
    Fv = np.asarray(F_mod, dtype=np.complex128).ravel()
    rhs = Fv - A @ u_t
    unk = np.asarray(dof_unknown, dtype=np.int64).ravel()
    if unk.size == 0:
        return u_t.copy()
    Auu = A[unk, :][:, unk]
    rhs_u = rhs[unk]
    u_sol = _solve_complex_sparse_system(
        Auu,
        rhs_u,
        linear_solver=linear_solver,
        sparsity_signature_cache=sparsity_signature_cache,
        splu_csc_workspace=splu_csc_workspace,
    )
    out = u_t.copy()
    out[unk] = u_sol
    return out


def sweep_displacements(
    K_mod: csr_matrix,
    M_mod: csr_matrix,
    C_mod: csr_matrix,
    F_mod: np.ndarray,
    frequencies_hz: np.ndarray,
    *,
    u_prescribed: Optional[np.ndarray] = None,
    dof_unknown: Optional[np.ndarray] = None,
    parallel: bool = False,
    max_workers: Optional[int] = None,
    mp_damping_reference: str = "geometric_mean",
    zeta_mp: float = 0.0,
    rayleigh_alpha_mp: float = 0.0,
    rayleigh_beta_mp: float = 0.0,
    linear_solver: str = "spsolve",
    zeta_per_frequency: Optional[np.ndarray] = None,
    geometric_omega_ref_rad: Optional[float] = None,
) -> np.ndarray:
    """
    Column k is u(2*pi*f_k) for each frequency in *frequencies_hz*.

    If *dof_unknown* is provided, *u_prescribed* must hold the template displacement
    (known nonzero motion DOFs; zeros elsewhere) and only *dof_unknown* rows are solved.

    Parameters
    ----------
    parallel : bool
        When True, frequency samples are solved concurrently (full or partitioned path).
    max_workers : int, optional
        Thread pool size when *parallel* is True (defaults to min(32, n_freq)).
    mp_damping_reference : str
        ``geometric_mean`` (default): use a single *C_mod* for all frequencies (assembled by caller).
        ``current_frequency``: rebuild **C** each sample using the **current** angular frequency
        in the mass-proportional term (``2*zeta*omega*M``) plus Rayleigh terms from *zeta_mp*,
        *rayleigh_alpha_mp*, *rayleigh_beta_mp* (same convention as ``harmonic_damping_matrix``).
    linear_solver : str
        ``spsolve`` or ``splu`` — passed to each frequency solve.
    zeta_per_frequency : ndarray, optional
        If set, length must match the number of frequency samples. Mass-proportional
        ``ζ`` at each sample (from a table). Rebuilds **C** every column; combine with
        ``mp_damping_reference`` as for scalar ``zeta_mp``.
    geometric_omega_ref_rad : float, optional
        Reference ``ω`` for the MP term when ``mp_damping_reference=geometric_mean`` and
        a per-frequency ζ table is used. Defaults to ``2π√(f[0]·f[-1])`` from *frequencies_hz*.
    """
    f = np.asarray(frequencies_hz, dtype=np.float64).ravel()
    n = int(f.size)
    n_dof = int(K_mod.shape[0])
    out = np.zeros((n_dof, n), dtype=np.complex128)
    F = np.asarray(F_mod, dtype=np.complex128).ravel()
    variable_mp = mp_damping_reference == "current_frequency"
    if mp_damping_reference not in ("geometric_mean", "current_frequency"):
        raise ValueError(
            "mp_damping_reference must be 'geometric_mean' or 'current_frequency'"
        )
    zpf = None
    if zeta_per_frequency is not None:
        zpf = np.asarray(zeta_per_frequency, dtype=np.float64).ravel()
        if zpf.size != n:
            raise ValueError("zeta_per_frequency must have length equal to number of frequency samples")
        if np.any(zpf < 0.0):
            raise ValueError("zeta_per_frequency values must be non-negative")
    omega_ref_geom = (
        float(geometric_omega_ref_rad)
        if geometric_omega_ref_rad is not None
        else geometric_mean_reference_omega_rad(f)
    )

    pat_cache: Optional[list] = None
    if os.environ.get("FEM_HARMONIC_VERIFY_A_PATTERN", "").strip().lower() in (
        "1",
        "true",
        "yes",
        "on",
    ):
        pat_cache = []

    splu_ws: Optional[dict] = None
    if linear_solver == "splu" and not parallel and _splu_csc_buffer_enabled():
        splu_ws = {}

    partitioned = dof_unknown is not None
    if partitioned:
        if u_prescribed is None:
            raise ValueError("u_prescribed is required when dof_unknown is set")
        u_t = np.asarray(u_prescribed, dtype=np.complex128).ravel()
        unk = np.asarray(dof_unknown, dtype=np.int64).ravel()

        def column(k: int) -> np.ndarray:
            omega = 2.0 * np.pi * float(f[k])
            if zpf is not None:
                zk = float(zpf[k])
                omega_mp = omega if variable_mp else omega_ref_geom
                C_use = harmonic_damping_matrix(
                    M_mod, K_mod, zk, omega_mp, rayleigh_alpha_mp, rayleigh_beta_mp
                )
            elif variable_mp:
                C_use = harmonic_damping_matrix(
                    M_mod, K_mod, zeta_mp, omega, rayleigh_alpha_mp, rayleigh_beta_mp
                )
            else:
                C_use = C_mod
            return solve_one_frequency_partitioned(
                K_mod,
                M_mod,
                C_use,
                F,
                omega,
                u_t,
                unk,
                linear_solver=linear_solver,
                sparsity_signature_cache=pat_cache,
            )

    else:

        def column(k: int) -> np.ndarray:
            omega = 2.0 * np.pi * float(f[k])
            if zpf is not None:
                zk = float(zpf[k])
                omega_mp = omega if variable_mp else omega_ref_geom
                C_use = harmonic_damping_matrix(
                    M_mod, K_mod, zk, omega_mp, rayleigh_alpha_mp, rayleigh_beta_mp
                )
            elif variable_mp:
                C_use = harmonic_damping_matrix(
                    M_mod, K_mod, zeta_mp, omega, rayleigh_alpha_mp, rayleigh_beta_mp
                )
            else:
                C_use = C_mod
            return solve_one_frequency(
                K_mod,
                M_mod,
                C_use,
                F,
                omega,
                linear_solver=linear_solver,
                sparsity_signature_cache=pat_cache,
                splu_csc_workspace=splu_ws,
            )

    if not parallel:
        for k in range(n):
            out[:, k] = column(k)
        return out

    workers = max_workers if max_workers is not None else min(32, max(1, n))
    with ThreadPoolExecutor(max_workers=workers) as pool:
        cols = list(pool.map(column, range(n)))
    for k in range(n):
        out[:, k] = cols[k]
    return out

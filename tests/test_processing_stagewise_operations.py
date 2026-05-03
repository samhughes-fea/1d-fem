"""Unit tests for stagewise processing.spectral / harmonic / dynamic operations."""

from __future__ import annotations

import tempfile

import numpy as np
import pytest
from scipy.sparse import coo_matrix, csr_matrix, eye

from processing.dynamic.operations import (
    AssembleDynamicGlobalSystem,
    IntegrateTransientSystem,
    ModifyDynamicGlobalSystem,
)
from processing.harmonic.operations import (
    AssembleHarmonicLoadVector,
    AssembleHarmonicStructuralMatrices,
    BuildHarmonicDampingMatrix,
    ModifyHarmonicStructuralMatrices,
)
from processing.spectral.operations import (
    AssembleSpectralGlobalSystem,
    ModifySpectralGlobalSystem,
    PrepareSpectralLocalMatrices,
    SolveGeneralizedEigenproblem,
)
from tests.test_modal_processing import _MockElement


def test_prepare_spectral_local_matrices_dense_to_coo():
    Ke = np.eye(4, dtype=np.float64)
    Me = np.eye(4, dtype=np.float64) * 0.5
    prep = PrepareSpectralLocalMatrices(
        element_stiffness_matrices=[Ke],
        element_mass_matrices=[Me],
        job_results_dir=None,
    )
    k_list, m_list = prep.run()
    assert len(k_list) == 1 and len(m_list) == 1
    assert isinstance(k_list[0], coo_matrix)
    assert k_list[0].shape == (4, 4)


def test_modify_spectral_global_system_penalty():
    n = 8
    K = csr_matrix(eye(n, format="csr", dtype=np.float64))
    M = csr_matrix(eye(n, format="csr", dtype=np.float64))
    mod = ModifySpectralGlobalSystem(job_results_dir=None, fixed_dofs=np.arange(2))
    Km, Mm, bc = mod.run(K, M)
    assert Km.shape == (n, n) and Mm.shape == (n, n)
    assert bc.size == 2


def test_solve_generalized_eigenproblem_2x2():
    n = 4
    K = csr_matrix(eye(n, format="csr", dtype=np.float64) * 4.0)
    M = csr_matrix(eye(n, format="csr", dtype=np.float64))
    ev, vec, fh = SolveGeneralizedEigenproblem(
        num_modes=2, context="test", dense_threshold=512, job_results_dir=None
    ).run(K, M)
    assert ev.shape == (2,) and vec.shape == (n, 2) and fh.shape == (2,)


def test_modify_harmonic_structural_matrices():
    n = 6
    K = csr_matrix(eye(n, format="csr", dtype=np.float64))
    M = csr_matrix(eye(n, format="csr", dtype=np.float64))
    F = np.ones(n, dtype=np.float64)
    with tempfile.TemporaryDirectory() as td:
        Km, Mm, bc, Fm = ModifyHarmonicStructuralMatrices(
            prescribed_displacements=None,
            job_results_dir=td,
        ).run(K, M, F)
    assert Fm.shape == (n,)
    if bc.size:
        assert np.all(Fm[bc] == 0.0)


def test_build_harmonic_damping_matrix():
    n = 5
    K = csr_matrix(eye(n, format="csr", dtype=np.float64))
    M = csr_matrix(eye(n, format="csr", dtype=np.float64))
    C = BuildHarmonicDampingMatrix(job_results_dir=None).run(
        M, K, zeta=0.02, omega_ref=1.0, rayleigh_alpha=0.0, rayleigh_beta=0.0
    )
    assert C.shape == (n, n)


def test_modify_dynamic_global_system():
    n = 4
    K = csr_matrix(eye(n, format="csr", dtype=np.float64))
    M = csr_matrix(eye(n, format="csr", dtype=np.float64))
    Km, Mm, Cm, bc = ModifyDynamicGlobalSystem(
        fixed_dofs=[0, 1],
        prescribed_displacements=None,
        job_results_dir=None,
    ).run(K, M, None)
    assert Cm is None
    assert bc.size == 2


def test_assemble_dynamic_global_system_requires_elements():
    with pytest.raises(ValueError, match="No elements"):
        AssembleDynamicGlobalSystem(
            elements=[],
            element_stiffness_matrices=None,
            element_mass_matrices=None,
            total_dof=3,
            job_results_dir=None,
        ).run()


def test_integrate_transient_system_sdof():
    n = 2
    K = csr_matrix(eye(n, format="csr", dtype=np.float64) * 100.0)
    M = csr_matrix(eye(n, format="csr", dtype=np.float64))
    t_grid = np.linspace(0.0, 0.02, 5)
    f = np.array([1.0, 0.0], dtype=np.float64)

    def F_func(_t: float):
        return f.copy()

    U, V, A = IntegrateTransientSystem(
        t_grid=t_grid, force_func=F_func, job_results_dir=None
    ).run(K, M, None, np.zeros(n), np.zeros(n))
    assert U.shape == (len(t_grid), n) and V.shape == U.shape and A.shape == U.shape


def test_assemble_spectral_global_system_single_element():
    total_dof = 12
    dof = np.arange(12, dtype=np.int32)
    elements = [_MockElement(0, dof)]
    Ke = np.eye(12, dtype=np.float64)
    Me = np.eye(12, dtype=np.float64) * 0.5
    ke_list, me_list = PrepareSpectralLocalMatrices(
        element_stiffness_matrices=[Ke],
        element_mass_matrices=[Me],
        job_results_dir=None,
    ).run()
    K, M, lm = AssembleSpectralGlobalSystem(
        elements=elements,
        element_stiffness_matrices=ke_list,
        element_mass_matrices=me_list,
        total_dof=total_dof,
        job_results_dir=None,
    ).run()
    assert K.shape == (12, 12) and M.shape == (12, 12)
    assert len(lm) == 1


def test_assemble_harmonic_structural_and_load_vector():
    total_dof = 12
    dof = np.arange(12, dtype=np.int32)
    elements = [_MockElement(0, dof)]
    Ke = np.eye(12, dtype=np.float64)
    Me = np.eye(12, dtype=np.float64) * 0.25
    Fe = np.ones(12, dtype=np.float64)
    ke_list, me_list = PrepareSpectralLocalMatrices(
        element_stiffness_matrices=[Ke],
        element_mass_matrices=[Me],
        job_results_dir=None,
    ).run()
    K, M, _ = AssembleHarmonicStructuralMatrices(
        elements=elements,
        element_stiffness_matrices=ke_list,
        element_mass_matrices=me_list,
        total_dof=total_dof,
        job_results_dir=None,
    ).run()
    F = AssembleHarmonicLoadVector(
        elements=elements,
        element_stiffness_matrices=ke_list,
        element_force_vectors=[Fe],
        total_dof=total_dof,
        job_results_dir=None,
    ).run()
    assert F.shape == (12,)
    assert K.nnz > 0 and M.nnz > 0

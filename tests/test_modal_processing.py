"""
Unit tests for processing/modal assembly and boundary conditions.
No imports from processing.static; tests use small hand-built K, M.
"""

import numpy as np
import pytest
from scipy.sparse import csr_matrix

from processing.eigen.assembly import assemble_global_matrices
from processing.eigen.boundary_conditions import apply_boundary_conditions


class _MockElement:
    """Minimal element with assemble_global_dof_indices for a single 2-node bar (12 DOFs)."""
    def __init__(self, element_id, dof_indices):
        self.element_id = element_id
        self._dof = np.asarray(dof_indices, dtype=np.int32)

    def assemble_global_dof_indices(self):
        return self._dof


class TestModalAssembly:
    """Tests for processing.eigen.assembly.assemble_global_matrices."""

    def test_assemble_single_element(self):
        """Single element: 12 DOFs, K and M 12x12 -> global 12x12."""
        total_dof = 12
        dof = np.arange(12, dtype=np.int32)
        elements = [_MockElement(0, dof)]
        Ke = np.eye(12, dtype=np.float64)
        Me = np.eye(12, dtype=np.float64) * 0.5
        K_global, M_global, dof_map = assemble_global_matrices(
            elements=elements,
            element_stiffness_matrices=[Ke],
            element_mass_matrices=[Me],
            total_dof=total_dof,
        )
        assert K_global.shape == (12, 12)
        assert M_global.shape == (12, 12)
        np.testing.assert_allclose(K_global.toarray(), np.eye(12))
        np.testing.assert_allclose(M_global.toarray(), np.eye(12) * 0.5)
        assert len(dof_map) == 1
        np.testing.assert_array_equal(dof_map[0], dof)

    def test_assemble_two_elements_no_overlap(self):
        """Two elements covering DOFs 0-11 and 12-23."""
        total_dof = 24
        elements = [
            _MockElement(0, np.arange(0, 12)),
            _MockElement(1, np.arange(12, 24)),
        ]
        Ke = np.eye(12, dtype=np.float64) * 2.0
        Me = np.eye(12, dtype=np.float64) * 0.25
        K_global, M_global, dof_map = assemble_global_matrices(
            elements=elements,
            element_stiffness_matrices=[Ke, Ke],
            element_mass_matrices=[Me, Me],
            total_dof=total_dof,
        )
        assert K_global.shape == (24, 24)
        assert M_global.shape == (24, 24)
        expected_k = np.zeros((24, 24))
        expected_k[:12, :12] = 2.0 * np.eye(12)
        expected_k[12:, 12:] = 2.0 * np.eye(12)
        expected_m = np.zeros((24, 24))
        expected_m[:12, :12] = 0.25 * np.eye(12)
        expected_m[12:, 12:] = 0.25 * np.eye(12)
        np.testing.assert_allclose(K_global.toarray(), expected_k)
        np.testing.assert_allclose(M_global.toarray(), expected_m)


class TestModalBoundaryConditions:
    """Tests for processing.eigen.boundary_conditions.apply_boundary_conditions."""

    def test_apply_fixed_dofs(self):
        """Fixed first 6 DOFs: K_mod and M_mod have penalty/1 on diagonal for those DOFs."""
        n = 12
        K = csr_matrix(np.eye(n) * 10.0)
        M = csr_matrix(np.eye(n) * 1.0)
        K_mod, M_mod, bc_dofs = apply_boundary_conditions(K, M, fixed_dofs=[0, 1, 2, 3, 4, 5])
        np.testing.assert_array_equal(bc_dofs, [0, 1, 2, 3, 4, 5])
        K_dense = K_mod.toarray()
        M_dense = M_mod.toarray()
        for i in range(6):
            assert K_dense[i, i] > 1e30
            assert abs(M_dense[i, i] - 1.0) < 1e-9
            assert np.all(K_dense[i, :] == 0) or (K_dense[i, i] > 1e30)
            assert np.all(K_dense[:, i] == 0) or (K_dense[i, i] > 1e30)
        for i in range(6, 12):
            assert abs(K_dense[i, i] - 10.0) < 1e-9
            assert abs(M_dense[i, i] - 1.0) < 1e-9

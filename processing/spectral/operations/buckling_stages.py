# processing/spectral/operations/buckling_stages.py
"""Staged Section 5 linear buckling matrix operations (after prestress displacement is known)."""

from __future__ import annotations

from pathlib import Path
from typing import List, Optional, Sequence, Tuple, Union

import numpy as np
from scipy.sparse import csr_matrix

from processing.buckling import (
    apply_buckling_boundary_conditions,
    assemble_global_geometric_stiffness,
    solve_linear_buckling_eigenpairs,
)
from processing.spectral.spectral_diagnostics import (
    log_spectral_constrained_dofs,
    log_spectral_diagnostics,
)

from processing.common.stage_logging import init_stage_logger


class AssembleBucklingGeometricStiffness:
    """Build global geometric stiffness from prestress displacements."""

    def __init__(
        self,
        elements: List,
        U_global: np.ndarray,
        total_dof: int,
        job_results_dir: Optional[Union[str, Path]] = None,
    ):
        self.elements = list(elements)
        self.U_global = np.asarray(U_global, dtype=np.float64).ravel()
        self.total_dof = int(total_dof)
        self.job_results_dir = Path(job_results_dir) if job_results_dir else None
        self._log = init_stage_logger(
            "AssembleBucklingGeometricStiffness", self.job_results_dir
        )

    def run(self) -> csr_matrix:
        self._log.info("Assembling global geometric stiffness (K_sigma)")
        Kg = assemble_global_geometric_stiffness(
            self.elements, self.U_global, self.total_dof
        )
        self._log.info("K_sigma nnz=%s", Kg.nnz)
        return Kg


class ModifyBucklingGlobalMatrices:
    """Apply buckling BCs to elastic K and geometric K_sigma."""

    def __init__(
        self,
        job_results_dir: Optional[Union[str, Path]] = None,
        prescribed_displacements: Optional[dict] = None,
        fixed_dofs: Optional[Sequence[int]] = None,
    ):
        self.job_results_dir = Path(job_results_dir) if job_results_dir else None
        self.prescribed_displacements = prescribed_displacements
        self.fixed_dofs = fixed_dofs
        self._log = init_stage_logger("ModifyBucklingGlobalMatrices", self.job_results_dir)

    def run(
        self,
        K_global: csr_matrix,
        Kg_global: csr_matrix,
    ) -> Tuple[csr_matrix, csr_matrix, np.ndarray]:
        self._log.info("Applying buckling boundary conditions")
        K_mod, Kg_mod, bc_dofs = apply_buckling_boundary_conditions(
            K_global,
            Kg_global,
            prescribed_displacements=self.prescribed_displacements,
            fixed_dofs=self.fixed_dofs,
        )
        log_spectral_diagnostics(
            K_mod,
            Kg_mod,
            str(self.job_results_dir) if self.job_results_dir else None,
        )
        log_spectral_constrained_dofs(
            bc_dofs,
            total_dof=int(K_mod.shape[0]),
            job_results_dir=str(self.job_results_dir) if self.job_results_dir else None,
        )
        self._log.info("Buckling BCs applied: %s constrained DOFs", bc_dofs.size)
        return K_mod, Kg_mod, bc_dofs


class SolveLinearBucklingEigenpairs:
    """Linear buckling eigenproblem on (K, K_sigma) pencil."""

    def __init__(
        self,
        num_modes: int,
        job_results_dir: Optional[Union[str, Path]] = None,
    ):
        self.num_modes = int(num_modes)
        self.job_results_dir = Path(job_results_dir) if job_results_dir else None
        self._log = init_stage_logger(
            "SolveLinearBucklingEigenpairs", self.job_results_dir
        )

    def run(
        self,
        K_mod: csr_matrix,
        Kg_mod: csr_matrix,
        constrained_dofs: np.ndarray,
    ) -> Tuple[np.ndarray, np.ndarray]:
        self._log.info("Solving linear buckling eigenpairs: num_modes=%s", self.num_modes)
        lambdas, modes = solve_linear_buckling_eigenpairs(
            K_mod, Kg_mod, self.num_modes, constrained_dofs=constrained_dofs
        )
        self._log.info("Buckling eigen solve finished")
        return lambdas, modes

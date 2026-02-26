# processing\static\results\compute_tertiary\tertiary_results_orchestrator.py

import numpy as np
from pathlib import Path
from typing import Optional
import logging

from .section_force import ComputeSectionForce
from .principal_stress import ComputePrincipalStress
from .integrated_elemental_results import ComputeIntegratedElementalResults
from ..containers.tertiary_results import TertiaryResults

logger = logging.getLogger(__name__)


class TertiaryResultsOrchestrator:
    """
    Orchestrates computation of tertiary (highly derived) results.

    Tertiary results are engineering quantities computed from secondary
    results, including:
    - Section force resultants [N, Vy, Vz, T, My, Mz]
    - Principal stresses [σ1, σ2, σ3]
    - Von Mises equivalent stress
    - Maximum shear stress

    These are typically used for design verification and detailed
    engineering analysis.

    Parameters
    ----------
    secondary_results : SecondaryResultSet
        Container with secondary results (strain, stress at Gauss points)
    job_results_dir : str or Path, optional
        Directory for logging and output
    """

    def __init__(
        self,
        secondary_results,
        formulation_cache=None,
        element_dictionary=None,
        grid_dictionary=None,
        job_results_dir: Optional[str | Path] = None,
    ):
        self.secondary_results = secondary_results
        self.formulation_cache = formulation_cache
        self.element_dictionary = element_dictionary
        self.grid_dictionary = grid_dictionary
        self.job_results_dir = Path(job_results_dir) if job_results_dir else None
        self.logger = self._init_logging()

    def compute(self) -> TertiaryResults:
        """
        Compute all tertiary results.

        Returns
        -------
        TertiaryResults
            Container with all tertiary result fields
        """
        self.logger.info("=" * 70)
        self.logger.info("TERTIARY RESULTS COMPUTATION")
        self.logger.info("=" * 70)

        # Extract stress data from secondary results
        stress_gauss = self.secondary_results.gaussian_results.stress

        if stress_gauss is None:
            self.logger.warning("⚠️  No Gaussian stress data available for tertiary computation")
            return TertiaryResults()

        # 1. Compute section forces
        self.logger.info("\n[1/3] Computing section force resultants...")
        section_force_computer = ComputeSectionForce(stress_gauss)
        section_forces = section_force_computer.compute()

        # 2. Compute principal stresses and stress invariants
        self.logger.info("\n[2/3] Computing principal stresses and invariants...")
        principal_computer = ComputePrincipalStress(stress_gauss)
        principal_stresses, von_mises, max_shear = principal_computer.compute_all()

        # 3. Compute integrated elemental results and nodal section forces
        total_strain_energy = None
        integrated_section_forces = None
        nodal_section_forces = None

        if self.formulation_cache is not None:
            self.logger.info("\n[3/3] Computing integrated elemental results...")
            integrated_computer = ComputeIntegratedElementalResults(
                secondary_results=self.secondary_results,
                formulation_cache=self.formulation_cache
            )
            
            # Compute total strain energy per element
            total_strain_energy = integrated_computer.compute_total_strain_energy()
            
            # Compute integrated section forces per element
            integrated_section_forces = integrated_computer.compute_integrated_section_forces(
                section_forces
            )
        else:
            self.logger.warning("⚠️  No formulation_cache provided, skipping integrated elemental results")

        # Nodal section forces (shape-function projection from GP data)
        if (
            self.formulation_cache is not None
            and self.element_dictionary is not None
            and self.grid_dictionary is not None
            and section_forces is not None
        ):
            from .nodal_section_forces_projector import NodalSectionForcesProjector
            projector = NodalSectionForcesProjector(
                section_forces_gauss=section_forces,
                formulation_cache=self.formulation_cache,
                element_dictionary=self.element_dictionary,
                grid_dictionary=self.grid_dictionary,
            )
            nodal_section_forces = projector.project()

        # Package results
        tertiary_results = TertiaryResults(
            section_forces=section_forces,
            principal_stresses=principal_stresses,
            von_mises_stress=von_mises,
            max_shear_stress=max_shear,
            total_strain_energy=total_strain_energy,
            integrated_section_forces=integrated_section_forces,
            nodal_section_forces=nodal_section_forces,
        )

        self.logger.info("\n" + "=" * 70)
        self.logger.info("✅ TERTIARY RESULTS COMPUTATION COMPLETE")
        self.logger.info("=" * 70)

        return tertiary_results

    def _init_logging(self) -> logging.Logger:
        """Initialize logger for tertiary results computation."""
        logger = logging.getLogger(f"TertiaryResultsOrchestrator.{id(self)}")
        logger.handlers.clear()
        logger.setLevel(logging.DEBUG)
        logger.propagate = False

        log_path = None
        if self.job_results_dir:
            logs_dir = self.job_results_dir.parent / "logs"
            logs_dir.mkdir(parents=True, exist_ok=True)
            log_path = logs_dir / "TertiaryResultsOrchestrator.log"
            try:
                file_handler = logging.FileHandler(log_path, mode="w", encoding="utf-8")
                file_handler.setFormatter(logging.Formatter(
                    "%(asctime)s [%(levelname)s] %(message)s "
                    "(Module: %(module)s, Line: %(lineno)d)"
                ))
                logger.addHandler(file_handler)
            except Exception as e:
                print(f"⚠️ Failed to create file handler for TertiaryResultsOrchestrator class log: {e}")

        stream_handler = logging.StreamHandler()
        stream_handler.setLevel(logging.INFO)
        stream_handler.setFormatter(logging.Formatter(
            "%(asctime)s - TERTIARY - %(levelname)s - %(message)s",
            datefmt="%H:%M:%S"
        ))
        logger.addHandler(stream_handler)

        if log_path:
            logger.debug(f"📁 Log file created at: {log_path}")

        return logger

    def compute_summary_statistics(self, tertiary_results: TertiaryResults) -> dict:
        """
        Compute summary statistics for tertiary results.

        Useful for quick assessment of maximum stresses, critical locations, etc.

        Parameters
        ----------
        tertiary_results : TertiaryResults
            Computed tertiary results

        Returns
        -------
        dict
            Summary statistics including max/min values and locations
        """
        summary = {}

        # Max Von Mises stress
        if tertiary_results.von_mises_stress:
            all_von_mises = [
                vm for elem in tertiary_results.von_mises_stress for vm in elem
            ]
            summary['max_von_mises'] = max(all_von_mises)
            summary['min_von_mises'] = min(all_von_mises)

        # Max principal stress
        if tertiary_results.principal_stresses:
            all_σ1 = [
                ps[0] for elem in tertiary_results.principal_stresses for ps in elem
            ]
            summary['max_principal_stress'] = max(all_σ1)
            summary['min_principal_stress'] = min(all_σ1)

        # Max shear stress
        if tertiary_results.max_shear_stress:
            all_shear = [
                τ for elem in tertiary_results.max_shear_stress for τ in elem
            ]
            summary['max_shear_stress'] = max(all_shear)

        # Log summary
        self.logger.info("\n" + "─" * 70)
        self.logger.info("TERTIARY RESULTS SUMMARY")
        self.logger.info("─" * 70)
        for key, value in summary.items():
            self.logger.info(f"  {key:.<50} {value:.6e}")
        self.logger.info("─" * 70)

        return summary

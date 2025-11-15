# processing_OOP\static\results\save_tertiary_container.py

import numpy as np
import pandas as pd
from pathlib import Path
from typing import Optional
import logging

logger = logging.getLogger(__name__)


class SaveTertiaryResults:
    """
    Saves tertiary results (section forces, principal stresses, failure indices) to CSV.

    Tertiary results are highly derived engineering quantities used for
    design verification and failure analysis:
    - Section force resultants [N, Vy, Vz, T, My, Mz]
    - Principal stresses [σ1, σ2, σ3]
    - Von Mises equivalent stress
    - Maximum shear stress
    - Material failure indices

    Parameters
    ----------
    tertiary_results : TertiaryResults
        Container with all tertiary results
    save_dir : str or Path
        Base directory where results should be saved
    """

    def __init__(
        self,
        tertiary_results,
        save_dir: str | Path
    ):
        self.results = tertiary_results
        self.save_dir = Path(save_dir)

        # Create directory structure
        self.tertiary_dir = self.save_dir / "tertiary_results"
        self.section_forces_dir = self.tertiary_dir / "section_forces"
        self.principal_stress_dir = self.tertiary_dir / "principal_stress"

        for directory in [
            self.section_forces_dir,
            self.principal_stress_dir
        ]:
            directory.mkdir(parents=True, exist_ok=True)

    def save_all(self):
        """
        Save tertiary results to disk.
        
        Saves only the core engineering resultants:
        - Section forces [N, Vy, Vz, T, My, Mz]
        - Principal stresses [σ1, σ2, σ3]
        
        Note: Von Mises, max shear, and failure indices are computed but
        not saved separately - they're included in the summary CSV.
        """
        logger.info("=" * 70)
        logger.info("SAVING TERTIARY RESULTS")
        logger.info("=" * 70)

        saved_count = 0

        # Save section forces
        if self.results.section_forces:
            logger.info("\n[1/2] Saving section force resultants...")
            self._save_section_forces()
            saved_count += 1

        # Save principal stresses
        if self.results.principal_stresses:
            logger.info("\n[2/2] Saving principal stresses...")
            self._save_principal_stresses()
            saved_count += 1

        logger.info("\n" + "=" * 70)
        logger.info(f"✅ TERTIARY RESULTS SAVED ({saved_count} categories)")
        logger.info(f"   Location: {self.tertiary_dir}")
        logger.info(f"   Note: Von Mises, shear stress, and failure data in summary CSV")
        logger.info("=" * 70)

    def _save_section_forces(self):
        """Save section force resultants [N, Vy, Vz, T, My, Mz] per Gauss point."""
        for elem_idx, elem_section_forces in enumerate(self.results.section_forces):
            # Stack all Gauss point section forces for this element
            section_force_array = np.array(elem_section_forces)  # shape: (n_gauss, 6)
            filename = self.section_forces_dir / f"section_forces_elem_{elem_idx:06d}.csv"
            header = "N,Vy,Vz,T,My,Mz"
            np.savetxt(filename, section_force_array, delimiter=",",
                      fmt="%.12e", header=header, comments='')

        logger.info(f"   ✓ Section forces saved: {len(self.results.section_forces)} elements")

    def _save_principal_stresses(self):
        """Save principal stresses [σ1, σ2, σ3] per Gauss point."""
        for elem_idx, elem_principals in enumerate(self.results.principal_stresses):
            # Stack all Gauss point principal stresses
            principal_array = np.array(elem_principals)  # shape: (n_gauss, 3)
            filename = self.principal_stress_dir / f"principal_stress_elem_{elem_idx:06d}.csv"
            header = "σ1,σ2,σ3"
            np.savetxt(filename, principal_array, delimiter=",",
                      fmt="%.12e", header=header, comments='')

        logger.info(f"   ✓ Principal stresses saved: {len(self.results.principal_stresses)} elements")


class SaveTertiaryResultsSummary:
    """
    Generates and saves summary statistics and critical locations for tertiary results.

    Creates comprehensive summary files including:
    - Maximum/minimum stress values
    - Critical element/Gauss point locations
    - Failure status assessment
    - Design margin statistics
    """

    def __init__(
        self,
        tertiary_results,
        save_dir: str | Path,
        element_info: Optional[dict] = None
    ):
        self.results = tertiary_results
        self.save_dir = Path(save_dir)
        self.element_info = element_info or {}

    def save(self):
        """Generate and save comprehensive tertiary results summary."""
        summary_data = []

        # Collect element-wise statistics
        n_elements = len(self.results.section_forces) if self.results.section_forces else 0

        for elem_idx in range(n_elements):
            elem_summary = {'element_id': elem_idx}

            # Section force statistics
            if self.results.section_forces:
                section_forces = np.array(self.results.section_forces[elem_idx])
                elem_summary['max_axial_force'] = np.max(np.abs(section_forces[:, 0]))
                elem_summary['max_shear_force'] = np.max(np.abs(section_forces[:, 1:3]))
                elem_summary['max_torque'] = np.max(np.abs(section_forces[:, 3]))
                elem_summary['max_moment'] = np.max(np.abs(section_forces[:, 4:6]))

            # Principal stress statistics
            if self.results.principal_stresses:
                principals = np.array(self.results.principal_stresses[elem_idx])
                elem_summary['max_principal_σ1'] = np.max(principals[:, 0])
                elem_summary['min_principal_σ3'] = np.min(principals[:, 2])

            # Von Mises stress
            if self.results.von_mises_stress:
                von_mises = self.results.von_mises_stress[elem_idx]
                elem_summary['max_von_mises'] = np.max(von_mises)
                elem_summary['mean_von_mises'] = np.mean(von_mises)

            # Failure index
            if self.results.failure_index:
                failure = self.results.failure_index[elem_idx]
                elem_summary['max_failure_index'] = np.max(failure)
                elem_summary['failure_status'] = 'FAIL' if np.max(failure) > 1.0 else 'PASS'

            summary_data.append(elem_summary)

        # Create DataFrame and save
        df = pd.DataFrame(summary_data)
        summary_file = self.save_dir / "tertiary_results" / "tertiary_summary.csv"
        df.to_csv(summary_file, index=False, float_format='%.6e')

        logger.info(f"✅ Tertiary results summary saved: {summary_file}")

        # Log critical findings
        if self.results.failure_index:
            critical_elements = df[df['max_failure_index'] > 1.0]
            if not critical_elements.empty:
                logger.warning(
                    f"⚠️  FAILURE DETECTED in {len(critical_elements)} elements! "
                    f"See {summary_file} for details"
                )
            else:
                logger.info("✅ All elements passed failure criteria")

        return df

    def save_critical_locations(self):
        """
        Identify and save critical locations (highest stress/failure).

        Returns a CSV with the top N most critical Gauss points.
        """
        critical_data = []

        # Find critical Von Mises stress locations
        if self.results.von_mises_stress:
            for elem_idx, elem_von_mises in enumerate(self.results.von_mises_stress):
                for gp_idx, vm_stress in enumerate(elem_von_mises):
                    critical_data.append({
                        'element_id': elem_idx,
                        'gauss_point': gp_idx,
                        'von_mises_stress': vm_stress,
                        'failure_index': (
                            self.results.failure_index[elem_idx][gp_idx]
                            if self.results.failure_index else None
                        )
                    })

        # Create DataFrame and sort by stress
        df = pd.DataFrame(critical_data)
        df_sorted = df.sort_values('von_mises_stress', ascending=False)

        # Save top 100 critical locations
        critical_file = self.save_dir / "tertiary_results" / "critical_locations.csv"
        df_sorted.head(100).to_csv(critical_file, index=False, float_format='%.6e')

        logger.info(f"✅ Critical locations saved: {critical_file}")
        logger.info(f"   Top stress: {df_sorted.iloc[0]['von_mises_stress']:.6e} Pa")
        logger.info(f"   Location: Element {df_sorted.iloc[0]['element_id']}, "
                   f"GP {df_sorted.iloc[0]['gauss_point']}")

        return df_sorted


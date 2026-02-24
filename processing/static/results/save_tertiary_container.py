 # processing\static\results\save_tertiary_container.py

import numpy as np
import pandas as pd
from pathlib import Path
from typing import Optional
import logging

logger = logging.getLogger(__name__)


class SaveTertiaryResults:
    """
    Saves tertiary results (section forces, principal stresses, etc.) to CSV.

    Tertiary results are highly derived engineering quantities used for
    design verification and analysis:
    - Section force resultants [N, Vy, Vz, T, My, Mz]
    - Principal stresses [σ1, σ2, σ3]
    - Von Mises equivalent stress
    - Maximum shear stress

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

        # Create directory structure (resolution-first, consistent with primary/secondary)
        self.tertiary_dir = self.save_dir / "tertiary_results"
        self.gaussian_dir = self.tertiary_dir / "gaussian"
        self.section_forces_dir = self.gaussian_dir / "section_forces"
        self.principal_stress_dir = self.gaussian_dir / "principal_stress"
        self.elemental_dir = self.tertiary_dir / "elemental"
        self.nodal_dir = self.tertiary_dir / "nodal"

        for directory in [
            self.gaussian_dir,
            self.section_forces_dir,
            self.principal_stress_dir,
            self.elemental_dir,
            self.nodal_dir,
        ]:
            directory.mkdir(parents=True, exist_ok=True)

    def save_all(self):
        """
        Save tertiary results to disk.
        
        Saves only the core engineering resultants:
        - Section forces [N, Vy, Vz, T, My, Mz]
        - Principal stresses [σ1, σ2, σ3]
        
        Note: Von Mises and max shear are computed but not saved separately;
        they are included in the summary CSV.
        """
        logger.info("=" * 70)
        logger.info("SAVING TERTIARY RESULTS")
        logger.info("=" * 70)

        saved_count = 0

        # Save section forces (Gaussian resolution)
        if self.results.section_forces:
            logger.info("\n[1/3] Saving section force resultants (Gaussian resolution)...")
            self._save_section_forces()
            saved_count += 1

        # Save principal stresses (Gaussian resolution)
        if self.results.principal_stresses:
            logger.info("\n[2/3] Saving principal stresses (Gaussian resolution)...")
            self._save_principal_stresses()
            saved_count += 1

        # Save integrated elemental results
        if self.results.total_strain_energy or self.results.integrated_section_forces:
            logger.info("\n[3/3] Saving integrated elemental results...")
            self._save_integrated_elemental_results()
            saved_count += 1

        # Save nodal section forces (projected from Gaussian)
        if self.results.nodal_section_forces is not None:
            logger.info("\n[4/4] Saving nodal section forces...")
            self._save_nodal_section_forces()
            saved_count += 1

        logger.info("\n" + "=" * 70)
        logger.info(f"✅ TERTIARY RESULTS SAVED ({saved_count} categories)")
        logger.info(f"   Location: {self.tertiary_dir}")
        logger.info(f"   Note: Von Mises and shear stress in summary CSV")
        logger.info("=" * 70)

    # Marker so readers can detect column order; legacy CSVs (no marker) were in formulation order
    _SECTION_FORCE_FORMAT_LINE = "# column_order=resultant\n"

    def _save_section_forces(self):
        """Save section force resultants [N, Vy, Vz, T, My, Mz] per Gauss point.
        Rows are in ascending natural coordinate xi (same order as formulation cache).
        """
        for elem_idx, elem_section_forces in enumerate(self.results.section_forces):
            section_force_array = np.array(elem_section_forces)  # shape: (n_gauss, 6)
            n_gp = section_force_array.shape[0]
            filename = self.section_forces_dir / f"section_forces_elem_{elem_idx:06d}.csv"
            xi_per_row = np.polynomial.legendre.leggauss(n_gp)[0]  # ascending
            with open(filename, "w", encoding="utf-8") as f:
                f.write(self._SECTION_FORCE_FORMAT_LINE)
                f.write(f"# xi_per_row={','.join(f'{x:.16g}' for x in xi_per_row)}\n")
                np.savetxt(f, section_force_array, delimiter=",",
                           fmt="%.12e", header="N,Vy,Vz,T,My,Mz", comments="")

        logger.info(f"   ✓ Section forces saved: {len(self.results.section_forces)} elements")

    def _save_principal_stresses(self):
        """Save principal stresses [sigma1, sigma2, sigma3] per Gauss point."""
        for elem_idx, elem_principals in enumerate(self.results.principal_stresses):
            # Stack all Gauss point principal stresses
            principal_array = np.array(elem_principals)  # shape: (n_gauss, 3)
            filename = self.principal_stress_dir / f"principal_stress_elem_{elem_idx:06d}.csv"
            header = "sigma1,sigma2,sigma3"
            np.savetxt(filename, principal_array, delimiter=",",
                      fmt="%.12e", header=header, comments='')

        logger.info(f"   ✓ Principal stresses saved: {len(self.results.principal_stresses)} elements")

    def _save_integrated_elemental_results(self):
        """Save integrated elemental results (total strain energy and integrated section forces)."""
        # Save total strain energy per element
        if self.results.total_strain_energy is not None:
            data = {
                'element_id': range(len(self.results.total_strain_energy)),
                'total_strain_energy': self.results.total_strain_energy
            }
            df = pd.DataFrame(data)
            filename = self.elemental_dir / "total_strain_energy.csv"
            df.to_csv(filename, index=False, float_format='%.12e')
            logger.info(f"   ✓ Total strain energy saved: {len(self.results.total_strain_energy)} elements")

        # Save integrated section forces per element
        if self.results.integrated_section_forces is not None:
            integrated_forces_array = np.array(self.results.integrated_section_forces)  # shape: (n_elements, 6)
            filename = self.elemental_dir / "integrated_section_forces.csv"
            with open(filename, "w", encoding="utf-8") as f:
                f.write(self._SECTION_FORCE_FORMAT_LINE)
                np.savetxt(f, integrated_forces_array, delimiter=",",
                           fmt="%.12e", header="N,Vy,Vz,T,My,Mz", comments="")
            logger.info(f"   ✓ Integrated section forces saved: {len(self.results.integrated_section_forces)} elements")

    def _save_nodal_section_forces(self):
        """Save nodal section force resultants [N, Vy, Vz, T, My, Mz], one row per node."""
        arr = np.asarray(self.results.nodal_section_forces)
        if arr.ndim != 2 or arr.shape[1] != 6:
            logger.warning("nodal_section_forces shape %s unexpected, skipping save", arr.shape)
            return
        filename = self.nodal_dir / "nodal_section_forces.csv"
        with open(filename, "w", encoding="utf-8") as f:
            f.write(self._SECTION_FORCE_FORMAT_LINE)
            np.savetxt(f, arr, delimiter=",", fmt="%.12e", header="N,Vy,Vz,T,My,Mz", comments="")
        logger.info(f"   ✓ Nodal section forces saved: {len(arr)} nodes")


class SaveTertiaryResultsSummary:
    """
    Generates and saves summary statistics and critical locations for tertiary results.

    Creates comprehensive summary files including:
    - Maximum/minimum stress values
    - Critical element/Gauss point locations
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
        self.tertiary_dir = self.save_dir / "tertiary_results"
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
                elem_summary['max_principal_s1'] = np.max(principals[:, 0])
                elem_summary['min_principal_s3'] = np.min(principals[:, 2])

            # Von Mises stress
            if self.results.von_mises_stress:
                von_mises = self.results.von_mises_stress[elem_idx]
                elem_summary['max_von_mises'] = np.max(von_mises)
                elem_summary['mean_von_mises'] = np.mean(von_mises)

            summary_data.append(elem_summary)

        # Create DataFrame and save
        df = pd.DataFrame(summary_data)
        summary_file = self.tertiary_dir / "tertiary_summary.csv"
        df.to_csv(summary_file, index=False, float_format='%.6e')

        logger.info(f"✅ Tertiary results summary saved: {summary_file}")

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
                    })

        # Create DataFrame and sort by stress
        df = pd.DataFrame(critical_data)
        df_sorted = df.sort_values('von_mises_stress', ascending=False)

        # Save top 100 critical locations
        critical_file = self.tertiary_dir / "critical_locations.csv"
        df_sorted.head(100).to_csv(critical_file, index=False, float_format='%.6e')

        logger.info(f"✅ Critical locations saved: {critical_file}")
        logger.info(f"   Top stress: {df_sorted.iloc[0]['von_mises_stress']:.6e} Pa")
        logger.info(f"   Location: Element {df_sorted.iloc[0]['element_id']}, "
                   f"GP {df_sorted.iloc[0]['gauss_point']}")

        return df_sorted


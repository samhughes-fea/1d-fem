# processing_OOP\static\results\save_secondary_container.py

import numpy as np
import pandas as pd
from pathlib import Path
from typing import Optional
import logging

logger = logging.getLogger(__name__)

# Beam formulation outputs stress conjugate to ε = [ε_x, κ_y, κ_z, γ_xy, γ_xz, φ_x]:
# stress = [N, M_y, M_z, V_y, V_z, T]. Section force / CSV convention is [N, Vy, Vz, T, My, Mz].
# Bar has 2 components (e.g. N_axial, T), truss has 3 (N_axial, V_trans, T).
STRESS_HEADER_RESULTANTS = "N,Vy,Vz,T,My,Mz"
STRESS_HEADER_BAR = "N_axial,T"
STRESS_HEADER_TRUSS = "N_axial,V_trans,T"
STRAIN_HEADER_6 = "ε_xx,ε_yy,ε_zz,γ_xy,γ_yz,γ_xz"


def _stress_formulation_to_resultants(stress: np.ndarray) -> np.ndarray:
    """Reorder 6-component beam stress [N, M_y, M_z, V_y, V_z, T] to [N, Vy, Vz, T, My, Mz]. Only for shape (..., 6)."""
    if stress.ndim == 1 and stress.size != 6:
        return stress
    if stress.ndim == 2 and stress.shape[1] != 6:
        return stress
    REORDER = (0, 3, 4, 5, 1, 2)
    if stress.ndim == 1:
        return np.array([stress[i] for i in REORDER], dtype=stress.dtype)
    return stress[:, REORDER]


def _strain_header(n_components: int) -> str:
    """Header for strain CSV based on number of components (bar 2, truss 3, beam 6)."""
    if n_components == 6:
        return STRAIN_HEADER_6
    if n_components == 2:
        return "ε_axial,φ_torsion"
    if n_components == 3:
        return "ε_axial,γ_transverse,φ_torsion"
    return ",".join(f"comp{i}" for i in range(n_components))


def _stress_header(n_components: int) -> str:
    """Header for stress CSV (bar 2, truss 3, beam 6)."""
    if n_components == 6:
        return STRESS_HEADER_RESULTANTS
    if n_components == 2:
        return STRESS_HEADER_BAR
    if n_components == 3:
        return STRESS_HEADER_TRUSS
    return ",".join(f"comp{i}" for i in range(n_components))


class SaveSecondaryResults:
    """
    Saves secondary results (strain, stress, energy) to CSV files.

    Every run saves all resolutions when data is present:
    - Gaussian resolution: strain, stress, energy density at integration points
    - Nodal resolution: extrapolated/interpolated field values at nodes
    - Element resolution: integrated quantities (when populated)

    Parameters
    ----------
    secondary_results : SecondaryResultSet
        Container with all secondary results
    save_dir : str or Path
        Base directory where results should be saved
    """

    def __init__(
        self,
        secondary_results,
        save_dir: str | Path,
    ):
        self.results = secondary_results
        self.save_dir = Path(save_dir)

        # Create directory structure
        # Note: save_dir is already secondary_results_dir, so no need to nest
        self.secondary_dir = self.save_dir
        self.gaussian_dir = self.secondary_dir / "gaussian"
        self.nodal_dir = self.secondary_dir / "nodal"
        self.elemental_dir = self.secondary_dir / "elemental"

        for directory in [self.gaussian_dir, self.nodal_dir, self.elemental_dir]:
            directory.mkdir(parents=True, exist_ok=True)

    def save_all(self):
        """Save all secondary results to disk."""
        logger.info("=" * 70)
        logger.info("SAVING SECONDARY RESULTS")
        logger.info("=" * 70)

        saved_count = 0

        # Save Gaussian resolution results
        if self.results.gaussian_results:
            logger.info("\n[1/3] Saving Gaussian resolution results...")
            self._save_gaussian_results()
            saved_count += 1

        # Save nodal resolution results
        if self.results.nodal_results:
            logger.info("\n[2/3] Saving nodal resolution results...")
            self._save_nodal_results()
            saved_count += 1

        # Save elemental resolution results (if populated)
        # Note: Total strain energy and integrated resultants are produced by the
        # tertiary pipeline, not secondary; secondary_results.elemental_results is
        # typically None unless populated elsewhere.
        if self.results.elemental_results:
            logger.info("\n[3/3] Saving elemental resolution results...")
            self._save_elemental_results()
            saved_count += 1

        logger.info("\n" + "=" * 70)
        logger.info(f"✅ SECONDARY RESULTS SAVED ({saved_count} categories)")
        logger.info(f"   Location: {self.secondary_dir}")
        logger.info("=" * 70)

    def _save_gaussian_results(self):
        """
        Save Gauss point resolution results.
        
        Only pure field quantities (strain, stress, energy density) are saved here.
        Section forces [N, Vy, Vz, T, My, Mz] are integrated stress resultants
        and are saved via save_tertiary_container.py instead.
        """
        gauss_res = self.results.gaussian_results

        # Create subdirectories for field quantities only
        strain_dir = self.gaussian_dir / "strain"
        stress_dir = self.gaussian_dir / "stress"
        energy_dir = self.gaussian_dir / "energy_density"

        for d in [strain_dir, stress_dir, energy_dir]:
            d.mkdir(parents=True, exist_ok=True)

        # Save strain at Gauss points (variable components: bar 2, truss 3, beam 6)
        if gauss_res.strain:
            for elem_idx, elem_strains in enumerate(gauss_res.strain):
                strain_array = np.array(elem_strains)
                n_comp = strain_array.shape[1] if strain_array.ndim > 1 else strain_array.size
                if strain_array.ndim == 1:
                    strain_array = strain_array.reshape(1, -1)
                filename = strain_dir / f"strain_elem_{elem_idx:06d}.csv"
                header = _strain_header(n_comp)
                with open(filename, 'w', encoding='utf-8') as f:
                    np.savetxt(f, strain_array, delimiter=",", fmt="%.12e",
                              header=header, comments='')
            logger.info(f"   ✓ Strain saved: {len(gauss_res.strain)} elements")

        # Save stress at Gauss points (variable components: bar 2, truss 3, beam 6)
        if gauss_res.stress:
            for elem_idx, elem_stresses in enumerate(gauss_res.stress):
                stress_array = np.array(elem_stresses)
                if stress_array.ndim == 1:
                    stress_array = stress_array.reshape(1, -1)
                n_comp = stress_array.shape[1]
                if n_comp == 6:
                    stress_array = _stress_formulation_to_resultants(stress_array)
                header = _stress_header(n_comp)
                filename = stress_dir / f"stress_elem_{elem_idx:06d}.csv"
                with open(filename, 'w', encoding='utf-8') as f:
                    np.savetxt(f, stress_array, delimiter=",", fmt="%.12e",
                              header=header, comments='')
            logger.info(f"   ✓ Stress saved: {len(gauss_res.stress)} elements")

        # Save energy density at Gauss points
        if gauss_res.internal_energy_density:
            for elem_idx, elem_energy in enumerate(gauss_res.internal_energy_density):
                energy_array = np.array(elem_energy).reshape(-1, 1)
                filename = energy_dir / f"energy_density_elem_{elem_idx:06d}.csv"
                header = "strain_energy_density"
                with open(filename, 'w', encoding='utf-8') as f:
                    np.savetxt(f, energy_array, delimiter=",", fmt="%.12e",
                              header=header, comments='')
            logger.info(f"   ✓ Energy density saved: {len(gauss_res.internal_energy_density)} elements")

    def _save_nodal_results(self):
        """
        Save nodal resolution results.
        
        Only pure field quantities (strain, stress, energy density) are saved.
        Section forces are integrated resultants, not nodal fields.
        """
        nodal_res = self.results.nodal_results

        # Save nodal strain (variable components: bar 2, truss 3, beam 6)
        if nodal_res.strain is not None:
            filename = self.nodal_dir / "nodal_strain.csv"
            n_comp = nodal_res.strain.shape[1] if nodal_res.strain.ndim > 1 else nodal_res.strain.size
            header = _strain_header(n_comp)
            with open(filename, 'w', encoding='utf-8') as f:
                np.savetxt(f, nodal_res.strain, delimiter=",",
                          fmt="%.12e", header=header, comments='')
            logger.info(f"   ✓ Nodal strain saved: {nodal_res.strain.shape}")

        # Save nodal stress (variable components: bar 2, truss 3, beam 6)
        if nodal_res.stress is not None:
            filename = self.nodal_dir / "nodal_stress.csv"
            stress_out = nodal_res.stress
            n_comp = stress_out.shape[1] if stress_out.ndim > 1 else stress_out.size
            if n_comp == 6:
                stress_out = _stress_formulation_to_resultants(nodal_res.stress)
            header = _stress_header(n_comp)
            with open(filename, 'w', encoding='utf-8') as f:
                np.savetxt(f, stress_out, delimiter=",",
                          fmt="%.12e", header=header, comments='')
            logger.info(f"   ✓ Nodal stress saved: {nodal_res.stress.shape}")

        # Save nodal strain energy density
        if nodal_res.strain_energy_density is not None:
            filename = self.nodal_dir / "nodal_strain_energy_density.csv"
            with open(filename, 'w', encoding='utf-8') as f:
                np.savetxt(f, nodal_res.strain_energy_density.reshape(-1, 1),
                          delimiter=",", fmt="%.12e", header="strain_energy_density", comments='')
            logger.info(f"   ✓ Nodal energy density saved: {nodal_res.strain_energy_density.shape}")

    def _save_elemental_results(self):
        """Save element-level integrated results."""
        elem_res = self.results.elemental_results

        # Save total strain energy per element
        if elem_res.total_strain_energy is not None:
            data = {
                'element_id': range(len(elem_res.total_strain_energy)),
                'strain_energy': elem_res.total_strain_energy
            }
            df = pd.DataFrame(data)
            filename = self.elemental_dir / "element_strain_energy.csv"
            df.to_csv(filename, index=False, float_format='%.12e')
            logger.info(f"   ✓ Element strain energy saved: {len(elem_res.total_strain_energy)} elements")


class SaveSecondaryResultsSummary:
    """
    Generates and saves summary statistics for secondary results.

    Creates summary files with min/max/mean values for quick assessment.
    """

    def __init__(
        self,
        secondary_results,
        save_dir: str | Path
    ):
        self.results = secondary_results
        self.save_dir = Path(save_dir)

    def save(self):
        """Generate and save summary statistics."""
        summary_data = {}

        # Gaussian results summary
        if self.results.gaussian_results:
            summary_data.update(self._summarize_gaussian())

        # Nodal results summary
        if self.results.nodal_results:
            summary_data.update(self._summarize_nodal())

        # Save summary to CSV (save_dir is already secondary_results dir)
        summary_file = self.save_dir / "secondary_summary.csv"
        df = pd.DataFrame([summary_data])
        df.to_csv(summary_file, index=False, float_format='%.6e')

        logger.info(f"✅ Secondary results summary saved: {summary_file}")
        return df

    def _summarize_gaussian(self) -> dict:
        """Summarize Gaussian resolution results."""
        gauss_res = self.results.gaussian_results
        summary = {}

        if gauss_res.stress:
            all_stresses = np.concatenate([
                np.array(elem_stress) for elem_stress in gauss_res.stress
            ])
            summary['max_stress'] = np.max(np.abs(all_stresses))
            summary['min_stress'] = np.min(np.abs(all_stresses))
            summary['mean_stress'] = np.mean(np.abs(all_stresses))

        if gauss_res.strain:
            all_strains = np.concatenate([
                np.array(elem_strain) for elem_strain in gauss_res.strain
            ])
            summary['max_strain'] = np.max(np.abs(all_strains))
            summary['min_strain'] = np.min(np.abs(all_strains))
            summary['mean_strain'] = np.mean(np.abs(all_strains))

        return summary

    def _summarize_nodal(self) -> dict:
        """Summarize nodal resolution results."""
        nodal_res = self.results.nodal_results
        summary = {}

        if nodal_res.stress is not None:
            summary['max_nodal_stress'] = np.max(np.abs(nodal_res.stress))
            summary['min_nodal_stress'] = np.min(np.abs(nodal_res.stress))

        if nodal_res.strain is not None:
            summary['max_nodal_strain'] = np.max(np.abs(nodal_res.strain))
            summary['min_nodal_strain'] = np.min(np.abs(nodal_res.strain))

        return summary


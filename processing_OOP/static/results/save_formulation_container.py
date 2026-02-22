# processing_OOP\static\results\save_formulation_container.py

import numpy as np
import pandas as pd
from pathlib import Path
from typing import List, Optional
import logging
from scipy.sparse import issparse

logger = logging.getLogger(__name__)


class SaveFormulationData:
    """
    Saves element formulation data including Gauss-level matrices.

    This module persists the cached element formulation data (FormulationResultSet)
    that includes:
    - Element stiffness matrices K_e
    - Element force vectors F_e
    - Gauss point data (B, D, J, weights, coordinates)
    - Element type tracking

    This enables post-processing of stresses, strains, and section forces
    without recomputing shape functions.

    Parameters
    ----------
    formulation_cache : FormulationResultSet
        Cached element formulation data with ElementObject and ForceObject lists
    save_dir : str or Path
        Base directory where formulation data should be saved
    save_gauss_data : bool
        Whether to save detailed Gauss point data (can be large)
    """

    def __init__(
        self,
        formulation_cache = None,
        save_dir: str | Path = None,
        save_gauss_data: bool = False
    ):
        self.element_objects = formulation_cache.element_objects if formulation_cache else []
        self.force_objects = formulation_cache.force_objects if formulation_cache else []
        self.save_dir = Path(save_dir)
        self.save_gauss_data = save_gauss_data

        # Create directory structure
        self.formulation_dir = self.save_dir / "formulation"
        self.stiffness_dir = self.formulation_dir / "stiffness"
        self.force_dir = self.formulation_dir / "force"
        self.gauss_dir = self.formulation_dir / "gauss_points"
        self.element_type_file = self.formulation_dir / "element_types.csv"

        for directory in [self.stiffness_dir, self.force_dir]:
            directory.mkdir(parents=True, exist_ok=True)

        if self.save_gauss_data:
            self.gauss_dir.mkdir(parents=True, exist_ok=True)

    def save_all(self):
        """Save all formulation data to disk."""
        logger.info("=" * 70)
        logger.info("SAVING ELEMENT FORMULATION DATA")
        logger.info("=" * 70)

        # Save element type mapping
        if self.element_objects:
            logger.info(f"\n[1/5] Saving element type mapping...")
            self._save_element_types()

        # Save element stiffness matrices
        if self.element_objects:
            logger.info(f"\n[2/5] Saving {len(self.element_objects)} element stiffness matrices...")
            self._save_stiffness_matrices()

        # Save element force vectors
        if self.force_objects:
            logger.info(f"\n[3/5] Saving {len(self.force_objects)} element force vectors...")
            self._save_force_vectors()

        # Save Gauss point data (optional, can be large)
        if self.save_gauss_data and self.element_objects:
            logger.info(f"\n[4/5] Saving Gauss point formulation data...")
            self._save_gauss_point_data()
        else:
            logger.info(f"\n[4/5] Skipping Gauss point data (save_gauss_data=False)")

        # Save B2 shape function coefficients when present (one set per element)
        if self.element_objects:
            n_saved = self._save_shape_function_coefficients()
            if n_saved:
                logger.info(f"\n[5/5] Saved B2 shape function coefficients for {n_saved} elements")
            else:
                logger.info(f"\n[5/5] No B2 shape function coefficients to save")

        logger.info("\n" + "=" * 70)
        logger.info("✅ FORMULATION DATA SAVED SUCCESSFULLY")
        logger.info(f"   Location: {self.formulation_dir}")
        logger.info("=" * 70)

    def _save_element_types(self):
        """Save element type mapping to CSV."""
        type_data = []
        for elem_obj in self.element_objects:
            type_data.append({
                'element_id': elem_obj.element_id,
                'element_type': elem_obj.element_type
            })
        
        df = pd.DataFrame(type_data)
        df.to_csv(self.element_type_file, index=False)
        logger.info(f"   ✓ Saved to: {self.element_type_file}")

    def _save_stiffness_matrices(self):
        """Save element stiffness matrices to CSV files."""
        for elem_obj in self.element_objects:
            elem_id = elem_obj.element_id
            K_e = elem_obj.K_e

            # Convert sparse to dense if needed
            if issparse(K_e):
                K_e = K_e.toarray()

            # Save to CSV
            filename = self.stiffness_dir / f"K_e_{elem_id:06d}.csv"
            np.savetxt(filename, K_e, delimiter=",", fmt="%.12e")

        logger.info(f"   ✓ Saved to: {self.stiffness_dir}")

    def _save_force_vectors(self):
        """Save element force vectors to CSV files."""
        for force_obj in self.force_objects:
            elem_id = force_obj.element_id
            F_e = force_obj.F_e

            # Save to CSV (as column vector)
            filename = self.force_dir / f"F_e_{elem_id:06d}.csv"
            np.savetxt(filename, F_e.reshape(-1, 1), delimiter=",", fmt="%.12e")

        logger.info(f"   ✓ Saved to: {self.force_dir}")

    def _save_gauss_point_data(self):
        """Save detailed Gauss point formulation data."""
        # Create subdirectories
        B_matrix_dir = self.gauss_dir / "B_matrices"
        D_matrix_dir = self.gauss_dir / "D_matrices"
        jacobian_dir = self.gauss_dir / "jacobians"
        coords_dir = self.gauss_dir / "coordinates"

        for d in [B_matrix_dir, D_matrix_dir, jacobian_dir, coords_dir]:
            d.mkdir(parents=True, exist_ok=True)

        for elem_obj in self.element_objects:
            elem_id = elem_obj.element_id

            # Save data for each Gauss point
            for gp_idx, gauss_data in enumerate(elem_obj.gauss_data):
                # B matrix (strain-displacement)
                B_file = B_matrix_dir / f"B_elem{elem_id:06d}_gp{gp_idx:02d}.csv"
                np.savetxt(B_file, gauss_data.B_matrix, delimiter=",", fmt="%.12e")

                # D matrix (constitutive)
                D_file = D_matrix_dir / f"D_elem{elem_id:06d}_gp{gp_idx:02d}.csv"
                np.savetxt(D_file, gauss_data.D_matrix, delimiter=",", fmt="%.12e")

                # Jacobian and weight
                jacob_file = jacobian_dir / f"J_elem{elem_id:06d}_gp{gp_idx:02d}.csv"
                data = np.array([
                    [gauss_data.xi, gauss_data.weight, gauss_data.jacobian]
                ])
                np.savetxt(jacob_file, data, delimiter=",", fmt="%.12e",
                          header="xi,weight,jacobian", comments='')

        logger.info(f"   ✓ Saved to: {self.gauss_dir}")

    def _save_shape_function_coefficients(self) -> int:
        """Save B2 shape function coefficient arrays (12, 6, 4) per element when present. Returns count saved."""
        coeff_dir = self.formulation_dir / "shape_function_coefficients"
        saved = 0
        for elem_obj in self.element_objects:
            n_c = getattr(elem_obj, "shape_function_N_coefficients", None)
            dn_c = getattr(elem_obj, "shape_function_dN_dxi_coefficients", None)
            d2n_c = getattr(elem_obj, "shape_function_d2N_dxi2_coefficients", None)
            if n_c is None or dn_c is None or d2n_c is None:
                continue
            coeff_dir.mkdir(parents=True, exist_ok=True)
            elem_id = elem_obj.element_id
            # Save as (72, 4) CSV: rows (dof*6+comp), cols ξ^0..ξ^3
            for name, arr in [
                ("N", n_c),
                ("dN_dxi", dn_c),
                ("d2N_dxi2", d2n_c),
            ]:
                flat = arr.reshape(-1, 4)
                path = coeff_dir / f"{name}_elem{elem_id:06d}.csv"
                np.savetxt(path, flat, delimiter=",", fmt="%.12e")
            saved += 1
        return saved


class SaveFormulationSummary:
    """
    Saves summary information about element formulations.

    Creates a summary CSV file with key formulation metadata:
    - Element ID
    - Element Type
    - Number of Gauss points
    - Integration scheme
    - Stiffness matrix condition number
    - Stiffness matrix norm
    """

    def __init__(
        self,
        formulation_cache,
        save_dir: str | Path
    ):
        self.element_objects = formulation_cache.element_objects if formulation_cache else []
        self.save_dir = Path(save_dir)

    def save(self):
        """Generate and save formulation summary."""
        summary_data = []

        for elem_obj in self.element_objects:
            K_e = elem_obj.K_e
            if issparse(K_e):
                K_e = K_e.toarray()

            # Compute summary statistics
            cond_number = np.linalg.cond(K_e)
            frobenius_norm = np.linalg.norm(K_e, 'fro')
            
            summary_data.append({
                'element_id': elem_obj.element_id,
                'element_type': elem_obj.element_type,
                'n_gauss_points': elem_obj.n_gauss_points,
                'integration_scheme': elem_obj.integration_scheme,
                'K_condition_number': cond_number,
                'K_frobenius_norm': frobenius_norm,
                'K_max_entry': np.max(np.abs(K_e)),
                'K_min_entry': np.min(np.abs(K_e[K_e != 0]))  # exclude zeros
            })

        # Create DataFrame and save
        df = pd.DataFrame(summary_data)
        summary_file = self.save_dir / "formulation" / "formulation_summary.csv"
        df.to_csv(summary_file, index=False, float_format='%.6e')

        logger.info(f"✅ Formulation summary saved: {summary_file}")
        return df


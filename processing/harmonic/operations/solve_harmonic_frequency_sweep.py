# processing/harmonic/operations/solve_harmonic_frequency_sweep.py
"""Modal superposition or direct frequency sweep for harmonic response."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable, Optional, Union

import numpy as np
from scipy.sparse import csr_matrix

from processing.harmonic.damping_table import interpolate_zeta_hz
from processing.harmonic.eigen_basis_io import load_modal_basis_from_modal_results_dir
from processing.harmonic.frequency_response import (
    geometric_mean_reference_omega_rad,
    sweep_displacements,
)
from processing.harmonic.modal_superposition import (
    harmonic_displacement_modal_superposition,
    undamped_natural_modes,
)

from processing.harmonic.operations._logging import init_stage_logger


@dataclass
class HarmonicSweepConfig:
    use_modal: bool
    num_modal: int
    basis_dir: Optional[str]
    basis_job: Optional[str]
    zeta: float
    use_parallel: bool
    mp_ref: str
    rayleigh_alpha: float
    rayleigh_beta: float
    lin_solver: str
    zpf: Optional[np.ndarray]
    ft_tab: Any
    zt_tab: Any


class SolveHarmonicFrequencySweep:
    def __init__(
        self,
        job_results_dir: Optional[Union[str, Path]] = None,
        resolve_job_path: Optional[Callable[[str], str]] = None,
    ):
        self.job_results_dir = Path(job_results_dir) if job_results_dir else None
        self.resolve_job_path = resolve_job_path
        self._log = init_stage_logger(
            "SolveHarmonicFrequencySweep", self.job_results_dir
        )

    def run(
        self,
        K_mod: csr_matrix,
        M_mod: csr_matrix,
        C_mod: csr_matrix,
        F_use: np.ndarray,
        freqs_hz: np.ndarray,
        u_template: Optional[np.ndarray],
        dof_unknown: Optional[np.ndarray],
        cfg: HarmonicSweepConfig,
    ) -> np.ndarray:
        if cfg.use_modal:
            if dof_unknown is not None:
                raise ValueError(
                    "harmonic: use_modal_superposition is incompatible with non-zero prescribed "
                    "displacements in this version; use the direct frequency sweep instead."
                )
            if cfg.rayleigh_alpha != 0.0 or cfg.rayleigh_beta != 0.0:
                self._log.warning(
                    "harmonic: Rayleigh damping (rayleigh_alpha / rayleigh_beta) is not applied on "
                    "the modal superposition path; use the direct frequency sweep for full Rayleigh C."
                )
            self._log.info("Harmonic solve: modal superposition, num_modes=%s", cfg.num_modal)
            if cfg.basis_dir:
                if self.resolve_job_path is None:
                    raise ValueError("harmonic_modal_basis_dir requires resolve_job_path")
                fh_hz, Phi_full = load_modal_basis_from_modal_results_dir(
                    self.resolve_job_path(cfg.basis_dir),
                    job_name=cfg.basis_job,
                )
                nm = min(cfg.num_modal, int(Phi_full.shape[1]))
                Phi = np.asarray(Phi_full[:, :nm], dtype=np.float64)
                omega_n = (2.0 * np.pi * np.asarray(fh_hz[:nm], dtype=np.float64)).ravel()
            else:
                omega_n, Phi = undamped_natural_modes(K_mod, M_mod, cfg.num_modal)
            zeta_pm = None
            if cfg.ft_tab is not None:
                zeta_pm = interpolate_zeta_hz(omega_n / (2.0 * np.pi), cfg.ft_tab, cfg.zt_tab)
            return harmonic_displacement_modal_superposition(
                omega_n,
                Phi,
                F_use,
                freqs_hz,
                cfg.zeta,
                zeta_per_mode=zeta_pm,
            )

        self._log.info("Harmonic solve: direct frequency sweep")
        sweep_kw: dict = dict(
            parallel=cfg.use_parallel,
            mp_damping_reference=cfg.mp_ref,
            zeta_mp=cfg.zeta,
            rayleigh_alpha_mp=cfg.rayleigh_alpha,
            rayleigh_beta_mp=cfg.rayleigh_beta,
            linear_solver=cfg.lin_solver,
        )
        if cfg.zpf is not None:
            sweep_kw["zeta_per_frequency"] = cfg.zpf
            sweep_kw["geometric_omega_ref_rad"] = geometric_mean_reference_omega_rad(freqs_hz)
        return sweep_displacements(
            K_mod,
            M_mod,
            C_mod,
            F_use,
            freqs_hz,
            u_prescribed=u_template,
            dof_unknown=dof_unknown,
            **sweep_kw,
        )

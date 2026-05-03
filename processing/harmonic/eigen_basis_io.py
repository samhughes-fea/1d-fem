# processing/harmonic/eigen_basis_io.py
"""Load undamped mode shapes / frequencies saved by §2 eigen jobs."""

from __future__ import annotations

import glob
import os

import numpy as np


def load_modal_basis_from_modal_results_dir(
    modal_results_dir: str,
    *,
    job_name: str | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Load natural frequencies (Hz) and mode-shape matrix from ``modal_results`` folder.

    Expected files: ``{job_name}_frequencies.txt`` and ``{job_name}_mode_shapes.txt``
    (same layout as :meth:`VibrationBucklingBackend._save_primary_results`).

    If *job_name* is omitted, glob ``*_frequencies.txt`` and require exactly one match.
    """
    if not os.path.isdir(modal_results_dir):
        raise NotADirectoryError(f"modal basis directory not found: {modal_results_dir}")

    if job_name is not None:
        base = os.path.join(modal_results_dir, f"{job_name}_frequencies.txt")
        if not os.path.isfile(base):
            raise FileNotFoundError(f"eigen frequencies file not found: {base}")
    else:
        pat = os.path.join(modal_results_dir, "*_frequencies.txt")
        matches = sorted(glob.glob(pat))
        if len(matches) != 1:
            raise FileNotFoundError(
                f"expected exactly one *_frequencies.txt under {modal_results_dir}, found {len(matches)}"
            )
        base = matches[0]
        job_name = os.path.basename(base).replace("_frequencies.txt", "")

    freq_hz = np.loadtxt(base, dtype=np.float64)
    if freq_hz.ndim == 0:
        freq_hz = np.array([float(freq_hz)], dtype=np.float64)
    else:
        freq_hz = np.asarray(freq_hz, dtype=np.float64).ravel()

    ms_path = os.path.join(modal_results_dir, f"{job_name}_mode_shapes.txt")
    if not os.path.isfile(ms_path):
        raise FileNotFoundError(f"mode shapes file not found: {ms_path}")
    Phi = np.loadtxt(ms_path, dtype=np.float64)
    if Phi.ndim == 1:
        Phi = Phi.reshape(-1, 1)
    nm_file = int(freq_hz.size)
    if Phi.shape[1] != nm_file:
        raise ValueError(
            f"mode shape columns ({Phi.shape[1]}) != frequency count ({nm_file})"
        )
    return freq_hz, Phi

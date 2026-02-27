# post_processing/validation_visualisers/abaqus/config.py
"""
Paths and element-type mapping for Abaqus validation.
Used by job_to_abaqus_script.py and run_abaqus_cae.py.
"""
import os
import re
from pathlib import Path

# Resolve validation_visualisers dir (parent of abaqus/)
_THIS_DIR = Path(__file__).resolve().parent
VALIDATION_DIR = _THIS_DIR.parent
# Project root: parent of post_processing
PROJECT_ROOT = VALIDATION_DIR.parent.parent

JOBS_DIR = PROJECT_ROOT / "jobs"
FEM_RESULTS_DIR = PROJECT_ROOT / "post_processing" / "results"
ABAQUS_GENERATED_DIR = _THIS_DIR / "generated"
ABAQUS_RESULTS_DIR = VALIDATION_DIR / "abaqus_results"

# Match timestamped FEM result dir: job_0000_n8_2026-02-22_..._pid123_abc
# Used by discovery (run_all_abaqus_jobs --from-results) and comparison scripts.
RESULT_DIR_PATTERN = re.compile(
    r"job_(?P<base_id>\d+)_n(?P<n>\d+)_[\d\-_]+_pid\d+_[a-f0-9]+"
)

# Abaqus CAE installation root (e.g. C:\SIMULIA\CAE). If set, run_abaqus_cae uses this to find abaqus.bat.
# Override with env ABAQUS_CAE_ROOT. If unset, ABAQUS_CAE_CMD falls back to "abaqus" on PATH.
ABAQUS_CAE_ROOT = os.environ.get("ABAQUS_CAE_ROOT", r"C:\SIMULIA\CAE")


def _resolve_abaqus_cmd() -> str:
    """Return full path to abaqus command if ABAQUS_CAE_ROOT is set and a launcher exists; else 'abaqus'."""
    root = Path(ABAQUS_CAE_ROOT) if ABAQUS_CAE_ROOT else None
    if not root or not root.is_dir():
        return "abaqus"
    # Common launcher locations under SIMULIA CAE
    candidates = [
        root / "Commands" / "abaqus.bat",
        root / "Commands" / "abaqus.exe",
    ]
    for sub in root.iterdir():
        if sub.is_dir():
            candidates.append(sub / "win_b64" / "code" / "bin" / "abaqus.bat")
            candidates.append(sub / "win_b64" / "code" / "bin" / "abaqus.exe")
    for p in candidates:
        if p.is_file():
            return str(p)
    return "abaqus"


# Abaqus command: from ABAQUS_CAE_ROOT if available, else "abaqus" (must be on PATH)
ABAQUS_CAE_CMD = _resolve_abaqus_cmd()

# Full path for abqpy's ABAQUS_BAT_PATH env var (so saveAs() launches correct Abaqus). None = use "abaqus" on PATH.
_launcher_path = Path(ABAQUS_CAE_CMD)
ABAQUS_LAUNCHER_PATH = str(ABAQUS_CAE_CMD) if (_launcher_path.is_absolute() and _launcher_path.is_file()) else None

# Map our element type names to Abaqus beam element type
# B33 = 2-node cubic beam (Euler-Bernoulli, no shear deformation)
# B31 = 2-node linear beam (Timoshenko, shear deformation)
ELEMENT_TYPE_MAP = {
    "EulerBernoulliBeamElement3D": "B33",
    "TimoshenkoBeamElement3D": "B31",
}
# Unsupported for this validation (e.g. Levinson) are ignored or raise
SUPPORTED_ELEMENT_TYPES = set(ELEMENT_TYPE_MAP.keys())

# Abaqus mesh size used as converged reference (source of truth) for validation.
# Comparison scripts use Abaqus at n=500 as the benchmark; FEM n128 (and GCI FEM n32,64,128) are compared to it.
# Job dirs job_XXXX_n500 must exist (create via mesh variant scripts including n500) and be run in Abaqus.
# To run the n500 batch: run_all_abaqus_jobs.py --n500-only (see README_VALIDATION_ABAQUS.md).
ABAQUS_REFERENCE_N: int = 500

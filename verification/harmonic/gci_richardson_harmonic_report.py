"""
GCI and Richardson extrapolation report for harmonic FRF peak amplitude.

Uses three mesh levels n=32, 64, 128 (refinement ratio R=2, safety factor F_s=1.25).
The QoI is peak amplitude (max |H(ω)|) over the frequency grid. Richardson extrapolation
provides the "exact" estimate; no external analytical reference is used because the
SDOF analytical model uses different parameters than the physical beam FEM.

Outputs:
  - harmonic_plots/gci_tables/gci_richardson_harmonic.csv
  Console summary.
"""
from __future__ import annotations

import csv
import glob
import math
import re
import sys
from pathlib import Path
from typing import Final

SCRIPT_DIR: Final[Path] = Path(__file__).resolve().parent
PROJECT_ROOT: Final[Path] = SCRIPT_DIR.parents[1]
sys.path.insert(0, str(PROJECT_ROOT))

_RESULTS_DIR: Final[Path] = PROJECT_ROOT / "results"
_N_RE: Final[re.Pattern] = re.compile(r"benchmark_harmonic_sdof_n(\d+)")

N_FINE, N_MED, N_COARSE = 128, 64, 32
R: Final[float] = 2.0
F_S: Final[float] = 1.25


def _read_peak_amplitude(csv_path: Path) -> float | None:
    try:
        max_amp = 0.0
        with open(csv_path, "r", encoding="utf-8", newline="") as f:
            for row in csv.DictReader(f):
                real = float(row["response_real"])
                imag = float(row["response_imag"])
                amp = math.sqrt(real * real + imag * imag)
                if amp > max_amp:
                    max_amp = amp
        return max_amp if max_amp > 0.0 else None
    except Exception:
        return None


def _get_amp_for_n(n: int, found: list[str]) -> float | None:
    for csv_path_str in found:
        csv_path = Path(csv_path_str)
        m = _N_RE.search(csv_path.parents[2].name)
        if m and int(m.group(1)) == n:
            return _read_peak_amplitude(csv_path)
    return None


def _richardson_gci(
    phi_fine: float, phi_med: float, phi_coarse: float
) -> tuple[float, float, float, float]:
    e21 = phi_med - phi_fine
    e32 = phi_coarse - phi_med
    if abs(e21) < 1e-30 and abs(e32) < 1e-30:
        return float("nan"), phi_fine, 0.0, 0.0
    ratio = e32 / e21 if abs(e21) >= 1e-30 else float("nan")
    if ratio <= 0 or not math.isfinite(ratio):
        return float("nan"), phi_fine, float("nan"), float("nan")
    p_obs = math.log(abs(ratio)) / math.log(R)
    denom = (R ** p_obs) - 1.0
    if denom <= 0 or not math.isfinite(denom):
        return float("nan"), phi_fine, float("nan"), float("nan")
    phi_ext = phi_fine + (phi_fine - phi_med) / denom
    gci_12 = 100.0 * F_S * abs(e21 / phi_fine) / denom if abs(phi_fine) >= 1e-30 else float("nan")
    gci_23 = 100.0 * F_S * abs(e32 / phi_med) / denom if abs(phi_med) >= 1e-30 else float("nan")
    return p_obs, phi_ext, gci_12, gci_23


def _run_report() -> None:
    out_dir = SCRIPT_DIR / "harmonic_plots" / "gci_tables"
    out_dir.mkdir(parents=True, exist_ok=True)

    pattern = str(
        _RESULTS_DIR
        / "*benchmark_harmonic_sdof_n*"
        / "primary_results"
        / "harmonic"
        / "frequency_response.csv"
    )
    found = sorted(glob.glob(pattern))

    missing = [n for n in (N_FINE, N_MED, N_COARSE) if _get_amp_for_n(n, found) is None]
    if missing:
        print(
            f"GCI/Richardson skipped: missing results for n={missing}.\n"
            f"Run job_benchmark_harmonic_sdof_n{missing[0]} (and others) first."
        )
        return

    amp_fine = _get_amp_for_n(N_FINE, found)
    amp_med = _get_amp_for_n(N_MED, found)
    amp_coarse = _get_amp_for_n(N_COARSE, found)

    p_obs, phi_ext, gci_12, gci_23 = _richardson_gci(amp_fine, amp_med, amp_coarse)

    r_p = (R ** p_obs) if math.isfinite(p_obs) else float("nan")
    gci_ratio = (
        gci_23 / gci_12
        if math.isfinite(gci_12) and abs(gci_12) > 1e-30 and math.isfinite(gci_23)
        else float("nan")
    )
    fine_vs_ext_pct = (
        100.0 * abs(amp_fine - phi_ext) / abs(phi_ext)
        if math.isfinite(phi_ext) and abs(phi_ext) >= 1e-30
        else float("nan")
    )

    csv_path = out_dir / "gci_richardson_harmonic.csv"
    with open(csv_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow([
            "QoI", "n_fine", "n_med", "n_coarse",
            "phi_fine", "phi_med", "phi_coarse",
            "p_obs", "r_p", "phi_ext",
            "fine_vs_ext_pct", "GCI_12_pct", "GCI_23_pct", "GCI_23_12_ratio",
        ])
        writer.writerow([
            "peak_amplitude_m_per_n", N_FINE, N_MED, N_COARSE,
            f"{amp_fine:.12e}", f"{amp_med:.12e}", f"{amp_coarse:.12e}",
            f"{p_obs:.6f}", f"{r_p:.6f}", f"{phi_ext:.12e}",
            f"{fine_vs_ext_pct:.6f}" if math.isfinite(fine_vs_ext_pct) else "nan",
            f"{gci_12:.6f}" if math.isfinite(gci_12) else "nan",
            f"{gci_23:.6f}" if math.isfinite(gci_23) else "nan",
            f"{gci_ratio:.6f}" if math.isfinite(gci_ratio) else "nan",
        ])
    print(f"Saved: {csv_path}")

    def _fmt(v: float, d: int = 4) -> str:
        return f"{v:.{d}f}" if math.isfinite(v) else "N/A"

    print(f"\n--- GCI / Richardson: harmonic peak amplitude ---")
    print(f"  Grids: n={N_COARSE}, {N_MED}, {N_FINE}   R={R}   F_s={F_S}")
    print(f"  phi_coarse (n={N_COARSE}) = {amp_coarse:.6e} m/N")
    print(f"  phi_med    (n={N_MED})  = {amp_med:.6e} m/N")
    print(f"  phi_fine   (n={N_FINE}) = {amp_fine:.6e} m/N")
    print(f"  phi_ext    (Richardson) = {phi_ext:.6e} m/N")
    print(f"  p_obs   = {_fmt(p_obs)}  (expected: 4 for cubic EB element)")
    print(f"  r^p     = {_fmt(r_p)}")
    print(f"  GCI_12  = {_fmt(gci_12)}%  (fine-grid GCI)")
    print(f"  GCI_23  = {_fmt(gci_23)}%  (medium-grid GCI)")
    print(f"  GCI_23/GCI_12 = {_fmt(gci_ratio)}  (expect ≈ r^p in asymptotic range)")
    print(f"  |phi_fine - phi_ext| / phi_ext = {_fmt(fine_vs_ext_pct)}%")
    print("")


if __name__ == "__main__":
    _run_report()
    sys.exit(0)

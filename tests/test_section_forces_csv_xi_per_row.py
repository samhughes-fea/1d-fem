"""Tests for section-forces element CSV read (xi_per_row vs Gauss–Legendre fallback)."""

from __future__ import annotations

import sys
import tempfile
from pathlib import Path

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))


def test_read_elem_section_forces_csv_uses_xi_per_row_when_present():
    from post_processing.graphical_visualisers.tertiary_visualisers.section_forces.section_forces_visualisation import (
        _read_elem_section_forces_csv,
    )

    xi_saved = np.array([-0.5, 0.0, 0.25], dtype=np.float64)
    rows = np.ones((3, 6)) * 100.0
    with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False, encoding="utf-8") as f:
        p = Path(f.name)
        f.write("# column_order=resultant\n")
        f.write("# xi_per_row=" + ",".join(str(x) for x in xi_saved) + "\n")
        f.write("N,Vy,Vz,T,My,Mz\n")
        for r in rows:
            f.write(",".join(str(x) for x in r) + "\n")
    try:
        data, xi_used = _read_elem_section_forces_csv(p)
        assert data.shape == (3, 6)
        np.testing.assert_allclose(xi_used, xi_saved)
    finally:
        p.unlink(missing_ok=True)


def test_read_elem_section_forces_csv_fallback_leggauss_without_xi_line():
    from post_processing.graphical_visualisers.tertiary_visualisers.section_forces.section_forces_visualisation import (
        _read_elem_section_forces_csv,
    )

    rows = np.ones((4, 6))
    with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False, encoding="utf-8") as f:
        p = Path(f.name)
        f.write("# column_order=resultant\n")
        f.write("N,Vy,Vz,T,My,Mz\n")
        for r in rows:
            f.write(",".join(str(x) for x in r) + "\n")
    try:
        _, xi_used = _read_elem_section_forces_csv(p)
        xi_exp = np.polynomial.legendre.leggauss(4)[0]
        np.testing.assert_allclose(xi_used, xi_exp)
    finally:
        p.unlink(missing_ok=True)

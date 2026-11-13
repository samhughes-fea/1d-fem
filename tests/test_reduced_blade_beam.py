"""A blade beam job can be built from stations and a reduced section."""

import ast
import math
from pathlib import Path

import numpy as np

from onedfem.pre_processing.mesh_library.reduced_blade_beam import reduced_blade_beam_job

SRC = Path(__file__).resolve().parents[1] / "src" / "onedfem" / "pre_processing" / "mesh_library" / "reduced_blade_beam.py"


class _Station:
    def __init__(self, z, twist, offset=(0.0, 0.0)):
        self.z = z
        self.twist = twist
        self.offset = offset


class _Geometry:
    """Stations only. Asking for a cut is a failure."""

    def __init__(self):
        self.stations = [_Station(0.0, 0.0, (0.1, -0.2)), _Station(4.0, 10.0, (0.1, -0.2))]

    def section_at(self, z):
        raise AssertionError(f"the beam job asked for a wall cut at {z}")


def _section(z):
    return {
        "area": 0.02 + 0.001 * z,
        "ixx": 4.0e-4,
        "iyy": 1.0e-4,
        "shear_centre": (0.05, -0.01),
        "warping_constant": 0.0,
        "torsion_constant": 0.0,
    }


def test_the_writer_never_names_a_wall_graph():
    tree = ast.parse(SRC.read_text())
    named = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Attribute):
            named.add(node.attr)
        elif isinstance(node, ast.Name):
            named.add(node.id)
    assert {"section_at", "walls", "shell_graph", "mid_surface"}.isdisjoint(named)


def test_stations_and_a_reduced_section_write_the_beam_job(tmp_path):
    job = tmp_path / "job"
    out = reduced_blade_beam_job(str(job), _Geometry(), _section, n_elements=2)
    assert not (job / "section_geometry.txt").exists()
    assert not (job / "laminate.txt").exists()
    assert not (job / "element_section.txt").exists()
    grid = (job / "grid.txt").read_text().splitlines()
    assert grid[2].split()[1:] == ["0.0", "0.1", "-0.2"]
    row = (job / "section.txt").read_text().splitlines()[2].split()
    # element 0, mid-span z = 1: A, I_x, I_y, I_z, J, kappa, alpha, y_sc, z_sc, Gamma
    assert math.isclose(float(row[1]), 0.021)
    assert float(row[3]) == 4.0e-4
    assert float(row[4]) == 1.0e-4
    assert float(row[8]) == 0.05
    assert float(row[9]) == -0.01
    assert float(row[10]) == 0.0
    length = 2.0
    assert np.allclose(out["k_x0"], -math.radians(5.0) / length)

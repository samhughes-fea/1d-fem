"""
Unit tests for the Levinson beam D-matrix (MaterialStiffnessOperator).

Ensures shear stiffness is GA (no shear correction factor κ), consistent with
higher-order kinematics in the B-matrix.
"""

import sys
from pathlib import Path

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from pre_processing.element_library.levinson.utilities.D_matrix import MaterialStiffnessOperator


def test_levinson_D_matrix_shear_terms_are_GA_no_kappa():
    """Levinson D-matrix shear diagonal is G*A; no κ (unlike Timoshenko)."""
    E = 2.1e11
    G = 8.1e10
    A = 0.00131
    I_y = 1.0e-6
    I_z = 2.08769e-06
    J_t = 1.0e-8

    op = MaterialStiffnessOperator(
        youngs_modulus=E,
        shear_modulus=G,
        cross_section_area=A,
        moment_inertia_y=I_y,
        moment_inertia_z=I_z,
        torsion_constant=J_t,
    )

    D = op.assembly_form()
    GA_expected = G * A

    assert D.shape == (6, 6), "D must be 6×6"
    assert D[3, 3] == GA_expected, "D[3,3] (shear xy) must be G*A, not κGA"
    assert D[4, 4] == GA_expected, "D[4,4] (shear xz) must be G*A, not κGA"

    # Shear diagonal only; no off-diagonal coupling into shear rows/cols (for standard case)
    assert D[3, 4] == 0.0 and D[4, 3] == 0.0, "Shear block must be diagonal"


def test_levinson_MaterialStiffnessOperator_has_no_shear_correction_factor():
    """Levinson formulation must not use κ; operator has no shear_correction_factor."""
    E = 2.1e11
    G = 8.1e10
    A = 0.00131
    I_y = 1.0e-6
    I_z = 2.08769e-06
    J_t = 1.0e-8

    op = MaterialStiffnessOperator(
        youngs_modulus=E,
        shear_modulus=G,
        cross_section_area=A,
        moment_inertia_y=I_y,
        moment_inertia_z=I_z,
        torsion_constant=J_t,
    )

    assert not hasattr(op, "shear_correction_factor"), (
        "Levinson MaterialStiffnessOperator must not have shear_correction_factor; "
        "shear is accounted for by higher-order kinematics (GA, no κ)."
    )

# pre_processing/element_library/nonlinear/large_rotations/gebt_unshearable/gebt_unshearable_3D.py
"""
Unshearable (Euler–Bernoulli) beam with finite-strain TL stack (registration ``GEBTUnshearableBeamElement3D``).

**Current implementation:** identical weak form and tangent to :class:`NonlinearEulerBernoulliBeamElement3D`
(Total Lagrangian). The ``GEBTUnshearable-3D`` label distinguishes logging from the generic nonlinear EB type.

**Roadmap:** specialized selective-integration layout or classical unshearable GEBT extras if required.

See ``docs/element_library/large_rotation_vs_total_lagrangian.md``.
"""

from __future__ import annotations

from pre_processing.element_library.nonlinear.euler_bernoulli.nonlinear_euler_bernoulli_3D import (
    NonlinearEulerBernoulliBeamElement3D,
)


class GEBTUnshearableBeamElement3D(NonlinearEulerBernoulliBeamElement3D):
    """Unshearable large-rotation beam using the nonlinear EB TL pipeline."""

    element_type_name = "GEBTUnshearable-3D"

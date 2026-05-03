"""
Removed or deregistered public element type strings.

Used by :class:`~pre_processing.element_library.element_factory.ElementFactory`
before resolving ``ELEMENT_CLASS_MAP`` so users get a clear migration error.

Emits :exc:`DeprecationWarning` immediately before raising :exc:`ValueError` for observability
in environments that surface deprecation warnings.
"""

from __future__ import annotations

import warnings
from typing import Dict

# Values are full error messages (raised as ValueError).
REMOVED_ELEMENT_TYPES: Dict[str, str] = {
    "LinearWarpingEulerBernoulliBeamElement3D": (
        "Element type 'LinearWarpingEulerBernoulliBeamElement3D' was removed. "
        "Use 'LinearEulerBernoulliBeamElement3D' with [warping]=1 in element.txt "
        "(or set element_dictionary['warping']). See docs/conventions/DEPRECATED_ELEMENT_TYPES.md."
    ),
    "LinearWarpingTimoshenkoBeamElement3D": (
        "Element type 'LinearWarpingTimoshenkoBeamElement3D' was removed. "
        "Use 'LinearTimoshenkoBeamElement3D' with [warping]=1 in element.txt "
        "(or set element_dictionary['warping']). See docs/conventions/DEPRECATED_ELEMENT_TYPES.md."
    ),
    "LinearWarpingLevinsonBeamElement3D": (
        "Element type 'LinearWarpingLevinsonBeamElement3D' was removed. "
        "Use 'LinearLevinsonBeamElement3D' with [warping]=1 in element.txt "
        "(or set element_dictionary['warping']). See docs/conventions/DEPRECATED_ELEMENT_TYPES.md."
    ),
    "LinearWarpingReddyBeamElement3D": (
        "Element type 'LinearWarpingReddyBeamElement3D' was removed. "
        "Use 'LinearReddyBeamElement3D' with [warping]=1 in element.txt "
        "(or set element_dictionary['warping']). See docs/conventions/DEPRECATED_ELEMENT_TYPES.md."
    ),
    "NonlinearWarpingTimoshenkoBeamElement3D": (
        "Element type 'NonlinearWarpingTimoshenkoBeamElement3D' is not registered (legacy alias removed). "
        "Use 'NonlinearTimoshenkoBeamElement3D' with [warping]=1 in element.txt "
        "(or set element_dictionary['warping']) for TL Timoshenko + Vlasov warping (14 local DOFs); "
        "omit warping for 12-DOF TL Timoshenko only. "
        "See docs/conventions/DEPRECATED_ELEMENT_TYPES.md."
    ),
    "NonlinearWarpingEulerBernoulliBeamElement3D": (
        "Element type 'NonlinearWarpingEulerBernoulliBeamElement3D' was removed. "
        "Use 'NonlinearEulerBernoulliBeamElement3D' with [warping]=1 in element.txt "
        "(or set element_dictionary['warping']). "
        "See docs/conventions/DEPRECATED_ELEMENT_TYPES.md."
    ),
}


def ensure_element_type_allowed(element_type: str) -> None:
    """
    Emit ``DeprecationWarning`` and raise ``ValueError`` with a migration message if
    *element_type* was removed or deregistered.
    """
    msg = REMOVED_ELEMENT_TYPES.get(element_type)
    if msg is not None:
        warnings.warn(msg, DeprecationWarning, stacklevel=2)
        raise ValueError(msg)

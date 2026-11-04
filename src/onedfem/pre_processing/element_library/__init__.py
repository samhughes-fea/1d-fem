# pre_processing/element_library/__init__.py
#
# Public API: ElementFactory, Element1DBase, and linear element classes (see __all__).
# Nonlinear element classes (NonlinearEulerBernoulliBeamElement3D, NonlinearTimoshenkoBeamElement3D)
# are not re-exported here; they are intended for use via ElementFactory only (ELEMENT_CLASS_MAP).
# Instantiate elements by type string through ElementFactory.create_elements_batch(...).

from .element_factory import ElementFactory
from .element_1D_base import Element1DBase

# Linear elements (re-export from linear subpackage)
from .linear import (
    LinearEulerBernoulliBeamElement3D,
    LinearTimoshenkoBeamElement3D,
    LinearLevinsonBeamElement3D,
    LinearTrussElement3D,
    LinearBarElement3D,
)

__all__ = [
    "ElementFactory",
    "Element1DBase",
    "LinearEulerBernoulliBeamElement3D",
    "LinearTimoshenkoBeamElement3D",
    "LinearLevinsonBeamElement3D",
    "LinearTrussElement3D",
    "LinearBarElement3D",
]
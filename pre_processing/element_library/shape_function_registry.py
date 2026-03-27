# pre_processing/element_library/shape_function_registry.py
"""
Registry mapping element type → ShapeFunctionOperator factory.

This is the single point of control for which shape-function implementation
each element type uses. Elements obtain their operator via get_shape_function_operator()
instead of hardcoding imports. See docs/element_library/shape_function_conventions.md.
"""

from typing import Any, Callable

# Type: factory that takes element_length and returns a ShapeFunctionOperator-like instance
ShapeFunctionOperatorFactory = Callable[[float], Any]


def _eb_operator(L: float) -> Any:
    from pre_processing.element_library.linear.beam.zero_order_shear_deformation_theory.euler_bernoulli.utilities.shape_functions import (
        ShapeFunctionOperator,
    )
    return ShapeFunctionOperator(element_length=L)


def _timoshenko_operator(L: float) -> Any:
    from pre_processing.element_library.linear.beam.first_order_shear_deformation_theory.timoshenko.utilities.shape_functions import (
        ShapeFunctionOperator,
    )
    return ShapeFunctionOperator(element_length=L)


def _levinson_operator(L: float) -> Any:
    from pre_processing.element_library.linear.beam.third_order_shear_deformation_theory.levinson.utilities.shape_functions import (
        ShapeFunctionOperator,
    )
    return ShapeFunctionOperator(element_length=L)


def _reddy_operator(L: float) -> Any:
    from pre_processing.element_library.linear.beam.third_order_shear_deformation_theory.reddy.utilities.shape_functions import (
        ShapeFunctionOperator,
    )
    return ShapeFunctionOperator(element_length=L)


def _bar_operator(L: float) -> Any:
    from pre_processing.element_library.linear.bar.utilities.shape_functions import (
        ShapeFunctionOperator,
    )
    return ShapeFunctionOperator(element_length=L)


def _truss_operator(L: float) -> Any:
    from pre_processing.element_library.linear.truss.utilities.shape_functions import (
        ShapeFunctionOperator,
    )
    return ShapeFunctionOperator(element_length=L)


# Keys must match class names in element_factory.ELEMENT_CLASS_MAP (plus nonlinear aliases below).
# Warping elements use same shape functions as base (12 DOF); warping DOFs use linear Lagrange in-element.
# GEBT uses Timoshenko 12-DOF shapes at the same quadrature points as linear Timoshenko.
SHAPE_FUNCTION_REGISTRY: dict[str, ShapeFunctionOperatorFactory] = {
    "LinearEulerBernoulliBeamElement3D": _eb_operator,
    "LinearWarpingEulerBernoulliBeamElement3D": _eb_operator,
    "NonlinearEulerBernoulliBeamElement3D": _eb_operator,
    "LinearTimoshenkoBeamElement3D": _timoshenko_operator,
    "LinearWarpingTimoshenkoBeamElement3D": _timoshenko_operator,
    "LinearCurvedTimoshenkoBeamElement3D": _timoshenko_operator,
    "NonlinearTimoshenkoBeamElement3D": _timoshenko_operator,
    "GEBTShearBeamElement3D": _timoshenko_operator,
    "LinearLevinsonBeamElement3D": _levinson_operator,
    "LinearReddyBeamElement3D": _reddy_operator,
    "LinearBarElement3D": _bar_operator,
    "LinearTrussElement3D": _truss_operator,
}


def get_shape_function_operator(element_type: str, element_length: float) -> Any:
    """
    Return a ShapeFunctionOperator for the given element type and length.

    Parameters
    ----------
    element_type : str
        Class name as in element_factory.ELEMENT_CLASS_MAP (e.g.
        "LinearEulerBernoulliBeamElement3D").
    element_length : float
        Element length L (physical length in reference configuration).

    Returns
    -------
    ShapeFunctionOperator (or equivalent)
        Instance for the canonical implementation of that element's convention.

    Raises
    ------
    KeyError
        If element_type is not registered.
    """
    if element_type not in SHAPE_FUNCTION_REGISTRY:
        raise KeyError(
            f"Unknown element type for shape functions: {element_type}. "
            f"Registered: {list(SHAPE_FUNCTION_REGISTRY.keys())}"
        )
    return SHAPE_FUNCTION_REGISTRY[element_type](element_length)

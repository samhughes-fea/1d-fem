"""
FormulationResultSet Container

Caches element formulation data with Gauss-level resolution for efficient
post-processing stress/strain recovery without recomputing shape functions.

This container stores:
- ElementObject: Cached K_e and Gauss point data (B, D, J, N) for each element
- ForceObject: Cached F_e and Gauss point data (N, J, distributed loads) for each element

Benefits:
- Performance: No shape function recomputation during results computation
- Traceability: Element type tracking for mixed-element meshes
- Accuracy: Use exact same B, D matrices as used in stiffness assembly
"""

import logging
import os
from dataclasses import dataclass, field
from typing import List
from pre_processing.element_library.gauss_point_data import ElementObject, ForceObject

logger = logging.getLogger(__name__)

STRICT_SHAPE_FUNCTIONS_ENV_VAR = "FEM_FORMULATION_CACHE_STRICT_SHAPE"


def strict_shape_functions_validation_from_env() -> bool:
    """
    If True, formulation-cache validation raises when any Gauss record lacks
    ``shape_functions`` / ``shape_derivatives`` (stiffness).

    Enable with environment variable ``FEM_FORMULATION_CACHE_STRICT_SHAPE=1``
    (also ``true`` / ``yes`` / ``on``, case-insensitive).
    """
    return os.environ.get(STRICT_SHAPE_FUNCTIONS_ENV_VAR, "").strip().lower() in (
        "1",
        "true",
        "yes",
        "on",
    )


def validate_shape_functions_populated(
    element_objects: List[ElementObject],
    force_objects: List[ForceObject],
    *,
    strict: bool = False,
) -> None:
    """
    Check that all element and force Gauss data have shape_functions (and
    shape_derivatives for stiffness) populated. Log warnings or raise if not.

    Parameters
    ----------
    element_objects : List[ElementObject]
        Element formulation cache.
    force_objects : List[ForceObject]
        Force formulation cache.
    strict : bool, optional
        If True, raise on first missing shape data; otherwise only log a warning.
        Default False. Runners may set this from
        :func:`strict_shape_functions_validation_from_env` (``FEM_FORMULATION_CACHE_STRICT_SHAPE``).
    """
    for obj in element_objects:
        for i, gp in enumerate(obj.gauss_data):
            if gp.shape_functions is None:
                msg = (
                    f"ElementObject element_id={obj.element_id} gauss_data[{i}]: "
                    "shape_functions is None (required for results pipeline)."
                )
                if strict:
                    raise ValueError(msg)
                logger.warning(msg)
            if gp.shape_derivatives is None:
                msg = (
                    f"ElementObject element_id={obj.element_id} gauss_data[{i}]: "
                    "shape_derivatives is None (required for results pipeline)."
                )
                if strict:
                    raise ValueError(msg)
                logger.warning(msg)
    for obj in force_objects:
        for i, gp in enumerate(obj.gauss_data):
            if gp.shape_functions is None:
                msg = (
                    f"ForceObject element_id={obj.element_id} gauss_data[{i}]: "
                    "shape_functions is None (required for results pipeline)."
                )
                if strict:
                    raise ValueError(msg)
                logger.warning(msg)


@dataclass
class FormulationResultSet:
    """
    Caches element formulation data with Gauss-level resolution.
    
    Attributes:
        element_objects: List of ElementObject containing K_e and Gauss data
        force_objects: List of ForceObject containing F_e and force Gauss data
    """
    element_objects: List[ElementObject] = field(default_factory=list)
    force_objects: List[ForceObject] = field(default_factory=list)


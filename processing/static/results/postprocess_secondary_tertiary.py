# processing/static/results/postprocess_secondary_tertiary.py
"""
Secondary and tertiary results pipelines driven by a FormulationResultSet.

Used by static runners and optionally by modal/dynamic runners when
``simulation_settings['post_processing']`` enables strain/stress recovery from a
reference displacement field ``U_global``.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import TYPE_CHECKING, Any, List, Tuple

from processing.static.results.compute_secondary.secondary_results_orchestrator import (
    SecondaryResultsOrchestrator,
)
from processing.static.results.compute_tertiary.tertiary_results_orchestrator import (
    TertiaryResultsOrchestrator,
)
from processing.static.results.save_secondary_container import (
    SaveSecondaryResults,
    SaveSecondaryResultsSummary,
)
from processing.static.results.save_tertiary_container import (
    SaveTertiaryResults,
    SaveTertiaryResultsSummary,
)

if TYPE_CHECKING:
    import numpy as np
    from processing.static.results.containers.formulation_results import (
        FormulationResultSet,
    )
    from processing.static.results.containers.container_hopper import SecondaryResultSet
    from processing.static.results.containers.tertiary_results import TertiaryResults

logger = logging.getLogger(__name__)


def compute_secondary_result_set(
    *,
    elements: List[Any],
    grid_dictionary: dict,
    element_dictionary: dict,
    material_dictionary: dict,
    section_dictionary: dict,
    U_global: "np.ndarray",
    formulation_cache: "FormulationResultSet",
    job_results_dir: Path | str | None,
) -> "SecondaryResultSet":
    orchestrator = SecondaryResultsOrchestrator(
        elements=elements,
        grid_dictionary=grid_dictionary,
        element_dictionary=element_dictionary,
        material_dictionary=material_dictionary,
        section_dictionary=section_dictionary,
        global_displacement=U_global,
        formulation_cache=formulation_cache,
        job_results_dir=Path(job_results_dir) if job_results_dir else None,
    )
    return orchestrator.compute_all()


def save_secondary_outputs(
    secondary_results_set: "SecondaryResultSet",
    secondary_results_dir: Path | str,
) -> None:
    saver = SaveSecondaryResults(
        secondary_results=secondary_results_set,
        save_dir=str(secondary_results_dir),
    )
    saver.save_all()
    SaveSecondaryResultsSummary(
        secondary_results=secondary_results_set,
        save_dir=str(secondary_results_dir),
    ).save()


def compute_tertiary_result_set(
    *,
    secondary_results_set: "SecondaryResultSet",
    formulation_cache: "FormulationResultSet",
    element_dictionary: dict,
    grid_dictionary: dict,
    tertiary_results_dir: Path | str,
) -> "TertiaryResults":
    orchestrator = TertiaryResultsOrchestrator(
        secondary_results=secondary_results_set,
        formulation_cache=formulation_cache,
        element_dictionary=element_dictionary,
        grid_dictionary=grid_dictionary,
        job_results_dir=str(tertiary_results_dir),
    )
    return orchestrator.compute()


def save_tertiary_outputs(
    tertiary_results: "TertiaryResults",
    results_root: Path | str,
) -> None:
    SaveTertiaryResults(
        tertiary_results=tertiary_results,
        save_dir=str(results_root),
    ).save_all()
    SaveTertiaryResultsSummary(
        tertiary_results=tertiary_results,
        save_dir=str(results_root),
    ).save()


def run_secondary_tertiary_from_formulation_cache(
    *,
    elements: List[Any],
    grid_dictionary: dict,
    element_dictionary: dict,
    material_dictionary: dict,
    section_dictionary: dict,
    U_global: "np.ndarray",
    formulation_cache: "FormulationResultSet",
    results_root: Path | str,
    secondary_results_dir: Path | str,
    tertiary_results_dir: Path | str,
) -> Tuple["SecondaryResultSet", "TertiaryResults"]:
    """
    Full secondary then tertiary compute + save (matches static runner layout).
    """
    logger.info("Computing secondary results from formulation cache (modal/dynamic/static-style path)...")
    secondary = compute_secondary_result_set(
        elements=elements,
        grid_dictionary=grid_dictionary,
        element_dictionary=element_dictionary,
        material_dictionary=material_dictionary,
        section_dictionary=section_dictionary,
        U_global=U_global,
        formulation_cache=formulation_cache,
        job_results_dir=secondary_results_dir,
    )
    save_secondary_outputs(secondary, secondary_results_dir)

    logger.info("Computing tertiary results...")
    tertiary = compute_tertiary_result_set(
        secondary_results_set=secondary,
        formulation_cache=formulation_cache,
        element_dictionary=element_dictionary,
        grid_dictionary=grid_dictionary,
        tertiary_results_dir=tertiary_results_dir,
    )
    save_tertiary_outputs(tertiary, results_root)
    logger.info("Secondary and tertiary results saved under %s", results_root)
    return secondary, tertiary

# pre_processing/element_library/element_factory.py

import importlib
import logging
import multiprocessing
import os
from typing import Any, Dict, List, Optional, TYPE_CHECKING

import numpy as np

from pre_processing.element_library.removed_element_types import ensure_element_type_allowed

if TYPE_CHECKING:  # avoid runtime circular import
    from pre_processing.element_library.element_1D_base import Element1DBase


def _corotational_tangent_mode_from_settings(
    simulation_settings: Optional[Dict[str, Any]],
) -> Optional[str]:
    """Resolve ``simulation_settings['nonlinear']['corotational_tangent_mode']`` for ``CorotationalBeamElement3D``."""
    if not simulation_settings:
        return None
    nl = simulation_settings.get("nonlinear")
    if not isinstance(nl, dict):
        return None
    raw = nl.get("corotational_tangent_mode")
    if raw is None:
        return None
    s = str(raw).strip().lower().replace("-", "_")
    if s in ("finite_difference", "fd", "finite_diff"):
        return "finite_difference"
    if s in ("elastic_material", "elastic", "material"):
        return "elastic_material"
    raise ValueError(
        "simulation_settings['nonlinear']['corotational_tangent_mode'] must be "
        "'finite_difference' or 'elastic_material', "
        f"got {raw!r}"
    )


def _inject_corotational_tangent_kwarg(
    etype: str,
    params: Dict[str, Any],
    corotational_tangent_mode: Optional[str],
) -> None:
    if etype == "CorotationalBeamElement3D" and corotational_tangent_mode is not None:
        params["tangent_stiffness_mode"] = corotational_tangent_mode


def _gesdb_tl_fallback_from_settings(
    simulation_settings: Optional[Dict[str, Any]],
) -> bool:
    """Resolve ``simulation_settings['nonlinear']['gesdb_tl_fallback']`` for GESDB elements."""
    if not simulation_settings:
        return False
    nl = simulation_settings.get("nonlinear")
    if not isinstance(nl, dict):
        return False
    raw = nl.get("gesdb_tl_fallback")
    if raw is None:
        return False
    if isinstance(raw, bool):
        return raw
    if isinstance(raw, str):
        return raw.strip().lower() in ("true", "1", "yes", "on")
    return bool(raw)


def _gesdb_kernel_from_settings(simulation_settings: Optional[Dict[str, Any]]) -> str:
    """Resolve ``simulation_settings['nonlinear']['gesdb_kernel']`` (tl_locked | native)."""
    if not simulation_settings:
        return "tl_locked"
    nl = simulation_settings.get("nonlinear")
    if not isinstance(nl, dict):
        return "tl_locked"
    raw = nl.get("gesdb_kernel")
    if raw is None:
        return "tl_locked"
    s = str(raw).strip().lower().replace("-", "_")
    if s in ("tl_locked", "tl", "locked", "chord_tl"):
        return "tl_locked"
    if s in ("native", "engineering", "gesdb_native"):
        return "native"
    raise ValueError(
        "simulation_settings['nonlinear']['gesdb_kernel'] must be 'tl_locked' or 'native', "
        f"got {raw!r}"
    )


def _inject_gesdb_kwargs(
    etype: str,
    params: Dict[str, Any],
    gesdb_tl_fallback: bool,
    gesdb_kernel: str,
) -> None:
    if etype == "GeometricallyExactShearDeformableBeam3D":
        params["gesdb_tl_fallback"] = bool(gesdb_tl_fallback)
        params["gesdb_kernel"] = gesdb_kernel


class ElementFactory:
    """
    Batch-builds 1-D element objects and validates their logging infrastructure.

    Notes
    -----
    * Concrete element classes (e.g. ``LinearEulerBernoulliBeamElement3D``) **must**
      inherit from
      :class:`pre_processing.element_library.element_1D_base.Element1DBase`.
    * Classes are discovered via :pyattr:`ELEMENT_CLASS_MAP` (combined from
      :pyattr:`LINEAR_ELEMENT_CLASS_MAP` and :pyattr:`NONLINEAR_ELEMENT_CLASS_MAP`).
    * Per-element logs are written to ``<job_results_dir>/logs``.
    * Optional ``precurvature_per_element`` (shape ``(N_e, 3)``) is merged into
      a copy of ``element_dictionary`` as key ``"precurvature_per_element"`` (rows
      ``[k_x0, k_y0, k_z0]`` in element-id order). Each :class:`Element1DBase` stores
      ``_E_0_voigt``; beam theories add ``B.T @ D @ E_0`` to ``F_e`` (and use ``E - E_0``
      in nonlinear internal force / ``K_sigma``) where the formulation uses that tensor layout.
    """

    LINEAR_ELEMENT_CLASS_MAP = {
        "LinearEulerBernoulliBeamElement3D":
            "pre_processing.element_library.linear.beam.zero_order_shear_deformation_theory.euler_bernoulli.linear_euler_bernoulli_3D",
        "LinearTimoshenkoBeamElement3D":
            "pre_processing.element_library.linear.beam.first_order_shear_deformation_theory.timoshenko.linear_timoshenko_3D",
        "LinearLevinsonBeamElement3D":
            "pre_processing.element_library.linear.beam.third_order_shear_deformation_theory.levinson.linear_levinson_3D",
        "LinearReddyBeamElement3D":
            "pre_processing.element_library.linear.beam.third_order_shear_deformation_theory.reddy.linear_reddy_3D",
        "LinearTrussElement3D":
            "pre_processing.element_library.linear.truss.linear_truss_3D",
        "LinearBarElement3D":
            "pre_processing.element_library.linear.bar.linear_bar_3D",
    }

    NONLINEAR_ELEMENT_CLASS_MAP = {
        "NonlinearEulerBernoulliBeamElement3D":
            "pre_processing.element_library.nonlinear.euler_bernoulli.nonlinear_euler_bernoulli_3D",
        "NonlinearTimoshenkoBeamElement3D":
            "pre_processing.element_library.nonlinear.timoshenko.nonlinear_timoshenko_3D",
        "NonlinearLevinsonBeamElement3D":
            "pre_processing.element_library.nonlinear.levinson.nonlinear_levinson_3D",
        "NonlinearReddyBeamElement3D":
            "pre_processing.element_library.nonlinear.reddy.nonlinear_reddy_3D",
        "UpdatedLagrangianEulerBernoulliBeamElement3D":
            "pre_processing.element_library.nonlinear.euler_bernoulli.updated_lagrangian_euler_bernoulli_3D",
        "UpdatedLagrangianTimoshenkoBeamElement3D":
            "pre_processing.element_library.nonlinear.timoshenko.updated_lagrangian_timoshenko_3D",
        "CorotationalBeamElement3D":
            "pre_processing.element_library.nonlinear.large_rotations.corotational.corotational_3D",
        "GeometricallyExactShearDeformableBeam3D":
            "pre_processing.element_library.nonlinear.large_rotations.geometrically_exact_shear_deformable_beam.geometrically_exact_shear_deformable_beam_3D",
        "GEBTUnshearableBeamElement3D":
            "pre_processing.element_library.nonlinear.large_rotations.gebt_unshearable.gebt_unshearable_3D",
    }

    ELEMENT_CLASS_MAP = {**LINEAR_ELEMENT_CLASS_MAP, **NONLINEAR_ELEMENT_CLASS_MAP}

    # ------------------------------------------------------------------ #
    def __init__(self, job_results_dir: str) -> None:
        """
        Parameters
        ----------
        job_results_dir
            Root directory where each element’s log files will be created.
        """
        self.job_results_dir = job_results_dir
        self.logger = self._setup_logger()

    # ------------------------------------------------------------------ #
    def create_elements_batch(
        self,
        *,
        element_ids: np.ndarray,
        element_dictionary: Dict[str, Any],
        grid_dictionary: Dict[str, Any],
        material_dictionary: Dict[str, Any],
        section_dictionary: Dict[str, Any],
        point_load_array: np.ndarray,
        distributed_load_array: np.ndarray,
        precurvature_per_element: Optional[np.ndarray] = None,
        simulation_settings: Optional[Dict[str, Any]] = None,
        enable_parallel: bool = False,
        num_processes: Optional[int] = None,
    ) -> List["Element1DBase"]:
        """
        Instantiate one element object for every ID in *element_ids*.

        Parameters
        ----------
        element_ids
            NumPy ``int64`` array of element IDs. Order must match
            ``element_dictionary["types"]`` and ``["connectivity"]``.
        element_dictionary
            Parsed output from :class:`ElementParser`.
            Required keys: ``"ids"``, ``"types"``, ``"connectivity"``.
        grid_dictionary
            Parsed output from :class:`GridParser`.
        material_dictionary
            Parsed output from :class:`MaterialParser`.
        section_dictionary
            Parsed output from :class:`SectionParser`.
        point_load_array
            ``(N, 4)`` array or empty array of point-load data.
        distributed_load_array
            ``(N, 5)`` array or empty array of distributed-load data.
        precurvature_per_element
            Optional ``(N_e, 3)`` array of ``[k_x0, k_y0, k_z0]`` (1/m) per row,
            aligned with *element_ids*. Omitted or ``None``: straight beams behave
            as zero reference curvature.
        simulation_settings
            Optional resolved settings dict. When ``simulation_settings['nonlinear']['corotational_tangent_mode']``
            is set, :class:`CorotationalBeamElement3D` receives ``tangent_stiffness_mode`` (``finite_difference`` or
            ``elastic_material``).             When ``simulation_settings['nonlinear']['gesdb_tl_fallback']`` is true,
            :class:`GeometricallyExactShearDeformableBeam3D` receives ``gesdb_tl_fallback=True``.
            ``gesdb_kernel`` (``tl_locked`` or ``native``) is passed when present.
        enable_parallel : bool, optional
            Enable parallel instantiation (default False).
        num_processes : int, optional
            Number of processes for parallel instantiation. If None or "auto",
            uses os.cpu_count(). Ignored if enable_parallel=False.

        Returns
        -------
        list[Element1DBase]
            One instance per ID in *element_ids*.
        """
        self.logger.info("🚩 Starting batch element creation.")

        element_ids_array = self._sanitize_element_ids(element_ids)
        self._validate_element_dictionary(element_dictionary)

        element_types_array = np.asarray(element_dictionary["types"])
        if len(element_types_array) != len(element_ids_array):
            self.logger.error(
                "Element ID / type count mismatch: IDs=%d, types=%d",
                len(element_ids_array), len(element_types_array),
            )
            raise ValueError("Mismatch between element IDs and types.")

        element_dictionary_use = self._merge_precurvature_dictionary(
            element_dictionary, precurvature_per_element
        )
        ctm = _corotational_tangent_mode_from_settings(simulation_settings)
        gtfb = _gesdb_tl_fallback_from_settings(simulation_settings)
        gkernel = _gesdb_kernel_from_settings(simulation_settings)

        # Determine if we should use parallel processing
        num_elements = len(element_ids_array)
        threshold = 50  # Minimum elements for parallel
        
        use_parallel = (
            enable_parallel and 
            num_elements >= threshold and
            num_processes != 1
        )
        
        if use_parallel:
            # Handle "auto" for num_processes
            if num_processes is None or num_processes == "auto":
                num_processes = os.cpu_count() or 1
            
            try:
                elements = self._instantiate_elements_parallel(
                    element_ids_array,
                    element_types_array,
                    element_dictionary_use,
                    grid_dictionary,
                    material_dictionary,
                    section_dictionary,
                    point_load_array,
                    distributed_load_array,
                    num_processes,
                    ctm,
                    gtfb,
                    gkernel,
                )
            except Exception as e:
                self.logger.warning(
                    f"Parallel instantiation failed ({type(e).__name__}: {e}), "
                    "falling back to sequential"
                )
                elements = self._instantiate_elements_sequential(
                    element_ids_array,
                    element_types_array,
                    element_dictionary_use,
                    grid_dictionary,
                    material_dictionary,
                    section_dictionary,
                    point_load_array,
                    distributed_load_array,
                    ctm,
                    gtfb,
                    gkernel,
                )
        else:
            elements = self._instantiate_elements_sequential(
                element_ids_array,
                element_types_array,
                element_dictionary_use,
                grid_dictionary,
                material_dictionary,
                section_dictionary,
                point_load_array,
                distributed_load_array,
                ctm,
                gtfb,
                gkernel,
            )

        self._validate_logging_infrastructure(elements)
        self.logger.info("✅ Batch element creation complete.")
        return elements

    # ------------------------------------------------------------------ #
    # PRIVATE HELPERS
    # ------------------------------------------------------------------ #
    def _setup_logger(self) -> logging.Logger:
        """
        Create a dedicated file logger ``ElementFactory.log`` under
        ``<job_results_dir>/logs``. Duplicate handlers are purged to avoid
        repeated messages in multi-process runs.
        """
        log_dir = os.path.join(self.job_results_dir, "logs")
        os.makedirs(log_dir, exist_ok=True)

        logger = logging.getLogger("ElementFactory")
        logger.setLevel(logging.DEBUG)
        if logger.handlers:
            logger.handlers.clear()

        fh = logging.FileHandler(
            os.path.join(log_dir, "ElementFactory.log"), mode="w", encoding="utf-8"
        )
        fh.setFormatter(logging.Formatter("%(asctime)s [%(levelname)s]: %(message)s"))
        logger.addHandler(fh)
        logger.propagate = False

        logger.debug("🟢 ElementFactory logger initialised.")
        return logger

    # .................................................................. #
    @staticmethod
    def _merge_precurvature_dictionary(
        element_dictionary: Dict[str, Any],
        precurvature_per_element: Optional[np.ndarray],
    ) -> Dict[str, Any]:
        """Return *element_dictionary* or a shallow copy with ``precurvature_per_element`` set."""
        if precurvature_per_element is None:
            return element_dictionary
        n_e = len(element_dictionary["ids"])
        arr = np.asarray(precurvature_per_element, dtype=np.float64).reshape(n_e, 3)
        if arr.shape != (n_e, 3):
            raise ValueError(
                f"precurvature_per_element must have shape ({n_e}, 3), got {arr.shape}"
            )
        return {**element_dictionary, "precurvature_per_element": arr}

    # .................................................................. #
    def _sanitize_element_ids(self, raw_ids) -> np.ndarray:
        """
        Convert *raw_ids* to ``int64`` NumPy array and verify positivity.

        Raises
        ------
        ValueError
            If any ID is negative or non-numeric.
        """
        try:
            ids = np.asarray(raw_ids, dtype=np.int64)
            if np.any(ids < 0):
                raise ValueError("Negative element IDs detected.")
            self.logger.debug("Element IDs sanitised: %s", ids.tolist())
            return ids
        except Exception as exc:
            self.logger.error("Invalid element IDs.", exc_info=True)
            raise ValueError("Invalid element IDs") from exc

    # .................................................................. #
    def _validate_element_dictionary(self, element_dict: Dict[str, Any]) -> None:
        """
        Ensure *element_dict* has the mandatory keys ``ids``, ``types``,
        ``connectivity``.
        """
        required = {"types", "ids", "connectivity"}
        missing = required - element_dict.keys()
        if missing:
            self.logger.error("Missing keys in element dictionary: %s", missing)
            raise KeyError(f"Missing keys: {missing}")

    # .................................................................. #
    def _load_element_modules(self, element_types: np.ndarray) -> Dict[str, Any]:
        """
        Dynamically import every unique element module referenced in
        *element_types*.

        Returns
        -------
        dict
            Mapping ``element_type → imported module``.
        """
        modules: Dict[str, Any] = {}
        for etype in np.unique(element_types):
            ensure_element_type_allowed(str(etype))
            if etype not in self.ELEMENT_CLASS_MAP:
                self.logger.error("Unregistered element type: %s", etype)
                raise ValueError(f"Unregistered element type: {etype}")

            module_path = self.ELEMENT_CLASS_MAP[etype]
            try:
                module = importlib.import_module(module_path)
                if not hasattr(module, etype):
                    raise AttributeError(
                        f"Module {module_path} missing class {etype}"
                    )
                modules[etype] = module
                self.logger.debug("Module for %s loaded.", etype)
            except Exception as exc:
                self.logger.error("Cannot load module %s", module_path, exc_info=True)
                raise ImportError(f"Cannot load module {module_path}") from exc
        return modules

    # .................................................................. #
    def _instantiate_elements_sequential(
        self,
        element_ids_array: np.ndarray,
        element_types_array: np.ndarray,
        element_dictionary: Dict[str, Any],
        grid_dictionary: Dict[str, Any],
        material_dictionary: Dict[str, Any],
        section_dictionary: Dict[str, Any],
        point_load_array: np.ndarray,
        distributed_load_array: np.ndarray,
        corotational_tangent_mode: Optional[str] = None,
        gesdb_tl_fallback: bool = False,
        gesdb_kernel: str = "tl_locked",
    ) -> List["Element1DBase"]:
        """Sequential element instantiation (original implementation)."""
        modules = self._load_element_modules(element_types_array)
        
        elements: List["Element1DBase"] = []
        for etype, eid in zip(element_types_array, element_ids_array):
            params = {
                "element_id":             int(eid),
                "element_dictionary":     element_dictionary,
                "grid_dictionary":        grid_dictionary,
                "material_dictionary":    material_dictionary,
                "section_dictionary":     section_dictionary,
                "point_load_array":       point_load_array,
                "distributed_load_array": distributed_load_array,
                "job_results_dir":        self.job_results_dir,
            }
            _inject_corotational_tangent_kwarg(etype, params, corotational_tangent_mode)
            _inject_gesdb_kwargs(etype, params, gesdb_tl_fallback, gesdb_kernel)
            elem = self._instantiate_element(etype, eid, params, modules)
            elements.append(elem)
            self.logger.debug("✅ Element %s (%s) instantiated.", eid, etype)
        
        return elements
    
    # .................................................................. #
    def _instantiate_elements_parallel(
        self,
        element_ids_array: np.ndarray,
        element_types_array: np.ndarray,
        element_dictionary: Dict[str, Any],
        grid_dictionary: Dict[str, Any],
        material_dictionary: Dict[str, Any],
        section_dictionary: Dict[str, Any],
        point_load_array: np.ndarray,
        distributed_load_array: np.ndarray,
        num_processes: int,
        corotational_tangent_mode: Optional[str] = None,
        gesdb_tl_fallback: bool = False,
        gesdb_kernel: str = "tl_locked",
    ) -> List["Element1DBase"]:
        """Parallel element instantiation using multiprocessing."""
        self.logger.info(f"Instantiating {len(element_ids_array)} elements in parallel ({num_processes} processes)")
        
        # Prepare arguments for worker function
        args = [
            (
                int(eid),
                str(etype),
                element_dictionary,
                grid_dictionary,
                material_dictionary,
                section_dictionary,
                point_load_array,
                distributed_load_array,
                self.job_results_dir,
                corotational_tangent_mode,
                gesdb_tl_fallback,
                gesdb_kernel,
                i  # index for ordering
            )
            for i, (eid, etype) in enumerate(zip(element_ids_array, element_types_array))
        ]
        
        try:
            with multiprocessing.Pool(processes=num_processes) as pool:
                results = pool.map(_worker_instantiate_element, args)
            
            # Sort by index and extract elements
            results.sort(key=lambda x: x[0])
            elements = [r[1] for r in results]
            
            self.logger.info("Parallel element instantiation completed successfully")
            return elements
            
        except (multiprocessing.PicklingError, AttributeError, TypeError) as e:
            raise RuntimeError(f"Parallel instantiation failed: {e}") from e
    
    # .................................................................. #
    def _instantiate_element(
        self,
        etype: str,
        eid: int,
        params: Dict[str, Any],
        modules: Dict[str, Any],
    ) -> "Element1DBase":
        """
        Instantiate a single element given its type *etype* and constructor
        *params* pulled from the various dictionaries.
        """
        from pre_processing.element_library.element_1D_base import Element1DBase

        cls = getattr(modules[etype], etype)
        if not issubclass(cls, Element1DBase):
            raise TypeError(f"Class {etype} must inherit from Element1DBase.")

        return cls(**params)

    # .................................................................. #
    def _validate_logging_infrastructure(self, elements: List["Element1DBase"]) -> None:
        """
        Smoke-test each element’s `logger_operator` by writing a dummy matrix
        and confirming a non-empty file on disk.

        Raises
        ------
        RuntimeError
            If any element fails the logging check.
        """
        import numpy as np  # local to keep top-level deps minimal

        for elem in elements:
            try:
                log_op = elem.logger_operator
                log_op.log_matrix("stiffness", np.zeros((2, 2)), {"test": True})
                log_op.flush_all()
                path = log_op._get_log_path("stiffness")
                if not os.path.isfile(path) or os.path.getsize(path) == 0:
                    raise IOError(f"Logging failed for element {elem.element_id}")
                self.logger.debug("Logging verified for element %s.", elem.element_id)
            except Exception as exc:
                self.logger.error(
                    "Logging infrastructure invalid for element %s",
                    elem.element_id, exc_info=True,
                )
                raise RuntimeError(
                    f"Logging validation failed for element {elem.element_id}"
                ) from exc


# Worker function for parallel element instantiation (must be at module level for pickling)
def _worker_instantiate_element(args):
    """Worker function for parallel element instantiation."""
    (
        eid,
        etype,
        element_dictionary,
        grid_dictionary,
        material_dictionary,
        section_dictionary,
        point_load_array,
        distributed_load_array,
        job_results_dir,
        corotational_tangent_mode,
        gesdb_tl_fallback,
        gesdb_kernel,
        index,
    ) = args
    
    try:
        # Re-import modules in worker process (required for multiprocessing)
        import importlib
        from pre_processing.element_library.element_1D_base import Element1DBase
        from pre_processing.element_library.removed_element_types import ensure_element_type_allowed

        ensure_element_type_allowed(etype)

        # Use same map as main class so Bar/Truss and future elements stay in sync
        module_path = ElementFactory.ELEMENT_CLASS_MAP[etype]
        module = importlib.import_module(module_path)
        cls = getattr(module, etype)
        
        if not issubclass(cls, Element1DBase):
            raise TypeError(f"Class {etype} must inherit from Element1DBase.")
        
        # Instantiate element
        params = {
            "element_id":             eid,
            "element_dictionary":     element_dictionary,
            "grid_dictionary":        grid_dictionary,
            "material_dictionary":    material_dictionary,
            "section_dictionary":     section_dictionary,
            "point_load_array":       point_load_array,
            "distributed_load_array": distributed_load_array,
            "job_results_dir":       job_results_dir,
        }
        _inject_corotational_tangent_kwarg(etype, params, corotational_tangent_mode)
        _inject_gesdb_kwargs(etype, params, gesdb_tl_fallback, gesdb_kernel)
        element = cls(**params)
        
        return (index, element)
    except Exception as e:
        import logging
        logger = logging.getLogger(__name__)
        logger.error(f"Error instantiating element {eid} ({etype}) at index {index}: {e}")
        return (index, None)
# Post-processing

Result files written after a solve: stresses, failure margins, sensitivities and field sampling for plots.

### `post_processing.failure.ply_failure`

[post_processing/failure/ply_failure.py](../../post_processing/failure/ply_failure.py)

Ply failure exposures and margins at the beam elements' Gauss points.

- **`write_ply_failure(job_dir: str, job_results_dir: str, elements: Sequence, element_objects: Sequence, U_global, thin_walled: Optional[Sequence] = None, solid: Optional[Sequence] = None, strengths: Optional[Mapping[int, PlyStrength]] = None, points_per_wall: Optional[int] = None) -> Optional[str]`**
  Write ``ply_failure.csv`` and ``ply_failure_summary.csv``; return the first path, or ``None`` when the job has no ``ply_strength.txt`` (and no ``strengths`` are given).

### `post_processing.sensitivity.layup_sensitivity`

[post_processing/sensitivity/layup_sensitivity.py](../../post_processing/sensitivity/layup_sensitivity.py)

Sensitivities of beam responses to ply angles and thicknesses (``layup_sensitivity.txt``).

- **`LAYUP_SENSITIVITY_COLUMNS`** (constant)
- **`element_derivatives(el, d_section, mass: bool = False) -> Tuple[np.ndarray, Optional[np.ndarray]]`**
  ``dK_e`` (and ``dM_e``) of element ``el`` for the derivative section ``d_section``.
- **`layup_sensitivities(elements, thin_walled, params: Sequence[LayupParameter], *, K, U=None, F=None, fixed_dofs=(), M=None, frequencies_hz=None, mode_shapes=None, responses: Sequence[Response] = (), dof_index=None, load_factors=None, prestress=None) -> List[dict]`**
  Rows of ``LAYUP_SENSITIVITY_COLUMNS`` for a solved static (``U``, ``F``), eigen (``M``, ``frequencies_hz``, ``mode_shapes``) or linear buckling (``load_factors``, ``mode_shapes``, ``prestress`` ``U_0``, ``fixed_dofs``) job.
- **`write_layup_sensitivity(job_dir: str, job_results_dir: str, element_dictionary, elements, runner, solver_type: str) -> Optional[str]`**
  Write ``layup_sensitivity.csv`` for a job with ``layup_sensitivity.txt``; the path, or ``None``.

### `post_processing.solid.gauss_projection`

[post_processing/solid/gauss_projection.py](../../post_processing/solid/gauss_projection.py)

Fields over 2D section elements from their Gauss-point values.

- **`GAUSS_1D`** (constant)
- **`gauss_lagrange_q8(xi: float, eta: float) -> npt.NDArray[np.float64]`**
  Weights (9,) of the nine Gauss values of a ``Q8`` at ``(xi, eta)``, in quadrature order.
- **`project_q8(gauss_values: npt.ArrayLike, xi: float, eta: float) -> npt.NDArray[np.float64]`**
  The projected field at ``(xi, eta)`` from the element's nine Gauss values (9,) or (9, k).
- **`element_patches(nodes: npt.ArrayLike, connectivity: Sequence[int], gauss_values: npt.ArrayLike, n_sub: int = 4) -> Tuple[List[npt.NDArray[np.float64]], npt.NDArray[np.float64]]`**
  Sub-cell polygons of one ``Q8`` and the projected value at each sub-cell centre.
- **`element_line(gauss_values: npt.ArrayLike, xi: float = 0.0, n: int = 21) -> Tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]`**
  ``(eta, value)`` of the projected field along ``xi`` = const, ``eta`` from -1 to 1.

### `post_processing.solid.section_stresses`

[post_processing/solid/section_stresses.py](../../post_processing/solid/section_stresses.py)

Write 3D section stresses of a solved job with sections meshed in 2D.

- **`element_section_rows(element, U_e, analysis, xi: float = 0.0, mode='') -> List[list]`**
  CSV rows (``SECTION_STRESS_COLUMNS``) for one element; ``x`` is the global ``X`` of the axis at ``xi``.
- **`element_interface_rows(element, U_e, analysis, xi: float = 0.0, mode='') -> List[list]`**
  CSV rows (``INTERFACE_STRESS_COLUMNS``): the interlaminar stresses on every ply interface of the section (``SolidSectionAnalysis.interface_stresses``), at the middle of each interface edge.
- **`write_section_stresses(job_results_dir: str, elements: Sequence, analyses: Sequence, U_global, stations: Optional[Sequence] = None, modes: Optional[np.ndarray] = None) -> str`**
  Write ``section_stresses.csv`` and ``interface_stresses.csv`` for the elements that have a 2D section; return the path of the first.

### `post_processing.graphical_visualisers.beam_field_sampling`

[post_processing/graphical_visualisers/beam_field_sampling.py](../../post_processing/graphical_visualisers/beam_field_sampling.py)

Displacement fields along a straight beam from its nodal solution, through the element shape functions.

- **`interpolation_kind(element_type: str) -> str`**
  ``"hermite"`` for Euler-Bernoulli element types, else ``"linear"`` (as ``_bending_kind``).
- **`beam_curve(x_nodes: Sequence[float], U: npt.ArrayLike, component: int, element_type: str, n_per_element: int = 11) -> Tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]`**
  ``(x, value)`` of one field component along a straight beam on the ``x`` axis.

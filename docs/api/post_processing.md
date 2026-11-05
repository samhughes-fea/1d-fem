# Post-processing

Result files written after a solve: stresses, failure margins, sensitivities and field sampling for plots.

### `onedfem.post_processing.failure.ply_failure`

[src/onedfem/post_processing/failure/ply_failure.py](../../src/onedfem/post_processing/failure/ply_failure.py)

Ply failure exposures and margins at the beam elements' Gauss points.

- **`write_ply_failure(job_dir: str, job_results_dir: str, elements: Sequence, element_objects: Sequence, U_global, thin_walled: Optional[Sequence] = None, solid: Optional[Sequence] = None, strengths: Optional[Mapping[int, PlyStrength]] = None, points_per_wall: Optional[int] = None) -> Optional[str]`**
  Write ``ply_failure.csv`` and ``ply_failure_summary.csv``; return the first path, or ``None`` when the job has no ``ply_strength.txt`` (and no ``strengths`` are given).

### `onedfem.post_processing.sensitivity.layup_sensitivity`

[src/onedfem/post_processing/sensitivity/layup_sensitivity.py](../../src/onedfem/post_processing/sensitivity/layup_sensitivity.py)

Sensitivities of beam responses to ply angles and thicknesses (``layup_sensitivity.txt``).

- **`LAYUP_SENSITIVITY_COLUMNS`** (constant)
- **`element_derivatives(el, d_section, mass: bool = False) -> Tuple[np.ndarray, Optional[np.ndarray]]`**
  ``dK_e`` (and ``dM_e``) of element ``el`` for the derivative section ``d_section``.
- **`layup_sensitivities(elements, thin_walled, params: Sequence[LayupParameter], *, K, U=None, F=None, fixed_dofs=(), M=None, frequencies_hz=None, mode_shapes=None, responses: Sequence[Response] = (), dof_index=None, load_factors=None, prestress=None) -> List[dict]`**
  Rows of ``LAYUP_SENSITIVITY_COLUMNS`` for a solved static (``U``, ``F``), eigen (``M``, ``frequencies_hz``, ``mode_shapes``) or linear buckling (``load_factors``, ``mode_shapes``, ``prestress`` ``U_0``, ``fixed_dofs``) job.
- **`write_layup_sensitivity(job_dir: str, job_results_dir: str, element_dictionary, elements, runner, solver_type: str) -> Optional[str]`**
  Write ``layup_sensitivity.csv`` for a job with ``layup_sensitivity.txt``; the path, or ``None``.

### `onedfem.post_processing.solid.element_recovery`

[src/onedfem/post_processing/solid/element_recovery.py](../../src/onedfem/post_processing/solid/element_recovery.py)

3D stresses over a 2D-meshed section at a point of a beam element.

- **`element_section_state(element, U_e: npt.ArrayLike, xi: float, analysis) -> Tuple[Optional[np.ndarray], Optional[np.ndarray]]`**
  ``(forces, benscoter)`` at ``xi``: right-handed section forces (6,) for Euler-Bernoulli and Timoshenko elements, or right-handed Benscoter strains (8,) for Benscoter elements (the other entry ``None``).
- **`element_section_stresses(element, U_e: npt.ArrayLike, xi: float, analysis) -> List`**
  ``PointStress`` at every quadrature point of the section at ``xi`` of an element.

### `onedfem.post_processing.solid.gauss_projection`

[src/onedfem/post_processing/solid/gauss_projection.py](../../src/onedfem/post_processing/solid/gauss_projection.py)

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

### `onedfem.post_processing.solid.section_stresses`

[src/onedfem/post_processing/solid/section_stresses.py](../../src/onedfem/post_processing/solid/section_stresses.py)

Write 3D section stresses of a solved job with sections meshed in 2D.

- **`element_section_rows(element, U_e, analysis, xi: float = 0.0, mode='') -> List[list]`**
  CSV rows (``SECTION_STRESS_COLUMNS``) for one element; ``x`` is the global ``X`` of the axis at ``xi``.
- **`element_interface_rows(element, U_e, analysis, xi: float = 0.0, mode='') -> List[list]`**
  CSV rows (``INTERFACE_STRESS_COLUMNS``): the interlaminar stresses on every ply interface of the section (``SolidSectionAnalysis.interface_stresses``), at the middle of each interface edge.
- **`write_section_stresses(job_results_dir: str, elements: Sequence, analyses: Sequence, U_global, stations: Optional[Sequence] = None, modes: Optional[np.ndarray] = None) -> str`**
  Write ``section_stresses.csv`` and ``interface_stresses.csv`` for the elements that have a 2D section; return the path of the first.

### `onedfem.post_processing.graphical_visualisers.beam_field_sampling`

[src/onedfem/post_processing/graphical_visualisers/beam_field_sampling.py](../../src/onedfem/post_processing/graphical_visualisers/beam_field_sampling.py)

Displacement fields along a straight beam from its nodal solution, through the element shape functions.

- **`interpolation_kind(element_type: str) -> str`**
  ``"hermite"`` for Euler-Bernoulli element types, else ``"linear"`` (as ``_bending_kind``).
- **`beam_curve(x_nodes: Sequence[float], U: npt.ArrayLike, component: int, element_type: str, n_per_element: int = 11) -> Tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]`**
  ``(x, value)`` of one field component along a straight beam on the ``x`` axis.

### `onedfem.post_processing.graphical_visualisers.figure_style`

[src/onedfem/post_processing/graphical_visualisers/figure_style.py](../../src/onedfem/post_processing/graphical_visualisers/figure_style.py)

One figure style for the examples and validation plots, and the marks that say what each value is.

- **`use() -> None`**
  Set the rcParams of the style (Agg backend, 200 dpi, quiet axes).
- **`save(fig, path: str) -> None`**
  Write ``fig`` to ``path`` at the style's resolution and close it.
- **`MODEL_COLOURS`** (constant)
- **`MATERIAL_COLOURS`** (constant)
- **`colour(name: str, i: int = 0) -> str`**
  The model's colour, or the ``i``-th colour of the cycle for any other series.
- **`field(ax, x: npt.ArrayLike, y: npt.ArrayLike, color: str = '#1b1f24', label: Optional[str] = None, **kw)`**
  An element field through its shape functions (``NaN`` breaks between elements).
- **`nodal(ax, x, y, color: str = '#1b1f24', label: Optional[str] = None, **kw)`**
  Explicit nodal values: solid triangles.
- **`gauss(ax, x, y, color: str = '#1b1f24', label: Optional[str] = None, **kw)`**
  Gauss point values: small solid circles.
- **`projected(ax, x, y, color: str = '#1b1f24', label: Optional[str] = None, **kw)`**
  Nodal values projected from Gauss points: hollow circles.
- **`solves(ax, x, y, color: str = '#1b1f24', label: Optional[str] = None, marker: str = 'o', **kw)`**
  One value per solve (load step, speed, length): filled markers, not joined.
- **`analytic(ax, x, y, label: Optional[str] = 'analytic', color: str = '#1b1f24', **kw)`**
  An analytic reference: black dashed line.
- **`reference(ax, x, y, label: Optional[str] = 'reference', color: str = '#1b1f24', **kw)`**
  A computed reference (shell, solid, finite strip): hollow black diamonds.
- **`marks_legend(ax, kinds: Iterable[str], loc: str = 'lower right', labels: Optional[dict] = None, **kw)`**
  A grey legend of what the marks mean, added beside the series legend.
- **`panel_title(ax, title: str, note: Optional[str] = None) -> None`**
  A bold panel title with an optional grey note under it.

### `onedfem.post_processing.graphical_visualisers.thin_walled_drawing`

[src/onedfem/post_processing/graphical_visualisers/thin_walled_drawing.py](../../src/onedfem/post_processing/graphical_visualisers/thin_walled_drawing.py)

Drawings of a thin-walled section and of wall fields around its contour.

- **`draw_section(ax, analysis, materials: Sequence[Tuple[str, object]] = (), nodes: bool = True, angle_shading: bool = True, edge: bool = True, walls: Optional[Iterable[int]] = None, legend: bool = True, zorder: float = 2.0, ply_lines: bool = False)`**
  Draw the walls at true thickness, ply by ply; returns the legend handles by material.
- **`contour_field(ax, analysis, value_at: ValueAt, scale: float, per_element: int = 9, pos_colour: str = '#d4570b', neg_colour: str = '#0b6fb0', alpha: float = 0.35, mark_gauss: bool = True, walls: Optional[Iterable[int]] = None, zorder: float = 3.0)`**
  Diagram of ``value_at(wall_id, s)`` normal to each wall, ``scale`` length per unit value.
- **`unrolled_field(ax, analysis, value_at: ValueAt, wall_order: Sequence[int], color: str = '#1b1f24', label: Optional[str] = None, per_element: int = 9, mark_gauss: bool = True, wall_labels: Optional[Dict[int, str]] = None, shade: bool = True, groups: Optional[Dict[str, Sequence[int]]] = None) -> Array`**
  ``value_at`` against the distance around the contour, walls in ``wall_order``.
- **`wall_samples(analysis, wall_id: int, per_element: int = 9) -> List[Array]`**
  Distances along a wall, ``per_element`` points inside each wall element (ends included).
- **`gauss_s(analysis, wall_id: int) -> Array`**
  Distances along a wall of the section analysis's Gauss points.

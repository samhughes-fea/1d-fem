# Elements

Beam, bar and truss elements and their shared utilities.

## Element types

The `element_type` names a job's `element.txt` accepts.

| Element type | Module | Summary |
|---|---|---|
| `CorotationalBeamElement3D` | `onedfem.pre_processing.element_library.nonlinear.large_rotations.corotational.corotational_3D` | Co-rotational 3D beam (``kernel`` = ``timoshenko`` or ``euler_bernoulli``). |
| `CorotationalBenscoterBeamElement3D` | `onedfem.pre_processing.element_library.nonlinear.large_rotations.corotational.corotational_section` | Corotational Benscoter beam: restrained warping and large twist on thin-walled sections. |
| `CorotationalEulerBernoulliBeamElement3D` | `onedfem.pre_processing.element_library.nonlinear.large_rotations.corotational.corotational_section` | Corotational Euler-Bernoulli beam (shear-rigid core, classical section stiffness). |
| `CorotationalGBTBeamElement3D` | `onedfem.pre_processing.element_library.nonlinear.large_rotations.corotational.corotational_section` | Corotational GBT beam: large rotations of the section, with its warping and distortional modes linear in the element's own frame. |
| `CorotationalLevinsonBeamElement3D` | `onedfem.pre_processing.element_library.nonlinear.large_rotations.corotational.corotational_section` | Corotational Levinson beam (non-variational third-order shear; unsymmetric tangent). |
| `CorotationalReddyBeamElement3D` | `onedfem.pre_processing.element_library.nonlinear.large_rotations.corotational.corotational_section` | Corotational Reddy beam: third-order shear (shear lag of the slaved warping) at large rotation. |
| `CorotationalTimoshenkoBeamElement3D` | `onedfem.pre_processing.element_library.nonlinear.large_rotations.corotational.corotational_section` | Corotational Timoshenko beam on the coupled Saint-Venant section stiffness. |
| `GEBTUnshearableBeamElement3D` | `onedfem.pre_processing.element_library.nonlinear.large_rotations.gebt_unshearable.gebt_unshearable_3D` | Unshearable GEBT: Kirchhoff stretch and curvature, linear EB section law. |
| `GeometricallyExactShearDeformableBeam3D` | `onedfem.pre_processing.element_library.nonlinear.large_rotations.geometrically_exact_shear_deformable_beam.geometrically_exact_shear_deformable_beam_3D` | Factory-facing **GeometricallyExactShearDeformableBeam3D** element type. |
| `LinearBarElement3D` | `onedfem.pre_processing.element_library.linear.bar.linear_bar_3D` | 2 nodes, 6 DOF/node; only axial and torsion strains are non-zero in ``B``. |
| `LinearBenscoterBeamElement3D` | `onedfem.pre_processing.element_library.linear.beam.fosdt.benscoter.linear_benscoter_3D` | Timoshenko beam with Benscoter restrained warping for thin-walled anisotropic sections. |
| `LinearEulerBernoulliBeamElement3D` | `onedfem.pre_processing.element_library.linear.beam.zosdt.euler_bernoulli.linear_euler_bernoulli_3D` | 2-node straight beam, 6 DOF per node, 12 total (see module docstring for ``U_e`` order). |
| `LinearGBTBeamElement3D` | `onedfem.pre_processing.element_library.linear.beam.fosdt.gbt.linear_gbt_3D` | Generalised Beam Theory element for thin-walled anisotropic sections. |
| `LinearLevinsonBeamElement3D` | `onedfem.pre_processing.element_library.linear.beam.tosdt.levinson.linear_levinson_3D` | Levinson third-order shear beam (isotropic or thin-walled sections). |
| `LinearReddyBeamElement3D` | `onedfem.pre_processing.element_library.linear.beam.tosdt.reddy.linear_reddy_3D` | Reddy third-order shear beam (isotropic or thin-walled sections). |
| `LinearTimoshenkoBeamElement3D` | `onedfem.pre_processing.element_library.linear.beam.fosdt.timoshenko.linear_timoshenko_3D` | 2-node straight beam, 6 DOF per node, 12 total; local ``x`` along the chord. |
| `LinearTrussElement3D` | `onedfem.pre_processing.element_library.linear.truss.linear_truss_3D` | 2 nodes, 6 DOF/node; transverse displacement follows the element transverse direction (utilities). |
| `NonlinearEulerBernoulliBeamElement3D` | `onedfem.pre_processing.element_library.nonlinear.euler_bernoulli.total_lagrangian.nonlinear_euler_bernoulli_3D` | TL EB: Green–Lagrange strain, ``S = D @ E``, material and geometric tangents from Gauss sums. |
| `NonlinearLevinsonBeamElement3D` | `onedfem.pre_processing.element_library.nonlinear.tosdt.levinson.nonlinear_levinson_3D` | TL Levinson: Petrov–Galerkin test on six resultants plus GL increments. |
| `NonlinearReddyBeamElement3D` | `onedfem.pre_processing.element_library.nonlinear.tosdt.reddy.nonlinear_reddy_3D` | TL Reddy: variational ``K = int B8.T D B8`` plus GL ``K_delta`` and Hessian ``K_sigma``. |
| `NonlinearTimoshenkoBeamElement3D` | `onedfem.pre_processing.element_library.nonlinear.timoshenko.total_lagrangian.nonlinear_timoshenko_3D` | TL Timoshenko: selective ``K_0`` plus ``K_delta`` and Hessian ``K_sigma``. |
| `UpdatedLagrangianEulerBernoulliBeamElement3D` | `onedfem.pre_processing.element_library.nonlinear.euler_bernoulli.updated_lagrangian.updated_lagrangian_euler_bernoulli_3D` | UL Euler–Bernoulli: linear ``B``, ``D`` on the current chord plus beam-column ``K_sigma``. |
| `UpdatedLagrangianLevinsonBeamElement3D` | `onedfem.pre_processing.element_library.nonlinear.tosdt.levinson.updated_lagrangian.updated_lagrangian_levinson_3D` | UL Levinson: Petrov–Galerkin ``K_0`` on the current chord plus beam-column ``K_sigma``. |
| `UpdatedLagrangianReddyBeamElement3D` | `onedfem.pre_processing.element_library.nonlinear.tosdt.reddy.updated_lagrangian.updated_lagrangian_reddy_3D` | UL Reddy: variational ``K_0`` on the current chord plus beam-column ``K_sigma``. |
| `UpdatedLagrangianTimoshenkoBeamElement3D` | `onedfem.pre_processing.element_library.nonlinear.timoshenko.updated_lagrangian.updated_lagrangian_timoshenko_3D` | UL Timoshenko: selective linear ``K_0`` on ``L_curr`` plus beam-column ``K_sigma``. |

## Modules

### `onedfem.pre_processing.element_library`

[src/onedfem/pre_processing/element_library/__init__.py](../../src/onedfem/pre_processing/element_library/__init__.py)

- **class `ElementFactory(job_results_dir: str) -> None`**
  Batch-builds 1-D element objects and validates their logging infrastructure.
- **class `Element1DBase(*, element_id: int, element_dictionary: Dict[str, Any], grid_dictionary: Dict[str, Any], material_dictionary: Dict[str, Any], section_dictionary: Dict[str, Any], point_load_array: numpy.ndarray, distributed_load_array: numpy.ndarray, job_results_dir: str, dof_per_node: int = 6)`**
  Abstract base class for 1-D structural finite elements.
- **class `LinearEulerBernoulliBeamElement3D(*, element_id: int, element_dictionary: dict, grid_dictionary: dict, section_dictionary: dict, material_dictionary: dict, point_load_array: numpy.ndarray, distributed_load_array: numpy.ndarray, job_results_dir: str, quadrature_order: Optional[int] = None)`**
  2-node straight beam, 6 DOF per node, 12 total (see module docstring for ``U_e`` order).
- **class `LinearTimoshenkoBeamElement3D(*, element_id: int, element_dictionary: dict, grid_dictionary: dict, section_dictionary: dict, material_dictionary: dict, point_load_array: numpy.ndarray, distributed_load_array: numpy.ndarray, job_results_dir: str, quadrature_order: int = 3)`**
  2-node straight beam, 6 DOF per node, 12 total; local ``x`` along the chord.
- **class `LinearLevinsonBeamElement3D(**kwargs)`**
  Levinson third-order shear beam (isotropic or thin-walled sections).
- **class `LinearTrussElement3D(*, element_id: int, element_dictionary: dict, grid_dictionary: dict, section_dictionary: dict, material_dictionary: dict, point_load_array: numpy.ndarray, distributed_load_array: numpy.ndarray, job_results_dir: str, quadrature_order: Optional[int] = None)`**
  2 nodes, 6 DOF/node; transverse displacement follows the element transverse direction (utilities).
- **class `LinearBarElement3D(*, element_id: int, element_dictionary: dict, grid_dictionary: dict, section_dictionary: dict, material_dictionary: dict, point_load_array: numpy.ndarray, distributed_load_array: numpy.ndarray, job_results_dir: str, quadrature_order: Optional[int] = None)`**
  2 nodes, 6 DOF/node; only axial and torsion strains are non-zero in ``B``.

### `onedfem.pre_processing.element_library.beam_load_interpolation`

[src/onedfem/pre_processing/element_library/beam_load_interpolation.py](../../src/onedfem/pre_processing/element_library/beam_load_interpolation.py)

Work-conjugate load interpolation for two-node, six-DOF-per-node beams.

- **`beam_field_interpolation(xi: float, L: float, bending: str = 'hermite') -> np.ndarray`**
  ``N_f`` (6, 12): ``[u, v, w, theta_x, theta_y, theta_z](xi) = N_f U_e``, right-handed.
- **`point_load_vector(xi: float, L: float, load6, bending: str = 'hermite') -> np.ndarray`**
  Nodal loads (12,) of a point load ``load6`` at ``xi``.
- **`distributed_load_vector(xi_g, w_g, q_g, L: float, bending: str = 'hermite') -> np.ndarray`**
  Nodal loads (12,) of line loads ``q_g`` (n_g, 6) at Gauss points ``xi_g``, weights ``w_g``.

### `onedfem.pre_processing.element_library.element_frame`

[src/onedfem/pre_processing/element_library/element_frame.py](../../src/onedfem/pre_processing/element_library/element_frame.py)

Local frames of beam elements: straight elements in any direction, and curved, twisted ones.

- **`rotation_exp(a: npt.ArrayLike) -> np.ndarray`**
  Rotation matrix of the rotation vector ``a`` (Rodrigues).
- **`straight_frame(p0: npt.ArrayLike, p1: npt.ArrayLike, twist_deg: float = 0.0) -> np.ndarray`**
  Columns: the element's ``x``, ``y``, ``z`` in global axes (see the module notes).
- **class `ArcGeometry(R_mid: np.ndarray, curvature: np.ndarray, length: float, chord: float) -> None`**
  An element arc: its frames and length, and the straight sub-chords that approximate it.
- **`arc_geometry(p0: npt.ArrayLike, p1: npt.ArrayLike, curvature: npt.ArrayLike, twist_deg: float = 0.0) -> ArcGeometry`**
  The arc of constant ``curvature`` (element axes) from ``p0`` to ``p1``.
- **`node_dof_rotation(R0: np.ndarray, R1: np.ndarray, n_dof: int, blocks: Optional[list] = None, pairs: Optional[list] = None) -> np.ndarray`**
  ``T`` with ``u_local = T @ u_global`` for an element whose nodes carry frames ``R0``, ``R1``.

### `onedfem.pre_processing.element_library.framed_element`

[src/onedfem/pre_processing/element_library/framed_element.py](../../src/onedfem/pre_processing/element_library/framed_element.py)

Element frames for the linear beams: rotation to global axes, and curved, twisted elements.

- **`setup_frame(element, element_dictionary, grid_dictionary, idx: int, ctor_kwargs: Dict[str, Any]) -> None`**
  Set the element's frame attributes and, on a framed mesh, its nodal point loads.
- **`frame_kernel(name: str, fn)`**
  Wrap a kernel so that it returns global-axis results for a rotated or curved element.
- **`FRAMED_KERNELS`** (constant)
- **`DEFAULT_SUB_ELEMENTS`** (constant)
- **`local_state(element, U_e, xi: float)`**
  ``(element, U_local, xi_local)`` for recovering the section state at ``xi`` of ``element``.
- **`chain(element) -> _Chain`**

### `onedfem.pre_processing.element_library.linear`

[src/onedfem/pre_processing/element_library/linear/__init__.py](../../src/onedfem/pre_processing/element_library/linear/__init__.py)

- **class `LinearEulerBernoulliBeamElement3D(*, element_id: int, element_dictionary: dict, grid_dictionary: dict, section_dictionary: dict, material_dictionary: dict, point_load_array: numpy.ndarray, distributed_load_array: numpy.ndarray, job_results_dir: str, quadrature_order: Optional[int] = None)`**
  2-node straight beam, 6 DOF per node, 12 total (see module docstring for ``U_e`` order).
- **class `LinearTimoshenkoBeamElement3D(*, element_id: int, element_dictionary: dict, grid_dictionary: dict, section_dictionary: dict, material_dictionary: dict, point_load_array: numpy.ndarray, distributed_load_array: numpy.ndarray, job_results_dir: str, quadrature_order: int = 3)`**
  2-node straight beam, 6 DOF per node, 12 total; local ``x`` along the chord.
- **class `LinearBenscoterBeamElement3D(**kwargs)`**
  Timoshenko beam with Benscoter restrained warping for thin-walled anisotropic sections.
- **class `LinearGBTBeamElement3D(**kwargs)`**
  Generalised Beam Theory element for thin-walled anisotropic sections.
- **class `LinearLevinsonBeamElement3D(**kwargs)`**
  Levinson third-order shear beam (isotropic or thin-walled sections).
- **class `LinearReddyBeamElement3D(**kwargs)`**
  Reddy third-order shear beam (isotropic or thin-walled sections).
- **class `LinearTrussElement3D(*, element_id: int, element_dictionary: dict, grid_dictionary: dict, section_dictionary: dict, material_dictionary: dict, point_load_array: numpy.ndarray, distributed_load_array: numpy.ndarray, job_results_dir: str, quadrature_order: Optional[int] = None)`**
  2 nodes, 6 DOF/node; transverse displacement follows the element transverse direction (utilities).
- **class `LinearBarElement3D(*, element_id: int, element_dictionary: dict, grid_dictionary: dict, section_dictionary: dict, material_dictionary: dict, point_load_array: numpy.ndarray, distributed_load_array: numpy.ndarray, job_results_dir: str, quadrature_order: Optional[int] = None)`**
  2 nodes, 6 DOF/node; only axial and torsion strains are non-zero in ``B``.

### `onedfem.pre_processing.element_library.linear.bar`

[src/onedfem/pre_processing/element_library/linear/bar/__init__.py](../../src/onedfem/pre_processing/element_library/linear/bar/__init__.py)

- **class `LinearBarElement3D(*, element_id: int, element_dictionary: dict, grid_dictionary: dict, section_dictionary: dict, material_dictionary: dict, point_load_array: numpy.ndarray, distributed_load_array: numpy.ndarray, job_results_dir: str, quadrature_order: Optional[int] = None)`**
  2 nodes, 6 DOF/node; only axial and torsion strains are non-zero in ``B``.

### `onedfem.pre_processing.element_library.linear.bar.linear_bar_3D`

[src/onedfem/pre_processing/element_library/linear/bar/linear_bar_3D.py](../../src/onedfem/pre_processing/element_library/linear/bar/linear_bar_3D.py)

2-node 3D bar: axial + torsion only (no bending/shear DOF coupling in B).

- **class `LinearBarElement3D(*, element_id: int, element_dictionary: dict, grid_dictionary: dict, section_dictionary: dict, material_dictionary: dict, point_load_array: numpy.ndarray, distributed_load_array: numpy.ndarray, job_results_dir: str, quadrature_order: Optional[int] = None)`**
  2 nodes, 6 DOF/node; only axial and torsion strains are non-zero in ``B``.

### `onedfem.pre_processing.element_library.linear.bar.utilities`

[src/onedfem/pre_processing/element_library/linear/bar/utilities/__init__.py](../../src/onedfem/pre_processing/element_library/linear/bar/utilities/__init__.py)

- **class `ShapeFunctionOperator(element_length: float) -> None`**
  Operator for evaluating 3D Bar element shape functions and their derivatives. Bar has axial and torsion only; shape functions are linear Lagrange in ξ.
- **class `StrainDisplacementOperator(element_length: float, axial: numpy.ndarray) -> None`**
  Strain–displacement matrix **B** for a 2-node 3-D Bar element (axial + torsion).
- **class `MaterialStiffnessOperator(youngs_modulus: float, shear_modulus: float, cross_section_area: float, torsion_constant: float) -> None`**
  Constitutive operator for 3-D Bar elements (axial + torsion only).
- **class `LoadInterpolationOperator(distributed_loads_array: numpy.ndarray, boundary_mode: str = 'error', interpolation_order: str = 'linear', n_gauss_points: int = 3) -> None`**
  Interpolates distributed loads for bar/truss 6-DOF elements. Same API as beam: array shape (N, 9) [x, y, z, Fx, Fy, Fz, Mx, My, Mz], interpolate(x_phys) returns (M, 6) [Fx, Fy, Fz, Mx, My, Mz].
- **`direction_cosines(node_coords: numpy.ndarray) -> numpy.ndarray`**
  Unit vector (cx, cy, cz) along element.
- **`build_L_matrix_4x12(axial: numpy.ndarray) -> numpy.ndarray`**
  Build 4×12 L: local DOF (u1_axial, u2_axial, θx1, θx2) from global 12 DOF. u_local = L @ u_global.

### `onedfem.pre_processing.element_library.linear.beam.fosdt.benscoter`

[src/onedfem/pre_processing/element_library/linear/beam/fosdt/benscoter/__init__.py](../../src/onedfem/pre_processing/element_library/linear/beam/fosdt/benscoter/__init__.py)

- **class `LinearBenscoterBeamElement3D(**kwargs)`**
  Timoshenko beam with Benscoter restrained warping for thin-walled anisotropic sections.

### `onedfem.pre_processing.element_library.linear.beam.fosdt.benscoter.utilities`

[src/onedfem/pre_processing/element_library/linear/beam/fosdt/benscoter/utilities/__init__.py](../../src/onedfem/pre_processing/element_library/linear/beam/fosdt/benscoter/utilities/__init__.py)

- **class `StrainDisplacementOperator(element_length: float, n_dof: int, psi_cols: tuple, parent: object) -> None`**
  Benscoter **B** at Gauss stations; composes Timoshenko 6×12 kinematics.
- **class `MaterialStiffnessOperator(D: np.ndarray) -> None`**
  Wraps the 8×8 Benscoter section stiffness for assembly.

### `onedfem.pre_processing.element_library.linear.beam.fosdt.gbt`

[src/onedfem/pre_processing/element_library/linear/beam/fosdt/gbt/__init__.py](../../src/onedfem/pre_processing/element_library/linear/beam/fosdt/gbt/__init__.py)

- **class `LinearGBTBeamElement3D(**kwargs)`**
  Generalised Beam Theory element for thin-walled anisotropic sections.

### `onedfem.pre_processing.element_library.linear.beam.fosdt.gbt.utilities`

[src/onedfem/pre_processing/element_library/linear/beam/fosdt/gbt/utilities/__init__.py](../../src/onedfem/pre_processing/element_library/linear/beam/fosdt/gbt/utilities/__init__.py)

- **class `StrainDisplacementOperator(element_length: float, n_dof: int, n_modes: int, n_warping: int, parent: object, warping_cols: Callable[[int], tuple], mode_cols: Callable[[int], tuple], reduced_rows: np.ndarray, natural_shape: Callable[[np.ndarray], tuple]) -> None`**
  GBT **B** at a station, with optional B-bar on reduced rows.
- **class `MaterialStiffnessOperator(D: np.ndarray, M: np.ndarray | None = None) -> None`**
  Wraps the GBT stiffness (and optional mass) for assembly.
- **`gbt_modes_of_row(element_dictionary: dict, idx: int, key: str = 'gbt_modes') -> int`**
- **`hermite_x(xi: float, L: float)`**
  Hermite shapes for ``(a_1, a_1', a_2, a_2')`` on ``x = (1 + xi) L / 2``: values, d/dx, d2/dx2.

### `onedfem.pre_processing.element_library.linear.beam.fosdt.timoshenko`

[src/onedfem/pre_processing/element_library/linear/beam/fosdt/timoshenko/__init__.py](../../src/onedfem/pre_processing/element_library/linear/beam/fosdt/timoshenko/__init__.py)

- **class `LinearTimoshenkoBeamElement3D(*, element_id: int, element_dictionary: dict, grid_dictionary: dict, section_dictionary: dict, material_dictionary: dict, point_load_array: numpy.ndarray, distributed_load_array: numpy.ndarray, job_results_dir: str, quadrature_order: int = 3)`**
  2-node straight beam, 6 DOF per node, 12 total; local ``x`` along the chord.

### `onedfem.pre_processing.element_library.linear.beam.fosdt.timoshenko.linear_timoshenko_3D`

[src/onedfem/pre_processing/element_library/linear/beam/fosdt/timoshenko/linear_timoshenko_3D.py](../../src/onedfem/pre_processing/element_library/linear/beam/fosdt/timoshenko/linear_timoshenko_3D.py)

2-node 3D Timoshenko beam (shear-deformable).

- **class `LinearTimoshenkoBeamElement3D(*, element_id: int, element_dictionary: dict, grid_dictionary: dict, section_dictionary: dict, material_dictionary: dict, point_load_array: numpy.ndarray, distributed_load_array: numpy.ndarray, job_results_dir: str, quadrature_order: int = 3)`**
  2-node straight beam, 6 DOF per node, 12 total; local ``x`` along the chord.

### `onedfem.pre_processing.element_library.linear.beam.second_order`

[src/onedfem/pre_processing/element_library/linear/beam/second_order.py](../../src/onedfem/pre_processing/element_library/linear/beam/second_order.py)

Second-order (geometric) stiffness of the linear beam elements from their section forces.

- **`geometric_section_matrix(c: np.ndarray) -> np.ndarray`**
  ``S`` (4, 4) on ``[v', w', phi', phi]`` from the section terms ``c`` (7,).
- **class `GeometricStiffnessMixin()`**
  ``linear_geometric_stiffness_matrix`` for an element that defines ``_geometric_state``.

### `onedfem.pre_processing.element_library.linear.beam.tosdt.levinson`

[src/onedfem/pre_processing/element_library/linear/beam/tosdt/levinson/__init__.py](../../src/onedfem/pre_processing/element_library/linear/beam/tosdt/levinson/__init__.py)

- **class `LinearLevinsonBeamElement3D(**kwargs)`**
  Levinson third-order shear beam (isotropic or thin-walled sections).

### `onedfem.pre_processing.element_library.linear.beam.tosdt.reddy`

[src/onedfem/pre_processing/element_library/linear/beam/tosdt/reddy/__init__.py](../../src/onedfem/pre_processing/element_library/linear/beam/tosdt/reddy/__init__.py)

- **class `LinearReddyBeamElement3D(**kwargs)`**
  Reddy third-order shear beam (isotropic or thin-walled sections).

### `onedfem.pre_processing.element_library.linear.beam.tosdt.utilities`

[src/onedfem/pre_processing/element_library/linear/beam/tosdt/utilities/__init__.py](../../src/onedfem/pre_processing/element_library/linear/beam/tosdt/utilities/__init__.py)

- **class `StrainDisplacementOperator(element_length: float, rate: np.ndarray) -> None`**
  StrainDisplacementOperator(element_length: 'float', rate: 'np.ndarray')
- **class `MaterialStiffnessOperator(D: np.ndarray) -> None`**
  MaterialStiffnessOperator(D: 'np.ndarray')
- **class `LoadInterpolationOperator(distributed_loads_array: numpy.ndarray, boundary_mode: str = 'error', interpolation_order: str = 'linear', n_gauss_points: int = 3) -> None`**
  Interpolates distributed loads for 6-DOF finite elements with: - Boundary-aware extrapolation control - Automatic monotonic sorting of position data - Adaptive fallback of interpolation scheme based on Gauss point resolution
- **class `ShapeFunctionOperator(element_length: float) -> None`**
  Field interpolation **N_q** (6, 16) at a natural coordinate.
- **`hermite_x(xi: float, L: float)`**
  Cubic Hermite on ``x = (1 + xi) L / 2``: values, d/dx, d2/dx2 for ``(f1, f1', f2, f2')``.

### `onedfem.pre_processing.element_library.linear.beam.tosdt.utilities.interpolate_loads`

[src/onedfem/pre_processing/element_library/linear/beam/tosdt/utilities/interpolate_loads.py](../../src/onedfem/pre_processing/element_library/linear/beam/tosdt/utilities/interpolate_loads.py)

Distributed-load interpolation for third-order beams — delegates to Timoshenko.

- **class `LoadInterpolationOperator(distributed_loads_array: numpy.ndarray, boundary_mode: str = 'error', interpolation_order: str = 'linear', n_gauss_points: int = 3) -> None`**
  Interpolates distributed loads for 6-DOF finite elements with: - Boundary-aware extrapolation control - Automatic monotonic sorting of position data - Adaptive fallback of interpolation scheme based on Gauss point resolution

### `onedfem.pre_processing.element_library.linear.beam.zosdt.euler_bernoulli`

[src/onedfem/pre_processing/element_library/linear/beam/zosdt/euler_bernoulli/__init__.py](../../src/onedfem/pre_processing/element_library/linear/beam/zosdt/euler_bernoulli/__init__.py)

- **class `LinearEulerBernoulliBeamElement3D(*, element_id: int, element_dictionary: dict, grid_dictionary: dict, section_dictionary: dict, material_dictionary: dict, point_load_array: numpy.ndarray, distributed_load_array: numpy.ndarray, job_results_dir: str, quadrature_order: Optional[int] = None)`**
  2-node straight beam, 6 DOF per node, 12 total (see module docstring for ``U_e`` order).

### `onedfem.pre_processing.element_library.linear.beam.zosdt.euler_bernoulli.linear_euler_bernoulli_3D`

[src/onedfem/pre_processing/element_library/linear/beam/zosdt/euler_bernoulli/linear_euler_bernoulli_3D.py](../../src/onedfem/pre_processing/element_library/linear/beam/zosdt/euler_bernoulli/linear_euler_bernoulli_3D.py)

2-node 3D Euler–Bernoulli beam.

- **class `LinearEulerBernoulliBeamElement3D(*, element_id: int, element_dictionary: dict, grid_dictionary: dict, section_dictionary: dict, material_dictionary: dict, point_load_array: numpy.ndarray, distributed_load_array: numpy.ndarray, job_results_dir: str, quadrature_order: Optional[int] = None)`**
  2-node straight beam, 6 DOF per node, 12 total (see module docstring for ``U_e`` order).

### `onedfem.pre_processing.element_library.linear.beam.zosdt.euler_bernoulli_with_warp.utilities`

[src/onedfem/pre_processing/element_library/linear/beam/zosdt/euler_bernoulli_with_warp/utilities/__init__.py](../../src/onedfem/pre_processing/element_library/linear/beam/zosdt/euler_bernoulli_with_warp/utilities/__init__.py)

- **`N_DOF`** (constant)
- **`N_STANDARD_DOF`** (constant)
- **`N_STRAIN`** (constant)
- **`N_WARPING_DOF`** (constant)
- **class `WarpingMaterialStiffnessOperator(base_material_operator: MaterialStiffnessOperator, youngs_modulus: float, warping_gamma: float) -> None`**
  Constitutive tensor ``D ∈ ℝ^{7×7}`` for 3-D Euler–Bernoulli beam elements with Vlasov warping.
- **class `WarpingStrainDisplacementOperator(element_length: float, base_strain_operator: StrainDisplacementOperator) -> None`**
  Strain–displacement tensor ``B ∈ ℝ^{7×14}`` for a 2-node 3-D Euler–Bernoulli beam with Vlasov warping.
- **`extend_natural_shape_to_warping(N12: numpy.ndarray, xi: float) -> numpy.ndarray`**
  Extend standard (12, 6) ``N`` at one Gauss point to (14, 6).

### `onedfem.pre_processing.element_library.linear.truss`

[src/onedfem/pre_processing/element_library/linear/truss/__init__.py](../../src/onedfem/pre_processing/element_library/linear/truss/__init__.py)

- **class `LinearTrussElement3D(*, element_id: int, element_dictionary: dict, grid_dictionary: dict, section_dictionary: dict, material_dictionary: dict, point_load_array: numpy.ndarray, distributed_load_array: numpy.ndarray, job_results_dir: str, quadrature_order: Optional[int] = None)`**
  2 nodes, 6 DOF/node; transverse displacement follows the element transverse direction (utilities).

### `onedfem.pre_processing.element_library.linear.truss.linear_truss_3D`

[src/onedfem/pre_processing/element_library/linear/truss/linear_truss_3D.py](../../src/onedfem/pre_processing/element_library/linear/truss/linear_truss_3D.py)

2-node 3D truss: axial + one transverse shear direction + torsion (six DOF/node for layout compatibility).

- **class `LinearTrussElement3D(*, element_id: int, element_dictionary: dict, grid_dictionary: dict, section_dictionary: dict, material_dictionary: dict, point_load_array: numpy.ndarray, distributed_load_array: numpy.ndarray, job_results_dir: str, quadrature_order: Optional[int] = None)`**
  2 nodes, 6 DOF/node; transverse displacement follows the element transverse direction (utilities).

### `onedfem.pre_processing.element_library.linear.truss.utilities`

[src/onedfem/pre_processing/element_library/linear/truss/utilities/__init__.py](../../src/onedfem/pre_processing/element_library/linear/truss/utilities/__init__.py)

- **class `ShapeFunctionOperator(element_length: float) -> None`**
  Operator for evaluating 3D Truss element shape functions and their derivatives. Truss has axial, transverse, and torsion; all use linear Lagrange in ξ.
- **class `StrainDisplacementOperator(element_length: float, axial: numpy.ndarray, transverse: numpy.ndarray) -> None`**
  Strain–displacement matrix **B** for a 2-node 3-D Truss element (axial + transverse + torsion).
- **class `MaterialStiffnessOperator(youngs_modulus: float, shear_modulus: float, cross_section_area: float, torsion_constant: float, shear_correction_factor: float = 0.8333333333333334) -> None`**
  Constitutive operator for 3-D Truss elements (axial + transverse shear + torsion).
- **class `LoadInterpolationOperator(distributed_loads_array: numpy.ndarray, boundary_mode: str = 'error', interpolation_order: str = 'linear', n_gauss_points: int = 3) -> None`**
  Interpolates distributed loads for truss 6-DOF elements. Same API as beam: array shape (N, 9) [x, y, z, Fx, Fy, Fz, Mx, My, Mz], interpolate(x_phys) returns (M, 6) [Fx, Fy, Fz, Mx, My, Mz].
- **`direction_cosines_and_transverse(node_coords: numpy.ndarray) -> tuple[numpy.ndarray, numpy.ndarray]`**
  Compute axial direction cosines and one transverse direction (local y) for the element.
- **`build_L_matrix_6x12(axial: numpy.ndarray, transverse: numpy.ndarray) -> numpy.ndarray`**
  Build the 6×12 matrix L that maps global 12 DOF to local 6 DOF.

### `onedfem.pre_processing.element_library.nonlinear`

[src/onedfem/pre_processing/element_library/nonlinear/__init__.py](../../src/onedfem/pre_processing/element_library/nonlinear/__init__.py)

- **class `NonlinearEulerBernoulliBeamElement3D(**kwargs)`**
  TL EB: Green–Lagrange strain, ``S = D @ E``, material and geometric tangents from Gauss sums.
- **class `NonlinearTimoshenkoBeamElement3D(**kwargs)`**
  TL Timoshenko: selective ``K_0`` plus ``K_delta`` and Hessian ``K_sigma``.
- **class `NonlinearLevinsonBeamElement3D(**kwargs)`**
  TL Levinson: Petrov–Galerkin test on six resultants plus GL increments.
- **class `NonlinearReddyBeamElement3D(**kwargs)`**
  TL Reddy: variational ``K = int B8.T D B8`` plus GL ``K_delta`` and Hessian ``K_sigma``.
- **class `UpdatedLagrangianEulerBernoulliBeamElement3D(**kwargs)`**
  UL Euler–Bernoulli: linear ``B``, ``D`` on the current chord plus beam-column ``K_sigma``.
- **class `UpdatedLagrangianTimoshenkoBeamElement3D(**kwargs)`**
  UL Timoshenko: selective linear ``K_0`` on ``L_curr`` plus beam-column ``K_sigma``.
- **class `UpdatedLagrangianLevinsonBeamElement3D(**kwargs)`**
  UL Levinson: Petrov–Galerkin ``K_0`` on the current chord plus beam-column ``K_sigma``.
- **class `UpdatedLagrangianReddyBeamElement3D(**kwargs)`**
  UL Reddy: variational ``K_0`` on the current chord plus beam-column ``K_sigma``.
- **class `CorotationalBeamElement3D(*, element_id: int, element_dictionary: dict, grid_dictionary: dict, section_dictionary: dict, material_dictionary: dict, point_load_array: np.ndarray, distributed_load_array: np.ndarray, job_results_dir: str, quadrature_order: int | None = None, kernel: str = 'timoshenko', tangent_stiffness_mode: "Literal['finite_difference', 'elastic_material']" = 'finite_difference')`**
  Co-rotational 3D beam (``kernel`` = ``timoshenko`` or ``euler_bernoulli``).
- **class `CorotationalEulerBernoulliBeamElement3D(**kwargs)`**
  Corotational Euler-Bernoulli beam (shear-rigid core, classical section stiffness).
- **class `CorotationalTimoshenkoBeamElement3D(**kwargs)`**
  Corotational Timoshenko beam on the coupled Saint-Venant section stiffness.
- **class `CorotationalBenscoterBeamElement3D(**kwargs)`**
  Corotational Benscoter beam: restrained warping and large twist on thin-walled sections.
- **class `CorotationalGBTBeamElement3D(**kwargs)`**
  Corotational GBT beam: large rotations of the section, with its warping and distortional modes linear in the element's own frame.
- **class `CorotationalReddyBeamElement3D(**kwargs)`**
  Corotational Reddy beam: third-order shear (shear lag of the slaved warping) at large rotation.
- **class `CorotationalLevinsonBeamElement3D(**kwargs)`**
  Corotational Levinson beam (non-variational third-order shear; unsymmetric tangent).
- **class `GeometricallyExactShearDeformableBeam3D(*, element_id: int, element_dictionary: dict, grid_dictionary: dict, section_dictionary: dict, material_dictionary: dict, point_load_array: np.ndarray, distributed_load_array: np.ndarray, job_results_dir: str, quadrature_order: int | None = 3, gesdb_tl_fallback: bool = False, gesdb_kernel: str = 'tl_locked')`**
  Factory-facing **GeometricallyExactShearDeformableBeam3D** element type.
- **class `GEBTUnshearableBeamElement3D(**kwargs)`**
  Unshearable GEBT: Kirchhoff stretch and curvature, linear EB section law.

### `onedfem.pre_processing.element_library.nonlinear.euler_bernoulli`

[src/onedfem/pre_processing/element_library/nonlinear/euler_bernoulli/__init__.py](../../src/onedfem/pre_processing/element_library/nonlinear/euler_bernoulli/__init__.py)

- **class `NonlinearEulerBernoulliBeamElement3D(**kwargs)`**
  TL EB: Green–Lagrange strain, ``S = D @ E``, material and geometric tangents from Gauss sums.
- **class `UpdatedLagrangianEulerBernoulliBeamElement3D(**kwargs)`**
  UL Euler–Bernoulli: linear ``B``, ``D`` on the current chord plus beam-column ``K_sigma``.

### `onedfem.pre_processing.element_library.nonlinear.euler_bernoulli.total_lagrangian`

[src/onedfem/pre_processing/element_library/nonlinear/euler_bernoulli/total_lagrangian/__init__.py](../../src/onedfem/pre_processing/element_library/nonlinear/euler_bernoulli/total_lagrangian/__init__.py)

- **class `NonlinearEulerBernoulliBeamElement3D(**kwargs)`**
  TL EB: Green–Lagrange strain, ``S = D @ E``, material and geometric tangents from Gauss sums.

### `onedfem.pre_processing.element_library.nonlinear.euler_bernoulli.total_lagrangian.utilities`

[src/onedfem/pre_processing/element_library/nonlinear/euler_bernoulli/total_lagrangian/utilities/__init__.py](../../src/onedfem/pre_processing/element_library/nonlinear/euler_bernoulli/total_lagrangian/utilities/__init__.py)

- **class `GreenLagrangeStrainOperator(element_length: float, include_shear: bool = False) -> None`**
  Green–Lagrange strain measures for a 2-node 3D Euler–Bernoulli beam (Total Lagrangian).
- **class `StressResultantOperator() -> None`**
  Section forces (**2PK-type resultants**) from Green–Lagrange strain \(\mathbf{E}\) in the **reference** configuration.
- **class `GeometricStiffnessOperator(element_length: float) -> None`**
  Assemble \(\mathbf{K}_\sigma\) from section forces and shape-function derivatives (**TL nonlinear**).

### `onedfem.pre_processing.element_library.nonlinear.euler_bernoulli.updated_lagrangian`

[src/onedfem/pre_processing/element_library/nonlinear/euler_bernoulli/updated_lagrangian/__init__.py](../../src/onedfem/pre_processing/element_library/nonlinear/euler_bernoulli/updated_lagrangian/__init__.py)

- **class `UpdatedLagrangianEulerBernoulliBeamElement3D(**kwargs)`**
  UL Euler–Bernoulli: linear ``B``, ``D`` on the current chord plus beam-column ``K_sigma``.

### `onedfem.pre_processing.element_library.nonlinear.euler_bernoulli.updated_lagrangian.utilities`

[src/onedfem/pre_processing/element_library/nonlinear/euler_bernoulli/updated_lagrangian/utilities/__init__.py](../../src/onedfem/pre_processing/element_library/nonlinear/euler_bernoulli/updated_lagrangian/utilities/__init__.py)

- **class `ConfigurationUpdateOperator() -> None`**
  Push nodal translations from ``U_e`` into the current chord endpoints.
- **class `CurrentChordFrameOperator() -> None`**
  Map increment DOFs between global axes and the current chord.
- **class `UpdatedLagrangianChordMixin()`**
  Local ``K_mat``/``K_σ`` on the current chord, rotated to global increment DOFs.
- **class `CurrentConfigurationStrainOperator() -> None`**
  ``eps = B @ U_e`` with linear ``B`` evaluated on the current length.
- **class `GeometricStiffnessOperator(element_length: float) -> None`**
  Assemble \(\mathbf{K}_\sigma\) from section forces and shape-function derivatives (**TL nonlinear**).
- **class `StressResultantOperator() -> None`**
  StressResultantOperator()

### `onedfem.pre_processing.element_library.nonlinear.euler_bernoulli.updated_lagrangian.utilities.geometric_stiffness`

[src/onedfem/pre_processing/element_library/nonlinear/euler_bernoulli/updated_lagrangian/utilities/geometric_stiffness.py](../../src/onedfem/pre_processing/element_library/nonlinear/euler_bernoulli/updated_lagrangian/utilities/geometric_stiffness.py)

- **class `GeometricStiffnessOperator(element_length: float) -> None`**
  Assemble \(\mathbf{K}_\sigma\) from section forces and shape-function derivatives (**TL nonlinear**).

### `onedfem.pre_processing.element_library.nonlinear.large_rotations`

[src/onedfem/pre_processing/element_library/nonlinear/large_rotations/__init__.py](../../src/onedfem/pre_processing/element_library/nonlinear/large_rotations/__init__.py)

Finite-rotation beam kinematics (co-rotational frame, GEBT-class registrations).

- **class `CorotationalBeamElement3D(*, element_id: int, element_dictionary: dict, grid_dictionary: dict, section_dictionary: dict, material_dictionary: dict, point_load_array: np.ndarray, distributed_load_array: np.ndarray, job_results_dir: str, quadrature_order: int | None = None, kernel: str = 'timoshenko', tangent_stiffness_mode: "Literal['finite_difference', 'elastic_material']" = 'finite_difference')`**
  Co-rotational 3D beam (``kernel`` = ``timoshenko`` or ``euler_bernoulli``).
- **class `GeometricallyExactShearDeformableBeam3D(*, element_id: int, element_dictionary: dict, grid_dictionary: dict, section_dictionary: dict, material_dictionary: dict, point_load_array: np.ndarray, distributed_load_array: np.ndarray, job_results_dir: str, quadrature_order: int | None = 3, gesdb_tl_fallback: bool = False, gesdb_kernel: str = 'tl_locked')`**
  Factory-facing **GeometricallyExactShearDeformableBeam3D** element type.
- **class `GEBTUnshearableBeamElement3D(**kwargs)`**
  Unshearable GEBT: Kirchhoff stretch and curvature, linear EB section law.

### `onedfem.pre_processing.element_library.nonlinear.large_rotations.corotational`

[src/onedfem/pre_processing/element_library/nonlinear/large_rotations/corotational/__init__.py](../../src/onedfem/pre_processing/element_library/nonlinear/large_rotations/corotational/__init__.py)

Co-rotational beam (12 DOF).

- **class `CorotationalBeamElement3D(*, element_id: int, element_dictionary: dict, grid_dictionary: dict, section_dictionary: dict, material_dictionary: dict, point_load_array: np.ndarray, distributed_load_array: np.ndarray, job_results_dir: str, quadrature_order: int | None = None, kernel: str = 'timoshenko', tangent_stiffness_mode: "Literal['finite_difference', 'elastic_material']" = 'finite_difference')`**
  Co-rotational 3D beam (``kernel`` = ``timoshenko`` or ``euler_bernoulli``).
- **class `CorotationalEulerBernoulliBeamElement3D(**kwargs)`**
  Corotational Euler-Bernoulli beam (shear-rigid core, classical section stiffness).
- **class `CorotationalTimoshenkoBeamElement3D(**kwargs)`**
  Corotational Timoshenko beam on the coupled Saint-Venant section stiffness.
- **class `CorotationalBenscoterBeamElement3D(**kwargs)`**
  Corotational Benscoter beam: restrained warping and large twist on thin-walled sections.
- **class `CorotationalGBTBeamElement3D(**kwargs)`**
  Corotational GBT beam: large rotations of the section, with its warping and distortional modes linear in the element's own frame.
- **class `CorotationalReddyBeamElement3D(**kwargs)`**
  Corotational Reddy beam: third-order shear (shear lag of the slaved warping) at large rotation.
- **class `CorotationalLevinsonBeamElement3D(**kwargs)`**
  Corotational Levinson beam (non-variational third-order shear; unsymmetric tangent).

### `onedfem.pre_processing.element_library.nonlinear.large_rotations.corotational.corotational_section`

[src/onedfem/pre_processing/element_library/nonlinear/large_rotations/corotational/corotational_section.py](../../src/onedfem/pre_processing/element_library/nonlinear/large_rotations/corotational/corotational_section.py)

Corotational beam elements on any SectionAnalysis, with third- and fourth-order local energy.

- **class `CorotationalMixin(**kwargs)`**
  Corotational large-displacement behaviour on top of a linear element core.
- **class `CorotationalEulerBernoulliBeamElement3D(**kwargs)`**
  Corotational Euler-Bernoulli beam (shear-rigid core, classical section stiffness).
- **class `CorotationalTimoshenkoBeamElement3D(**kwargs)`**
  Corotational Timoshenko beam on the coupled Saint-Venant section stiffness.
- **class `CorotationalBenscoterBeamElement3D(**kwargs)`**
  Corotational Benscoter beam: restrained warping and large twist on thin-walled sections.
- **class `CorotationalGBTBeamElement3D(**kwargs)`**
  Corotational GBT beam: large rotations of the section, with its warping and distortional modes linear in the element's own frame.
- **`brazier_energy_derivatives(T: np.ndarray, Q: np.ndarray, z: np.ndarray) -> Tuple[np.ndarray, np.ndarray]`**
  Gradient and Hessian of ``W = a^T (kappa . T E) + 1/2 a a kappa kappa . Q`` in ``z = [a, kappa, E]``.
- **class `CorotationalReddyBeamElement3D(**kwargs)`**
  Corotational Reddy beam: third-order shear (shear lag of the slaved warping) at large rotation.
- **class `CorotationalLevinsonBeamElement3D(**kwargs)`**
  Corotational Levinson beam (non-variational third-order shear; unsymmetric tangent).

### `onedfem.pre_processing.element_library.nonlinear.large_rotations.corotational.utilities`

[src/onedfem/pre_processing/element_library/nonlinear/large_rotations/corotational/utilities/__init__.py](../../src/onedfem/pre_processing/element_library/nonlinear/large_rotations/corotational/utilities/__init__.py)

- **`block_twelve(R: np.ndarray) -> np.ndarray`**
  Map global nodal vectors to corotated axes: four R.T blocks (trans/rot at each node).
- **`orthogonal_basis_from_chord(d: np.ndarray) -> np.ndarray`**
  Rotation matrix R (3, 3) with columns e1, e2, e3; e1 aligns with chord ``d``.
- **`reference_frame(d: np.ndarray) -> np.ndarray`**
  Undeformed element frame: ``e1`` along the chord, ``e2`` towards global ``y`` (or ``z``).
- **`K_local(elem, L: float) -> np.ndarray`**
- **`cubic_quartic_force_and_tangent(p: np.ndarray, K_core: np.ndarray, C: np.ndarray, axial_pts: Optional[np.ndarray], G_gauss: Sequence[np.ndarray], gauss_w: Iterable[float], l0: float, f_E0: Optional[np.ndarray] = None, K_sigma0: Optional[np.ndarray] = None) -> Tuple[np.ndarray, np.ndarray]`**
  Local force ``f`` and Hessian ``H`` of ``½ pᵀ K p + ½ pᵀ K_σ(p) p + Π₄``.
- **`battini_p_batch(U: np.ndarray, X1: np.ndarray, X2: np.ndarray, E0: np.ndarray, l0: float, extras_fn: Callable[[np.ndarray, np.ndarray], np.ndarray], element_id: int) -> np.ndarray`**
  Local deformational DOFs ``p`` (m, n) from global states ``U`` (m, n).
- **`vector_norm(a: np.ndarray) -> np.ndarray`**
  sqrt(a . a) along the last axis (complex-step safe).
- **`exp_batch(th: np.ndarray) -> np.ndarray`**
  Rotation matrices (m, 3, 3) of rotation vectors (m, 3); complex-step safe.
- **`log_batch(R: np.ndarray) -> np.ndarray`**
  Rotation vectors (m, 3) of rotations (m, 3, 3) below pi/2 in angle.
- **`rotation_exp(theta: np.ndarray) -> np.ndarray`**
  Rotation matrix of the rotation vector ``theta`` (Rodrigues).
- **`rotation_log(R: np.ndarray) -> np.ndarray`**
  Rotation vector of ``R`` (through the quaternion, robust up to pi).
- **`skew(a: np.ndarray) -> np.ndarray`**

### `onedfem.pre_processing.element_library.nonlinear.large_rotations.gebt_unshearable`

[src/onedfem/pre_processing/element_library/nonlinear/large_rotations/gebt_unshearable/__init__.py](../../src/onedfem/pre_processing/element_library/nonlinear/large_rotations/gebt_unshearable/__init__.py)

GEBT unshearable (Euler-Bernoulli) beam (12 DOF).

- **class `GEBTUnshearableBeamElement3D(**kwargs)`**
  Unshearable GEBT: Kirchhoff stretch and curvature, linear EB section law.

### `onedfem.pre_processing.element_library.nonlinear.large_rotations.gebt_unshearable.utilities`

[src/onedfem/pre_processing/element_library/nonlinear/large_rotations/gebt_unshearable/utilities/__init__.py](../../src/onedfem/pre_processing/element_library/nonlinear/large_rotations/gebt_unshearable/utilities/__init__.py)

- **`hessian_contraction(U_e: np.ndarray, f: np.ndarray, X1, X2, E0, l0, eid) -> np.ndarray`**
- **`jacobian_p(U_e: np.ndarray, X1, X2, E0, l0, eid) -> np.ndarray`**
- **`local_force(p: np.ndarray, K_core: np.ndarray, f_E0: np.ndarray) -> Tuple[np.ndarray, np.ndarray]`**
- **`local_p(U: np.ndarray, X1: np.ndarray, X2: np.ndarray, E0: np.ndarray, l0: float, eid: int) -> np.ndarray`**

### `onedfem.pre_processing.element_library.nonlinear.large_rotations.geometrically_exact_shear_deformable_beam`

[src/onedfem/pre_processing/element_library/nonlinear/large_rotations/geometrically_exact_shear_deformable_beam/__init__.py](../../src/onedfem/pre_processing/element_library/nonlinear/large_rotations/geometrically_exact_shear_deformable_beam/__init__.py)

- **class `GeometricallyExactShearDeformableBeam3D(*, element_id: int, element_dictionary: dict, grid_dictionary: dict, section_dictionary: dict, material_dictionary: dict, point_load_array: np.ndarray, distributed_load_array: np.ndarray, job_results_dir: str, quadrature_order: int | None = 3, gesdb_tl_fallback: bool = False, gesdb_kernel: str = 'tl_locked')`**
  Factory-facing **GeometricallyExactShearDeformableBeam3D** element type.

### `onedfem.pre_processing.element_library.nonlinear.large_rotations.geometrically_exact_shear_deformable_beam.utilities`

[src/onedfem/pre_processing/element_library/nonlinear/large_rotations/geometrically_exact_shear_deformable_beam/utilities/__init__.py](../../src/onedfem/pre_processing/element_library/nonlinear/large_rotations/geometrically_exact_shear_deformable_beam/utilities/__init__.py)

- **`native_internal_force(elem, U_e: np.ndarray) -> np.ndarray`**
  ``B_engᵀ S`` Gauss sum on the 12-DOF native kernel.
- **`native_tangent_stiffness(elem, U_e: np.ndarray) -> np.ndarray`**
  ``K_mat + K_sigma`` for the native kernel (no TL ``K_delta`` split).

### `onedfem.pre_processing.element_library.nonlinear.timoshenko`

[src/onedfem/pre_processing/element_library/nonlinear/timoshenko/__init__.py](../../src/onedfem/pre_processing/element_library/nonlinear/timoshenko/__init__.py)

- **class `NonlinearTimoshenkoBeamElement3D(**kwargs)`**
  TL Timoshenko: selective ``K_0`` plus ``K_delta`` and Hessian ``K_sigma``.
- **class `UpdatedLagrangianTimoshenkoBeamElement3D(**kwargs)`**
  UL Timoshenko: selective linear ``K_0`` on ``L_curr`` plus beam-column ``K_sigma``.

### `onedfem.pre_processing.element_library.nonlinear.timoshenko.total_lagrangian`

[src/onedfem/pre_processing/element_library/nonlinear/timoshenko/total_lagrangian/__init__.py](../../src/onedfem/pre_processing/element_library/nonlinear/timoshenko/total_lagrangian/__init__.py)

- **class `NonlinearTimoshenkoBeamElement3D(**kwargs)`**
  TL Timoshenko: selective ``K_0`` plus ``K_delta`` and Hessian ``K_sigma``.

### `onedfem.pre_processing.element_library.nonlinear.timoshenko.total_lagrangian.utilities`

[src/onedfem/pre_processing/element_library/nonlinear/timoshenko/total_lagrangian/utilities/__init__.py](../../src/onedfem/pre_processing/element_library/nonlinear/timoshenko/total_lagrangian/utilities/__init__.py)

- **class `GreenLagrangeStrainOperator(element_length: float, include_shear: bool = False) -> None`**
  Green–Lagrange strain measures for a 2-node 3D beam in Total Lagrangian formulation.
- **class `StressResultantOperator() -> None`**
  Section forces (**2PK-type resultants**) from Green–Lagrange strain \(\mathbf{E}\) in the **reference** configuration.
- **class `GeometricStiffnessOperator(element_length: float) -> None`**
  Assemble \(\mathbf{K}_\sigma\) from section forces and shape-function derivatives (**TL nonlinear**).

### `onedfem.pre_processing.element_library.nonlinear.timoshenko.total_lagrangian.utilities.geometric_stiffness`

[src/onedfem/pre_processing/element_library/nonlinear/timoshenko/total_lagrangian/utilities/geometric_stiffness.py](../../src/onedfem/pre_processing/element_library/nonlinear/timoshenko/total_lagrangian/utilities/geometric_stiffness.py)

\(\mathbf{K}_\sigma\) for nonlinear Timoshenko — beam-column Gauss form (**re-export of EB operator**).

- **class `GeometricStiffnessOperator(element_length: float) -> None`**
  Assemble \(\mathbf{K}_\sigma\) from section forces and shape-function derivatives (**TL nonlinear**).

### `onedfem.pre_processing.element_library.nonlinear.timoshenko.updated_lagrangian`

[src/onedfem/pre_processing/element_library/nonlinear/timoshenko/updated_lagrangian/__init__.py](../../src/onedfem/pre_processing/element_library/nonlinear/timoshenko/updated_lagrangian/__init__.py)

- **class `UpdatedLagrangianTimoshenkoBeamElement3D(**kwargs)`**
  UL Timoshenko: selective linear ``K_0`` on ``L_curr`` plus beam-column ``K_sigma``.

### `onedfem.pre_processing.element_library.nonlinear.timoshenko.updated_lagrangian.utilities`

[src/onedfem/pre_processing/element_library/nonlinear/timoshenko/updated_lagrangian/utilities/__init__.py](../../src/onedfem/pre_processing/element_library/nonlinear/timoshenko/updated_lagrangian/utilities/__init__.py)

- **class `ConfigurationUpdateOperator() -> None`**
  Push nodal translations from ``U_e`` into the current chord endpoints.
- **class `CurrentConfigurationStrainOperator() -> None`**
  ``eps = B @ U_e`` with linear ``B`` evaluated on the current length.
- **class `GeometricStiffnessOperator(element_length: float) -> None`**
  Assemble \(\mathbf{K}_\sigma\) from section forces and shape-function derivatives (**TL nonlinear**).
- **class `StressResultantOperator() -> None`**
  StressResultantOperator()

### `onedfem.pre_processing.element_library.nonlinear.timoshenko.updated_lagrangian.utilities.config_update`

[src/onedfem/pre_processing/element_library/nonlinear/timoshenko/updated_lagrangian/utilities/config_update.py](../../src/onedfem/pre_processing/element_library/nonlinear/timoshenko/updated_lagrangian/utilities/config_update.py)

UL Timoshenko configuration update (same chord-push as UL Euler–Bernoulli).

- **class `ConfigurationUpdateOperator() -> None`**
  Push nodal translations from ``U_e`` into the current chord endpoints.

### `onedfem.pre_processing.element_library.nonlinear.timoshenko.updated_lagrangian.utilities.current_strain`

[src/onedfem/pre_processing/element_library/nonlinear/timoshenko/updated_lagrangian/utilities/current_strain.py](../../src/onedfem/pre_processing/element_library/nonlinear/timoshenko/updated_lagrangian/utilities/current_strain.py)

UL Timoshenko current-chord strain ``eps = B U_e``.

- **class `CurrentConfigurationStrainOperator() -> None`**
  ``eps = B @ U_e`` with linear ``B`` evaluated on the current length.

### `onedfem.pre_processing.element_library.nonlinear.timoshenko.updated_lagrangian.utilities.geometric_stiffness`

[src/onedfem/pre_processing/element_library/nonlinear/timoshenko/updated_lagrangian/utilities/geometric_stiffness.py](../../src/onedfem/pre_processing/element_library/nonlinear/timoshenko/updated_lagrangian/utilities/geometric_stiffness.py)

UL Timoshenko current-config beam-column ``K_σ``.

- **class `GeometricStiffnessOperator(element_length: float) -> None`**
  Assemble \(\mathbf{K}_\sigma\) from section forces and shape-function derivatives (**TL nonlinear**).

### `onedfem.pre_processing.element_library.nonlinear.timoshenko.updated_lagrangian.utilities.stress_resultant`

[src/onedfem/pre_processing/element_library/nonlinear/timoshenko/updated_lagrangian/utilities/stress_resultant.py](../../src/onedfem/pre_processing/element_library/nonlinear/timoshenko/updated_lagrangian/utilities/stress_resultant.py)

UL Timoshenko current-config stress resultants ``S = D eps``.

- **class `StressResultantOperator() -> None`**
  StressResultantOperator()

### `onedfem.pre_processing.element_library.nonlinear.tosdt`

[src/onedfem/pre_processing/element_library/nonlinear/tosdt/__init__.py](../../src/onedfem/pre_processing/element_library/nonlinear/tosdt/__init__.py)

- **class `NonlinearLevinsonBeamElement3D(**kwargs)`**
  TL Levinson: Petrov–Galerkin test on six resultants plus GL increments.
- **class `NonlinearReddyBeamElement3D(**kwargs)`**
  TL Reddy: variational ``K = int B8.T D B8`` plus GL ``K_delta`` and Hessian ``K_sigma``.
- **class `UpdatedLagrangianLevinsonBeamElement3D(**kwargs)`**
  UL Levinson: Petrov–Galerkin ``K_0`` on the current chord plus beam-column ``K_sigma``.
- **class `UpdatedLagrangianReddyBeamElement3D(**kwargs)`**
  UL Reddy: variational ``K_0`` on the current chord plus beam-column ``K_sigma``.

### `onedfem.pre_processing.element_library.nonlinear.tosdt.levinson`

[src/onedfem/pre_processing/element_library/nonlinear/tosdt/levinson/__init__.py](../../src/onedfem/pre_processing/element_library/nonlinear/tosdt/levinson/__init__.py)

Nonlinear (Total Lagrangian) Levinson 3rd-order beam (16 DOF).

- **class `NonlinearLevinsonBeamElement3D(**kwargs)`**
  TL Levinson: Petrov–Galerkin test on six resultants plus GL increments.

### `onedfem.pre_processing.element_library.nonlinear.tosdt.levinson.updated_lagrangian`

[src/onedfem/pre_processing/element_library/nonlinear/tosdt/levinson/updated_lagrangian/__init__.py](../../src/onedfem/pre_processing/element_library/nonlinear/tosdt/levinson/updated_lagrangian/__init__.py)

- **class `UpdatedLagrangianLevinsonBeamElement3D(**kwargs)`**
  UL Levinson: Petrov–Galerkin ``K_0`` on the current chord plus beam-column ``K_sigma``.

### `onedfem.pre_processing.element_library.nonlinear.tosdt.levinson.utilities`

[src/onedfem/pre_processing/element_library/nonlinear/tosdt/levinson/utilities/__init__.py](../../src/onedfem/pre_processing/element_library/nonlinear/tosdt/levinson/utilities/__init__.py)

TL operators for nonlinear Levinson — shared ``nonlinear.tosdt.utilities`` trio.

- **class `GeometricStiffnessOperator(element_length: float) -> None`**
  Assemble \(\mathbf{K}_\sigma\) from section forces and shape-function derivatives (**TL nonlinear**).
- **class `GreenLagrangeStrainOperator(element_length: float, include_shear: bool = False) -> None`**
  Green–Lagrange strain measures for a 2-node 3D beam in Total Lagrangian formulation.
- **class `StressResultantOperator() -> None`**
  Section forces (**2PK-type resultants**) from Green–Lagrange strain \(\mathbf{E}\) in the **reference** configuration.

### `onedfem.pre_processing.element_library.nonlinear.tosdt.reddy`

[src/onedfem/pre_processing/element_library/nonlinear/tosdt/reddy/__init__.py](../../src/onedfem/pre_processing/element_library/nonlinear/tosdt/reddy/__init__.py)

Nonlinear (Total Lagrangian) Reddy 3rd-order beam (16 DOF).

- **class `NonlinearReddyBeamElement3D(**kwargs)`**
  TL Reddy: variational ``K = int B8.T D B8`` plus GL ``K_delta`` and Hessian ``K_sigma``.

### `onedfem.pre_processing.element_library.nonlinear.tosdt.reddy.updated_lagrangian`

[src/onedfem/pre_processing/element_library/nonlinear/tosdt/reddy/updated_lagrangian/__init__.py](../../src/onedfem/pre_processing/element_library/nonlinear/tosdt/reddy/updated_lagrangian/__init__.py)

- **class `UpdatedLagrangianReddyBeamElement3D(**kwargs)`**
  UL Reddy: variational ``K_0`` on the current chord plus beam-column ``K_sigma``.

### `onedfem.pre_processing.element_library.nonlinear.tosdt.reddy.utilities`

[src/onedfem/pre_processing/element_library/nonlinear/tosdt/reddy/utilities/__init__.py](../../src/onedfem/pre_processing/element_library/nonlinear/tosdt/reddy/utilities/__init__.py)

TL operators for nonlinear Reddy — shared ``nonlinear.tosdt.utilities`` trio.

- **class `GeometricStiffnessOperator(element_length: float) -> None`**
  Assemble \(\mathbf{K}_\sigma\) from section forces and shape-function derivatives (**TL nonlinear**).
- **class `GreenLagrangeStrainOperator(element_length: float, include_shear: bool = False) -> None`**
  Green–Lagrange strain measures for a 2-node 3D beam in Total Lagrangian formulation.
- **class `StressResultantOperator() -> None`**
  Section forces (**2PK-type resultants**) from Green–Lagrange strain \(\mathbf{E}\) in the **reference** configuration.

### `onedfem.pre_processing.element_library.nonlinear.tosdt.utilities`

[src/onedfem/pre_processing/element_library/nonlinear/tosdt/utilities/__init__.py](../../src/onedfem/pre_processing/element_library/nonlinear/tosdt/utilities/__init__.py)

Shared TL operators for nonlinear tosdt (Reddy / Levinson).

- **class `GeometricStiffnessOperator(element_length: float) -> None`**
  Assemble \(\mathbf{K}_\sigma\) from section forces and shape-function derivatives (**TL nonlinear**).
- **class `GreenLagrangeStrainOperator(element_length: float, include_shear: bool = False) -> None`**
  Green–Lagrange strain measures for a 2-node 3D beam in Total Lagrangian formulation.
- **class `StressResultantOperator() -> None`**
  Section forces (**2PK-type resultants**) from Green–Lagrange strain \(\mathbf{E}\) in the **reference** configuration.

### `onedfem.pre_processing.element_library.rotating_frame`

[src/onedfem/pre_processing/element_library/rotating_frame.py](../../src/onedfem/pre_processing/element_library/rotating_frame.py)

Beams in a frame rotating at a constant rate: centrifugal load, spin softening and Coriolis terms.

- **class `Rotation(omega: float, axis: Tuple[float, float, float] = (0.0, 0.0, 1.0), origin: Tuple[float, float, float] = (0.0, 0.0, 0.0)) -> None`**
  Constant rotation ``omega`` (rad/s) about ``axis`` through ``origin`` (global axes).
- **class `RotatingFrameGaussPointData(xi: float, weight: float, jacobian: float, positions: npt.NDArray[np.float64], offsets: npt.NDArray[np.float64], masses: npt.NDArray[np.float64], Phi: npt.NDArray[np.float64], Theta: npt.NDArray[np.float64], forces: npt.NDArray[np.float64]) -> None`**
  One Gauss station of an element in the rotating frame.
- **class `RotatingFrameObject(element_id: int, rotation: Rotation, f_c: npt.NDArray[np.float64], K_soft: npt.NDArray[np.float64], K_rot: npt.NDArray[np.float64], G: npt.NDArray[np.float64], gauss_data: List[RotatingFrameGaussPointData] = <factory>) -> None`**
  Rotating-frame element terms (global DOFs), with their Gauss data.
- **`section_point_masses(mass: Dict[str, float]) -> Tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]`**
  Four point masses ``(yz (4, 2), m (4,))`` with the section's mass, mass centre and second moments.
- **`element_rotating_frame(element, p0: npt.ArrayLike, p1: npt.ArrayLike, rotation: Rotation, R0: Optional[np.ndarray] = None, R1: Optional[np.ndarray] = None, n_dof: Optional[int] = None, n_gauss: int = 4) -> RotatingFrameObject`**
  Rotating-frame terms of one straight element from ``p0`` to ``p1`` (global), node frames ``R0``, ``R1``.
- **`rotation_from_settings(simulation_settings: Dict[str, Any]) -> Optional[Rotation]`**
  The ``[Rotation]`` settings as a ``Rotation``, or None when absent or ``omega`` is zero.

### `onedfem.pre_processing.element_library.utilities`

[src/onedfem/pre_processing/element_library/utilities/__init__.py](../../src/onedfem/pre_processing/element_library/utilities/__init__.py)

- **`expand_stiffness_matrix(reduced_Ke, full_size=12, dof_map_binary=None)`**
  Expands a reduced stiffness matrix (Ke) to fit a full DOF system using efficient vectorized indexing.
- **`expand_force_vector(reduced_Fe, full_size=12, dof_map_binary=None)`**
  Expands a reduced force vector (Fe) to fit a full DOF system using vectorized mapping.
- **`get_gauss_points(n, dim=1)`**
  Get Gauss-Legendre quadrature points and weights for 1D, 2D, or 3D integration in natural coordinates.
- **`integrate_vector(n_gauss, integrand_func, jacobian_func, dim=1)`**
  Perform numerical integration of a vector over the element in natural coordinates using Gauss quadrature.
- **`integrate_matrix(n_gauss, integrand_func, jacobian_func, dim=1)`**
  Perform numerical integration of a matrix over the element in natural coordinates using Gauss quadrature.
- **`compute_jacobian_matrix(shape_function_derivatives, node_coordinates)`**
  Compute the Jacobian matrix for an element based on shape function derivatives and nodal coordinates.
- **`compute_jacobian_determinant(jacobian_matrix)`**
  Compute the determinant of the Jacobian matrix.
- **`general_jacobian_and_determinant(shape_function_derivatives, node_coordinates)`**
  Compute both the Jacobian matrix and its determinant in one step.
- **class `GreenLagrangeStrainOperator(element_length: float, include_shear: bool = False) -> None`**
  Minimal Green–Lagrange strain helper (legacy).
- **class `StressResultantOperator() -> None`**
  Section forces from ``S = D @ E`` at a point (same as EB/Timoshenko utilities).
- **class `GeometricStiffnessOperator(element_length: float) -> None`**
  Assemble \(\mathbf{K}_\sigma\) from section forces and shape-function derivatives (**TL nonlinear**).

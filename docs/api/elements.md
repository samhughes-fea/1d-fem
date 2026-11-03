# Elements

Beam, bar and truss elements and their shared utilities.

## Element types

The `element_type` names a job's `element.txt` accepts.

| Element type | Module | Summary |
|---|---|---|
| `CorotationalBeamElement3D` | `pre_processing.element_library.nonlinear.large_rotations.corotational.corotational_3D` | Co-rotational 3D beam (``kernel`` = ``timoshenko`` or ``euler_bernoulli``). |
| `CorotationalBenscoterBeamElement3D` | `pre_processing.element_library.nonlinear.large_rotations.corotational.corotational_section` | Corotational Benscoter beam: restrained warping and large twist on thin-walled sections. |
| `CorotationalEulerBernoulliBeamElement3D` | `pre_processing.element_library.nonlinear.large_rotations.corotational.corotational_section` | Corotational Euler-Bernoulli beam (shear-rigid core, classical section stiffness). |
| `CorotationalGBTBeamElement3D` | `pre_processing.element_library.nonlinear.large_rotations.corotational.corotational_section` | Corotational GBT beam: large rotations of the section, with its warping and distortional modes linear in the element's own frame. |
| `CorotationalLevinsonBeamElement3D` | `pre_processing.element_library.nonlinear.large_rotations.corotational.corotational_section` | Corotational Levinson beam (non-variational third-order shear; unsymmetric tangent). |
| `CorotationalReddyBeamElement3D` | `pre_processing.element_library.nonlinear.large_rotations.corotational.corotational_section` | Corotational Reddy beam: third-order shear (shear lag of the slaved warping) at large rotation. |
| `CorotationalTimoshenkoBeamElement3D` | `pre_processing.element_library.nonlinear.large_rotations.corotational.corotational_section` | Corotational Timoshenko beam on the coupled Saint-Venant section stiffness. |
| `GEBTUnshearableBeamElement3D` | `pre_processing.element_library.nonlinear.large_rotations.gebt_unshearable.gebt_unshearable_3D` | Unshearable large-rotation beam using the nonlinear EB TL pipeline. |
| `GeometricallyExactShearDeformableBeam3D` | `pre_processing.element_library.nonlinear.large_rotations.geometrically_exact_shear_deformable_beam.geometrically_exact_shear_deformable_beam_3D` | Factory-facing **GeometricallyExactShearDeformableBeam3D** element type. |
| `LinearBarElement3D` | `pre_processing.element_library.linear.bar.linear_bar_3D` | 2 nodes, 6 DOF/node; only axial and torsion strains are non-zero in ``B``. |
| `LinearBenscoterBeamElement3D` | `pre_processing.element_library.linear.beam.first_order_shear_deformation_theory.benscoter.linear_benscoter_3D` | Timoshenko beam with Benscoter restrained warping for thin-walled anisotropic sections. |
| `LinearEulerBernoulliBeamElement3D` | `pre_processing.element_library.linear.beam.zero_order_shear_deformation_theory.euler_bernoulli.linear_euler_bernoulli_3D` | 2-node straight beam, 6 DOF per node, 12 total (see module docstring for ``U_e`` order). |
| `LinearGBTBeamElement3D` | `pre_processing.element_library.linear.beam.first_order_shear_deformation_theory.benscoter.linear_gbt_3D` | Generalised Beam Theory element for thin-walled anisotropic sections. |
| `LinearLevinsonBeamElement3D` | `pre_processing.element_library.linear.beam.third_order_shear_deformation_theory.levinson.linear_levinson_3D` | Levinson third-order shear beam (isotropic or thin-walled sections). |
| `LinearReddyBeamElement3D` | `pre_processing.element_library.linear.beam.third_order_shear_deformation_theory.reddy.linear_reddy_3D` | Reddy third-order shear beam (isotropic or thin-walled sections). |
| `LinearTimoshenkoBeamElement3D` | `pre_processing.element_library.linear.beam.first_order_shear_deformation_theory.timoshenko.linear_timoshenko_3D` | 2-node straight beam, 6 DOF per node, 12 total; local ``x`` along the chord. |
| `LinearTrussElement3D` | `pre_processing.element_library.linear.truss.linear_truss_3D` | 2 nodes, 6 DOF/node; transverse displacement follows the element transverse direction (utilities). |
| `NonlinearEulerBernoulliBeamElement3D` | `pre_processing.element_library.nonlinear.euler_bernoulli.nonlinear_euler_bernoulli_3D` | TL EB element: Green-Lagrange strain, ``S = D @ E``, material and geometric tangents from Gauss sums. |
| `NonlinearLevinsonBeamElement3D` | `pre_processing.element_library.nonlinear.levinson.nonlinear_levinson_3D` | Total-Lagrangian 3D Levinson 3rd-order beam element (12 DOF) — STUB. |
| `NonlinearReddyBeamElement3D` | `pre_processing.element_library.nonlinear.reddy.nonlinear_reddy_3D` | Total-Lagrangian 3D Reddy 3rd-order beam element (12 DOF) — STUB. |
| `NonlinearTimoshenkoBeamElement3D` | `pre_processing.element_library.nonlinear.timoshenko.nonlinear_timoshenko_3D` | TL Timoshenko: ``K_0`` (linear selective assembly); ``E_nl`` and ``B_tot`` for ``F_int`` and ``K_delta``; ``K_sigma`` from ``N``, ``M_y``, ``M_z``. |
| `UpdatedLagrangianEulerBernoulliBeamElement3D` | `pre_processing.element_library.nonlinear.euler_bernoulli.updated_lagrangian_euler_bernoulli_3D` | Updated-Lagrangian 3D Euler-Bernoulli beam element (12 DOF) — STUB. |
| `UpdatedLagrangianTimoshenkoBeamElement3D` | `pre_processing.element_library.nonlinear.timoshenko.updated_lagrangian_timoshenko_3D` | Updated-Lagrangian 3D Timoshenko beam element (12 DOF) — STUB. |

## Modules

### `pre_processing.element_library`

[pre_processing/element_library/__init__.py](../../pre_processing/element_library/__init__.py)

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

### `pre_processing.element_library.beam_load_interpolation`

[pre_processing/element_library/beam_load_interpolation.py](../../pre_processing/element_library/beam_load_interpolation.py)

Work-conjugate load interpolation for two-node, six-DOF-per-node beams.

- **`beam_field_interpolation(xi: float, L: float, bending: str = 'hermite') -> np.ndarray`**
  ``N_f`` (6, 12): ``[u, v, w, theta_x, theta_y, theta_z](xi) = N_f U_e``, right-handed.
- **`point_load_vector(xi: float, L: float, load6, bending: str = 'hermite') -> np.ndarray`**
  Nodal loads (12,) of a point load ``load6`` at ``xi``.
- **`distributed_load_vector(xi_g, w_g, q_g, L: float, bending: str = 'hermite') -> np.ndarray`**
  Nodal loads (12,) of line loads ``q_g`` (n_g, 6) at Gauss points ``xi_g``, weights ``w_g``.

### `pre_processing.element_library.element_frame`

[pre_processing/element_library/element_frame.py](../../pre_processing/element_library/element_frame.py)

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

### `pre_processing.element_library.framed_element`

[pre_processing/element_library/framed_element.py](../../pre_processing/element_library/framed_element.py)

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

### `pre_processing.element_library.linear`

[pre_processing/element_library/linear/__init__.py](../../pre_processing/element_library/linear/__init__.py)

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

### `pre_processing.element_library.linear.bar`

[pre_processing/element_library/linear/bar/__init__.py](../../pre_processing/element_library/linear/bar/__init__.py)

- **class `LinearBarElement3D(*, element_id: int, element_dictionary: dict, grid_dictionary: dict, section_dictionary: dict, material_dictionary: dict, point_load_array: numpy.ndarray, distributed_load_array: numpy.ndarray, job_results_dir: str, quadrature_order: Optional[int] = None)`**
  2 nodes, 6 DOF/node; only axial and torsion strains are non-zero in ``B``.

### `pre_processing.element_library.linear.bar.linear_bar_3D`

[pre_processing/element_library/linear/bar/linear_bar_3D.py](../../pre_processing/element_library/linear/bar/linear_bar_3D.py)

2-node 3D bar: axial + torsion only (no bending/shear DOF coupling in B).

- **class `LinearBarElement3D(*, element_id: int, element_dictionary: dict, grid_dictionary: dict, section_dictionary: dict, material_dictionary: dict, point_load_array: numpy.ndarray, distributed_load_array: numpy.ndarray, job_results_dir: str, quadrature_order: Optional[int] = None)`**
  2 nodes, 6 DOF/node; only axial and torsion strains are non-zero in ``B``.

### `pre_processing.element_library.linear.bar.utilities`

[pre_processing/element_library/linear/bar/utilities/__init__.py](../../pre_processing/element_library/linear/bar/utilities/__init__.py)

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

### `pre_processing.element_library.linear.beam.first_order_shear_deformation_theory.benscoter`

[pre_processing/element_library/linear/beam/first_order_shear_deformation_theory/benscoter/__init__.py](../../pre_processing/element_library/linear/beam/first_order_shear_deformation_theory/benscoter/__init__.py)

- **class `LinearBenscoterBeamElement3D(**kwargs)`**
  Timoshenko beam with Benscoter restrained warping for thin-walled anisotropic sections.
- **class `LinearGBTBeamElement3D(**kwargs)`**
  Generalised Beam Theory element for thin-walled anisotropic sections.

### `pre_processing.element_library.linear.beam.first_order_shear_deformation_theory.timoshenko`

[pre_processing/element_library/linear/beam/first_order_shear_deformation_theory/timoshenko/__init__.py](../../pre_processing/element_library/linear/beam/first_order_shear_deformation_theory/timoshenko/__init__.py)

- **class `LinearTimoshenkoBeamElement3D(*, element_id: int, element_dictionary: dict, grid_dictionary: dict, section_dictionary: dict, material_dictionary: dict, point_load_array: numpy.ndarray, distributed_load_array: numpy.ndarray, job_results_dir: str, quadrature_order: int = 3)`**
  2-node straight beam, 6 DOF per node, 12 total; local ``x`` along the chord.

### `pre_processing.element_library.linear.beam.first_order_shear_deformation_theory.timoshenko.linear_timoshenko_3D`

[pre_processing/element_library/linear/beam/first_order_shear_deformation_theory/timoshenko/linear_timoshenko_3D.py](../../pre_processing/element_library/linear/beam/first_order_shear_deformation_theory/timoshenko/linear_timoshenko_3D.py)

2-node 3D Timoshenko beam (shear-deformable).

- **class `LinearTimoshenkoBeamElement3D(*, element_id: int, element_dictionary: dict, grid_dictionary: dict, section_dictionary: dict, material_dictionary: dict, point_load_array: numpy.ndarray, distributed_load_array: numpy.ndarray, job_results_dir: str, quadrature_order: int = 3)`**
  2-node straight beam, 6 DOF per node, 12 total; local ``x`` along the chord.

### `pre_processing.element_library.linear.beam.second_order`

[pre_processing/element_library/linear/beam/second_order.py](../../pre_processing/element_library/linear/beam/second_order.py)

Second-order (geometric) stiffness of the linear beam elements from their section forces.

- **`geometric_section_matrix(c: np.ndarray) -> np.ndarray`**
  ``S`` (4, 4) on ``[v', w', phi', phi]`` from the section terms ``c`` (7,).
- **class `GeometricStiffnessMixin()`**
  ``linear_geometric_stiffness_matrix`` for an element that defines ``_geometric_state``.

### `pre_processing.element_library.linear.beam.third_order_shear_deformation_theory.reddy`

[pre_processing/element_library/linear/beam/third_order_shear_deformation_theory/reddy/__init__.py](../../pre_processing/element_library/linear/beam/third_order_shear_deformation_theory/reddy/__init__.py)

- **class `LinearReddyBeamElement3D(**kwargs)`**
  Reddy third-order shear beam (isotropic or thin-walled sections).

### `pre_processing.element_library.linear.beam.zero_order_shear_deformation_theory.euler_bernoulli`

[pre_processing/element_library/linear/beam/zero_order_shear_deformation_theory/euler_bernoulli/__init__.py](../../pre_processing/element_library/linear/beam/zero_order_shear_deformation_theory/euler_bernoulli/__init__.py)

- **class `LinearEulerBernoulliBeamElement3D(*, element_id: int, element_dictionary: dict, grid_dictionary: dict, section_dictionary: dict, material_dictionary: dict, point_load_array: numpy.ndarray, distributed_load_array: numpy.ndarray, job_results_dir: str, quadrature_order: Optional[int] = None)`**
  2-node straight beam, 6 DOF per node, 12 total (see module docstring for ``U_e`` order).

### `pre_processing.element_library.linear.beam.zero_order_shear_deformation_theory.euler_bernoulli.linear_euler_bernoulli_3D`

[pre_processing/element_library/linear/beam/zero_order_shear_deformation_theory/euler_bernoulli/linear_euler_bernoulli_3D.py](../../pre_processing/element_library/linear/beam/zero_order_shear_deformation_theory/euler_bernoulli/linear_euler_bernoulli_3D.py)

2-node 3D Euler–Bernoulli beam.

- **class `LinearEulerBernoulliBeamElement3D(*, element_id: int, element_dictionary: dict, grid_dictionary: dict, section_dictionary: dict, material_dictionary: dict, point_load_array: numpy.ndarray, distributed_load_array: numpy.ndarray, job_results_dir: str, quadrature_order: Optional[int] = None)`**
  2-node straight beam, 6 DOF per node, 12 total (see module docstring for ``U_e`` order).

### `pre_processing.element_library.linear.beam.zero_order_shear_deformation_theory.euler_bernoulli_with_warp.utilities`

[pre_processing/element_library/linear/beam/zero_order_shear_deformation_theory/euler_bernoulli_with_warp/utilities/__init__.py](../../pre_processing/element_library/linear/beam/zero_order_shear_deformation_theory/euler_bernoulli_with_warp/utilities/__init__.py)

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

### `pre_processing.element_library.linear.truss`

[pre_processing/element_library/linear/truss/__init__.py](../../pre_processing/element_library/linear/truss/__init__.py)

- **class `LinearTrussElement3D(*, element_id: int, element_dictionary: dict, grid_dictionary: dict, section_dictionary: dict, material_dictionary: dict, point_load_array: numpy.ndarray, distributed_load_array: numpy.ndarray, job_results_dir: str, quadrature_order: Optional[int] = None)`**
  2 nodes, 6 DOF/node; transverse displacement follows the element transverse direction (utilities).

### `pre_processing.element_library.linear.truss.linear_truss_3D`

[pre_processing/element_library/linear/truss/linear_truss_3D.py](../../pre_processing/element_library/linear/truss/linear_truss_3D.py)

2-node 3D truss: axial + one transverse shear direction + torsion (six DOF/node for layout compatibility).

- **class `LinearTrussElement3D(*, element_id: int, element_dictionary: dict, grid_dictionary: dict, section_dictionary: dict, material_dictionary: dict, point_load_array: numpy.ndarray, distributed_load_array: numpy.ndarray, job_results_dir: str, quadrature_order: Optional[int] = None)`**
  2 nodes, 6 DOF/node; transverse displacement follows the element transverse direction (utilities).

### `pre_processing.element_library.linear.truss.utilities`

[pre_processing/element_library/linear/truss/utilities/__init__.py](../../pre_processing/element_library/linear/truss/utilities/__init__.py)

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

### `pre_processing.element_library.nonlinear`

[pre_processing/element_library/nonlinear/__init__.py](../../pre_processing/element_library/nonlinear/__init__.py)

- **class `NonlinearEulerBernoulliBeamElement3D(*, element_id: int, element_dictionary: dict, grid_dictionary: dict, section_dictionary: dict, material_dictionary: dict, point_load_array: numpy.ndarray, distributed_load_array: numpy.ndarray, job_results_dir: str, quadrature_order: Optional[int] = None)`**
  TL EB element: Green-Lagrange strain, ``S = D @ E``, material and geometric tangents from Gauss sums.
- **class `NonlinearTimoshenkoBeamElement3D(*, element_id: int, element_dictionary: dict, grid_dictionary: dict, section_dictionary: dict, material_dictionary: dict, point_load_array: numpy.ndarray, distributed_load_array: numpy.ndarray, job_results_dir: str, quadrature_order: int | None = 3)`**
  TL Timoshenko: ``K_0`` (linear selective assembly); ``E_nl`` and ``B_tot`` for ``F_int`` and ``K_delta``; ``K_sigma`` from ``N``, ``M_y``, ``M_z``.
- **class `NonlinearLevinsonBeamElement3D(*, element_id: int, element_dictionary: dict, grid_dictionary: dict, section_dictionary: dict, material_dictionary: dict, point_load_array: np.ndarray, distributed_load_array: np.ndarray, job_results_dir: str, quadrature_order: int | None = None)`**
  Total-Lagrangian 3D Levinson 3rd-order beam element (12 DOF) — STUB.
- **class `NonlinearReddyBeamElement3D(*, element_id: int, element_dictionary: dict, grid_dictionary: dict, section_dictionary: dict, material_dictionary: dict, point_load_array: np.ndarray, distributed_load_array: np.ndarray, job_results_dir: str, quadrature_order: int | None = None)`**
  Total-Lagrangian 3D Reddy 3rd-order beam element (12 DOF) — STUB.
- **class `UpdatedLagrangianEulerBernoulliBeamElement3D(*, element_id: int, element_dictionary: dict, grid_dictionary: dict, section_dictionary: dict, material_dictionary: dict, point_load_array: np.ndarray, distributed_load_array: np.ndarray, job_results_dir: str, quadrature_order: int | None = None)`**
  Updated-Lagrangian 3D Euler-Bernoulli beam element (12 DOF) — STUB.
- **class `UpdatedLagrangianTimoshenkoBeamElement3D(*, element_id: int, element_dictionary: dict, grid_dictionary: dict, section_dictionary: dict, material_dictionary: dict, point_load_array: np.ndarray, distributed_load_array: np.ndarray, job_results_dir: str, quadrature_order: int | None = None)`**
  Updated-Lagrangian 3D Timoshenko beam element (12 DOF) — STUB.
- **class `CorotationalBeamElement3D(*, element_id: int, element_dictionary: dict, grid_dictionary: dict, section_dictionary: dict, material_dictionary: dict, point_load_array: np.ndarray, distributed_load_array: np.ndarray, job_results_dir: str, quadrature_order: int | None = None, kernel: str = 'timoshenko', tangent_stiffness_mode: "Literal['finite_difference', 'elastic_material']" = 'finite_difference')`**
  Co-rotational 3D beam (``kernel`` = ``timoshenko`` or ``euler_bernoulli``).
- **class `GeometricallyExactShearDeformableBeam3D(*, element_id: int, element_dictionary: dict, grid_dictionary: dict, section_dictionary: dict, material_dictionary: dict, point_load_array: np.ndarray, distributed_load_array: np.ndarray, job_results_dir: str, quadrature_order: int | None = 3, gesdb_tl_fallback: bool = False, gesdb_kernel: str = 'tl_locked')`**
  Factory-facing **GeometricallyExactShearDeformableBeam3D** element type.
- **class `GEBTUnshearableBeamElement3D(*, element_id: int, element_dictionary: dict, grid_dictionary: dict, section_dictionary: dict, material_dictionary: dict, point_load_array: numpy.ndarray, distributed_load_array: numpy.ndarray, job_results_dir: str, quadrature_order: Optional[int] = None)`**
  Unshearable large-rotation beam using the nonlinear EB TL pipeline.

### `pre_processing.element_library.nonlinear.euler_bernoulli`

[pre_processing/element_library/nonlinear/euler_bernoulli/__init__.py](../../pre_processing/element_library/nonlinear/euler_bernoulli/__init__.py)

- **class `NonlinearEulerBernoulliBeamElement3D(*, element_id: int, element_dictionary: dict, grid_dictionary: dict, section_dictionary: dict, material_dictionary: dict, point_load_array: numpy.ndarray, distributed_load_array: numpy.ndarray, job_results_dir: str, quadrature_order: Optional[int] = None)`**
  TL EB element: Green-Lagrange strain, ``S = D @ E``, material and geometric tangents from Gauss sums.
- **class `UpdatedLagrangianEulerBernoulliBeamElement3D(*, element_id: int, element_dictionary: dict, grid_dictionary: dict, section_dictionary: dict, material_dictionary: dict, point_load_array: np.ndarray, distributed_load_array: np.ndarray, job_results_dir: str, quadrature_order: int | None = None)`**
  Updated-Lagrangian 3D Euler-Bernoulli beam element (12 DOF) — STUB.

### `pre_processing.element_library.nonlinear.euler_bernoulli.updated_lagrangian_euler_bernoulli_3D`

[pre_processing/element_library/nonlinear/euler_bernoulli/updated_lagrangian_euler_bernoulli_3D.py](../../pre_processing/element_library/nonlinear/euler_bernoulli/updated_lagrangian_euler_bernoulli_3D.py)

Updated-Lagrangian (UL) 3D Euler-Bernoulli beam element (12 DOF) — STUB.

- **class `UpdatedLagrangianEulerBernoulliBeamElement3D(*, element_id: int, element_dictionary: dict, grid_dictionary: dict, section_dictionary: dict, material_dictionary: dict, point_load_array: np.ndarray, distributed_load_array: np.ndarray, job_results_dir: str, quadrature_order: int | None = None)`**
  Updated-Lagrangian 3D Euler-Bernoulli beam element (12 DOF) — STUB.

### `pre_processing.element_library.nonlinear.euler_bernoulli.utilities`

[pre_processing/element_library/nonlinear/euler_bernoulli/utilities/__init__.py](../../pre_processing/element_library/nonlinear/euler_bernoulli/utilities/__init__.py)

- **class `GreenLagrangeStrainOperator(element_length: float, include_shear: bool = False) -> None`**
  Green–Lagrange strain measures for a 2-node 3D Euler–Bernoulli beam (Total Lagrangian).
- **class `StressResultantOperator() -> None`**
  Section forces (**2PK-type resultants**) from Green–Lagrange strain \(\mathbf{E}\) in the **reference** configuration.
- **class `GeometricStiffnessOperator(element_length: float) -> None`**
  Assemble \(\mathbf{K}_\sigma\) from section forces and shape-function derivatives (**TL nonlinear**).

### `pre_processing.element_library.nonlinear.large_rotations`

[pre_processing/element_library/nonlinear/large_rotations/__init__.py](../../pre_processing/element_library/nonlinear/large_rotations/__init__.py)

Finite-rotation beam kinematics (co-rotational frame, GEBT-class registrations).

- **class `CorotationalBeamElement3D(*, element_id: int, element_dictionary: dict, grid_dictionary: dict, section_dictionary: dict, material_dictionary: dict, point_load_array: np.ndarray, distributed_load_array: np.ndarray, job_results_dir: str, quadrature_order: int | None = None, kernel: str = 'timoshenko', tangent_stiffness_mode: "Literal['finite_difference', 'elastic_material']" = 'finite_difference')`**
  Co-rotational 3D beam (``kernel`` = ``timoshenko`` or ``euler_bernoulli``).
- **class `GeometricallyExactShearDeformableBeam3D(*, element_id: int, element_dictionary: dict, grid_dictionary: dict, section_dictionary: dict, material_dictionary: dict, point_load_array: np.ndarray, distributed_load_array: np.ndarray, job_results_dir: str, quadrature_order: int | None = 3, gesdb_tl_fallback: bool = False, gesdb_kernel: str = 'tl_locked')`**
  Factory-facing **GeometricallyExactShearDeformableBeam3D** element type.
- **class `GEBTUnshearableBeamElement3D(*, element_id: int, element_dictionary: dict, grid_dictionary: dict, section_dictionary: dict, material_dictionary: dict, point_load_array: numpy.ndarray, distributed_load_array: numpy.ndarray, job_results_dir: str, quadrature_order: Optional[int] = None)`**
  Unshearable large-rotation beam using the nonlinear EB TL pipeline.

### `pre_processing.element_library.nonlinear.large_rotations.corotational`

[pre_processing/element_library/nonlinear/large_rotations/corotational/__init__.py](../../pre_processing/element_library/nonlinear/large_rotations/corotational/__init__.py)

Co-rotational beam (12 DOF).

- **class `CorotationalBeamElement3D(*, element_id: int, element_dictionary: dict, grid_dictionary: dict, section_dictionary: dict, material_dictionary: dict, point_load_array: np.ndarray, distributed_load_array: np.ndarray, job_results_dir: str, quadrature_order: int | None = None, kernel: str = 'timoshenko', tangent_stiffness_mode: "Literal['finite_difference', 'elastic_material']" = 'finite_difference')`**
  Co-rotational 3D beam (``kernel`` = ``timoshenko`` or ``euler_bernoulli``).

### `pre_processing.element_library.nonlinear.large_rotations.corotational.corotational_section`

[pre_processing/element_library/nonlinear/large_rotations/corotational/corotational_section.py](../../pre_processing/element_library/nonlinear/large_rotations/corotational/corotational_section.py)

Corotational beam elements on any SectionAnalysis, with third- and fourth-order local energy.

- **`rotation_exp(theta: np.ndarray) -> np.ndarray`**
  Rotation matrix of the rotation vector ``theta`` (Rodrigues).
- **`rotation_log(R: np.ndarray) -> np.ndarray`**
  Rotation vector of ``R`` (through the quaternion, robust up to pi).
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
- **class `CorotationalReddyBeamElement3D(**kwargs)`**
  Corotational Reddy beam: third-order shear (shear lag of the slaved warping) at large rotation.
- **class `CorotationalLevinsonBeamElement3D(**kwargs)`**
  Corotational Levinson beam (non-variational third-order shear; unsymmetric tangent).

### `pre_processing.element_library.nonlinear.large_rotations.gebt_unshearable`

[pre_processing/element_library/nonlinear/large_rotations/gebt_unshearable/__init__.py](../../pre_processing/element_library/nonlinear/large_rotations/gebt_unshearable/__init__.py)

GEBT unshearable (Euler-Bernoulli) beam (12 DOF).

- **class `GEBTUnshearableBeamElement3D(*, element_id: int, element_dictionary: dict, grid_dictionary: dict, section_dictionary: dict, material_dictionary: dict, point_load_array: numpy.ndarray, distributed_load_array: numpy.ndarray, job_results_dir: str, quadrature_order: Optional[int] = None)`**
  Unshearable large-rotation beam using the nonlinear EB TL pipeline.

### `pre_processing.element_library.nonlinear.large_rotations.geometrically_exact_shear_deformable_beam`

[pre_processing/element_library/nonlinear/large_rotations/geometrically_exact_shear_deformable_beam/__init__.py](../../pre_processing/element_library/nonlinear/large_rotations/geometrically_exact_shear_deformable_beam/__init__.py)

- **class `GeometricallyExactShearDeformableBeam3D(*, element_id: int, element_dictionary: dict, grid_dictionary: dict, section_dictionary: dict, material_dictionary: dict, point_load_array: np.ndarray, distributed_load_array: np.ndarray, job_results_dir: str, quadrature_order: int | None = 3, gesdb_tl_fallback: bool = False, gesdb_kernel: str = 'tl_locked')`**
  Factory-facing **GeometricallyExactShearDeformableBeam3D** element type.

### `pre_processing.element_library.nonlinear.levinson`

[pre_processing/element_library/nonlinear/levinson/__init__.py](../../pre_processing/element_library/nonlinear/levinson/__init__.py)

Nonlinear (Total Lagrangian) Levinson 3rd-order beam (12 DOF).

- **class `NonlinearLevinsonBeamElement3D(*, element_id: int, element_dictionary: dict, grid_dictionary: dict, section_dictionary: dict, material_dictionary: dict, point_load_array: np.ndarray, distributed_load_array: np.ndarray, job_results_dir: str, quadrature_order: int | None = None)`**
  Total-Lagrangian 3D Levinson 3rd-order beam element (12 DOF) — STUB.

### `pre_processing.element_library.nonlinear.levinson.nonlinear_levinson_3D`

[pre_processing/element_library/nonlinear/levinson/nonlinear_levinson_3D.py](../../pre_processing/element_library/nonlinear/levinson/nonlinear_levinson_3D.py)

Nonlinear (Total Lagrangian) 3D Levinson third-order shear deformation beam (12 DOF) — STUB.

- **class `NonlinearLevinsonBeamElement3D(*, element_id: int, element_dictionary: dict, grid_dictionary: dict, section_dictionary: dict, material_dictionary: dict, point_load_array: np.ndarray, distributed_load_array: np.ndarray, job_results_dir: str, quadrature_order: int | None = None)`**
  Total-Lagrangian 3D Levinson 3rd-order beam element (12 DOF) — STUB.

### `pre_processing.element_library.nonlinear.reddy`

[pre_processing/element_library/nonlinear/reddy/__init__.py](../../pre_processing/element_library/nonlinear/reddy/__init__.py)

Nonlinear (Total Lagrangian) Reddy 3rd-order beam (12 DOF).

- **class `NonlinearReddyBeamElement3D(*, element_id: int, element_dictionary: dict, grid_dictionary: dict, section_dictionary: dict, material_dictionary: dict, point_load_array: np.ndarray, distributed_load_array: np.ndarray, job_results_dir: str, quadrature_order: int | None = None)`**
  Total-Lagrangian 3D Reddy 3rd-order beam element (12 DOF) — STUB.

### `pre_processing.element_library.nonlinear.reddy.nonlinear_reddy_3D`

[pre_processing/element_library/nonlinear/reddy/nonlinear_reddy_3D.py](../../pre_processing/element_library/nonlinear/reddy/nonlinear_reddy_3D.py)

Nonlinear (Total Lagrangian) 3D Reddy third-order shear deformation beam (12 DOF) — STUB.

- **class `NonlinearReddyBeamElement3D(*, element_id: int, element_dictionary: dict, grid_dictionary: dict, section_dictionary: dict, material_dictionary: dict, point_load_array: np.ndarray, distributed_load_array: np.ndarray, job_results_dir: str, quadrature_order: int | None = None)`**
  Total-Lagrangian 3D Reddy 3rd-order beam element (12 DOF) — STUB.

### `pre_processing.element_library.nonlinear.timoshenko`

[pre_processing/element_library/nonlinear/timoshenko/__init__.py](../../pre_processing/element_library/nonlinear/timoshenko/__init__.py)

- **class `NonlinearTimoshenkoBeamElement3D(*, element_id: int, element_dictionary: dict, grid_dictionary: dict, section_dictionary: dict, material_dictionary: dict, point_load_array: numpy.ndarray, distributed_load_array: numpy.ndarray, job_results_dir: str, quadrature_order: int | None = 3)`**
  TL Timoshenko: ``K_0`` (linear selective assembly); ``E_nl`` and ``B_tot`` for ``F_int`` and ``K_delta``; ``K_sigma`` from ``N``, ``M_y``, ``M_z``.
- **class `UpdatedLagrangianTimoshenkoBeamElement3D(*, element_id: int, element_dictionary: dict, grid_dictionary: dict, section_dictionary: dict, material_dictionary: dict, point_load_array: np.ndarray, distributed_load_array: np.ndarray, job_results_dir: str, quadrature_order: int | None = None)`**
  Updated-Lagrangian 3D Timoshenko beam element (12 DOF) — STUB.

### `pre_processing.element_library.nonlinear.timoshenko.updated_lagrangian_timoshenko_3D`

[pre_processing/element_library/nonlinear/timoshenko/updated_lagrangian_timoshenko_3D.py](../../pre_processing/element_library/nonlinear/timoshenko/updated_lagrangian_timoshenko_3D.py)

Updated-Lagrangian (UL) 3D Timoshenko beam element (12 DOF) — STUB.

- **class `UpdatedLagrangianTimoshenkoBeamElement3D(*, element_id: int, element_dictionary: dict, grid_dictionary: dict, section_dictionary: dict, material_dictionary: dict, point_load_array: np.ndarray, distributed_load_array: np.ndarray, job_results_dir: str, quadrature_order: int | None = None)`**
  Updated-Lagrangian 3D Timoshenko beam element (12 DOF) — STUB.

### `pre_processing.element_library.nonlinear.timoshenko.utilities`

[pre_processing/element_library/nonlinear/timoshenko/utilities/__init__.py](../../pre_processing/element_library/nonlinear/timoshenko/utilities/__init__.py)

- **class `GreenLagrangeStrainOperator(element_length: float, include_shear: bool = False) -> None`**
  Green–Lagrange strain measures for a 2-node 3D beam in Total Lagrangian formulation.
- **class `StressResultantOperator() -> None`**
  Section forces (**2PK-type resultants**) from Green–Lagrange strain \(\mathbf{E}\) in the **reference** configuration.
- **class `GeometricStiffnessOperator(element_length: float) -> None`**
  Assemble \(\mathbf{K}_\sigma\) from section forces and shape-function derivatives (**TL nonlinear**).

### `pre_processing.element_library.nonlinear.timoshenko.utilities.geometric_stiffness`

[pre_processing/element_library/nonlinear/timoshenko/utilities/geometric_stiffness.py](../../pre_processing/element_library/nonlinear/timoshenko/utilities/geometric_stiffness.py)

\(\mathbf{K}_\sigma\) for nonlinear Timoshenko — same weak-form implementation as EB (**re-export**).

- **class `GeometricStiffnessOperator(element_length: float) -> None`**
  Assemble \(\mathbf{K}_\sigma\) from section forces and shape-function derivatives (**TL nonlinear**).

### `pre_processing.element_library.rotating_frame`

[pre_processing/element_library/rotating_frame.py](../../pre_processing/element_library/rotating_frame.py)

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

### `pre_processing.element_library.utilities`

[pre_processing/element_library/utilities/__init__.py](../../pre_processing/element_library/utilities/__init__.py)

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

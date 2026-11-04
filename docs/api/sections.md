# Section analysis

Thin-walled and 2D (solid) cross-section analysis: stiffness, mass, warping, GBT modes, stress recovery, ply failure and layup sensitivities.

### `pre_processing.section_integration`

[pre_processing/section_integration/__init__.py](../../pre_processing/section_integration/__init__.py)

Section integration: geometry to A, I_y, I_z, J_t, kappa, shear centre (y_sc, z_sc), optional Gamma.

- **`rectangle_properties(width: float, height: float) -> Dict[str, float]`**
  Section properties for a rectangle (centroid at origin, width along y, height along z).
- **`i_section_properties(width_flange: float, height_web: float, thickness_flange: float, thickness_web: float) -> Dict[str, float]`**
  Section properties for a doubly symmetric I-section (centroid at origin).
- **`channel_properties(width_flange: float, height_web: float, thickness_flange: float, thickness_web: float) -> Dict[str, float]`**
  Section properties for a channel (open thin-walled): U-section with flanges at top and bottom.
- **`section_properties_from_geometry(section_type: str, **kwargs: float) -> Dict[str, float]`**
  Compute section properties from geometry (closed-form).

### `pre_processing.section_integration.gbt_alignment`

[pre_processing/section_integration/gbt_alignment.py](../../pre_processing/section_integration/gbt_alignment.py)

Consistent GBT bases on a beam whose sections change along it.

- **`MIN_MATCH`** (constant)
- **`align_gbt_sections(per_element: Sequence, connectivity: npt.ArrayLike, gbt_elements: Optional[Sequence[bool]] = None) -> Tuple[tuple, Dict[int, dict]]`**
  Turn the GBT bases of each element's section to its neighbour's (module docstring).
- **`gbt_signature(sa) -> Tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]`**
  ``(S_w, S_m)``: one row per warping mode and per distortional mode of ``sa``.
- **`mass_transform(sa, Rw: np.ndarray, Rm: np.ndarray) -> np.ndarray`**
  ``T`` on the GBT mass order: a generalised load ``f`` of the old basis is ``T^T f`` in the new.
- **`rotate_gbt_basis(sa, Rw: np.ndarray, Rm: np.ndarray)`**
  ``sa`` with its warping basis turned by ``Rw`` and its mode basis by ``Rm``.

### `pre_processing.section_integration.ply_failure`

[pre_processing/section_integration/ply_failure.py](../../pre_processing/section_integration/ply_failure.py)

Ply failure criteria: maximum stress, Tsai-Wu, Hashin and Puck, as exposures.

- **class `PlyStrength(Xt: float, Xc: float, Yt: float, Yc: float, S12: float, S23: Optional[float] = None, Zt: Optional[float] = None, Zc: Optional[float] = None, F_star: float = -0.5, p_TL_t: float = 0.3, p_TL_c: float = 0.25, p_TT: float = 0.25) -> None`**
  Ply strengths (all positive, compressive ones as magnitudes) and criterion parameters.
- **`CRITERIA`** (constant)
- **`max_stress(s, st: PlyStrength) -> Tuple[np.ndarray, np.ndarray]`**
  Maximum-stress exposure and its governing component (``"s11t"``, ``"s22c"``, ``"t12"`` ...).
- **`tsai_wu(s, st: PlyStrength) -> np.ndarray`**
  Tsai-Wu exposure (3D, ``F_ij = F_star sqrt(F_ii F_jj)`` for the normal-stress pairs).
- **`hashin(s, st: PlyStrength) -> Tuple[np.ndarray, np.ndarray]`**
  Hashin (1980) exposure and mode: ``"FT"``, ``"FC"``, ``"MT"`` or ``"MC"``.
- **`puck(s, st: PlyStrength) -> Tuple[np.ndarray, np.ndarray, np.ndarray]`**
  Puck exposure, mode (``"FFT"``, ``"FFC"``, ``"IFF-A"``, ``"IFF-B"``, ``"IFF-C"``) and the fracture-plane angle in degrees (for inter-fibre failure).
- **`exposures(s, st: PlyStrength) -> Dict[str, np.ndarray]`**
  Every criterion's exposure and mode: ``{"max_stress": f, "max_stress_mode": ..., ..., "puck_angle": degrees, "exposure": max over criteria, "margin": 1 / exposure - 1}``.

### `pre_processing.section_integration.section_analysis`

[pre_processing/section_integration/section_analysis.py](../../pre_processing/section_integration/section_analysis.py)

Section constitutive data shared by every beam element: ``SectionAnalysis``.

- **`to_element_voigt(K_rh: npt.ArrayLike) -> npt.NDArray[np.float64]`**
  Map a right-handed 6 x 6 section stiffness to the element sign convention.
- **class `SectionAnalysis(stiffness: npt.NDArray[np.float64], classical_stiffness: npt.NDArray[np.float64], warping_stiffness: float = 0.0, shear_centre: Tuple[float, float] = (0.0, 0.0), tension_centre: Tuple[float, float] = (0.0, 0.0), mass: Dict[str, float] = <factory>, source: str = 'isotropic', section_id: Optional[int] = None, benscoter_stiffness: Optional[npt.NDArray[np.float64]] = None, warping_inertia: float = 0.0, gbt_modes: int = 0, gbt_stiffness: Optional[npt.NDArray[np.float64]] = None, gbt_mass: Optional[npt.NDArray[np.float64]] = None, gbt_shear_modes: int = 0, gbt_sv_modes: int = 0, third_order_stiffness: Optional[npt.NDArray[np.float64]] = None, third_order_rate: Optional[npt.NDArray[np.float64]] = None, geometric_terms: Optional[npt.NDArray[np.float64]] = None, benscoter_geometric_terms: Optional[npt.NDArray[np.float64]] = None, third_order_geometric_terms: Optional[npt.NDArray[np.float64]] = None, axial_points: Optional[npt.NDArray[np.float64]] = None, gbt_geometric: Optional[Tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]] = None, gbt_moments: Optional[npt.NDArray[np.float64]] = None) -> None`**
  Section stiffness and properties in the element convention.

### `pre_processing.section_integration.solid`

[pre_processing/section_integration/solid/__init__.py](../../pre_processing/section_integration/solid/__init__.py)

2D finite-element analysis of beam cross-sections (thick, solid and laminated sections).

- **class `SolidMaterial(E1: float, E2: float, E3: float, G12: float, G13: float, G23: float, nu12: float, nu13: float, nu23: float, rho: float = 0.0) -> None`**
  Orthotropic material in its own axes ``1, 2, 3`` (engineering constants).
- **`bond_matrix(R: npt.ArrayLike) -> npt.NDArray[np.float64]`**
  Stress transformation ``sigma_section = T sigma_material`` for the rotation ``R`` (columns: material axes in section axes), Voigt order ``VOIGT``.
- **`material_axes(theta_deg: float, alpha_deg: float) -> npt.NDArray[np.float64]`**
  Columns ``e1, e2, e3`` in section axes ``(x, y, z)`` (see the module notes).
- **class `SectionMesh(nodes: npt.NDArray[np.float64], elements: List[npt.NDArray[np.int64]], kinds: List[str], material_ids: List[int], theta: List[float], alpha: List[float], materials: Dict[int, SolidMaterial] = <factory>, section_id: Optional[int] = None, ties: Dict[int, List[Tuple[int, float]]] = <factory>) -> None`**
  SectionMesh(nodes: 'npt.NDArray[np.float64]', elements: 'List[npt.NDArray[np.int64]]', kinds: 'List[str]', material_ids: 'List[int]', theta: 'List[float]', alpha: 'List[float]', materials: 'Dict[int, SolidMaterial]' = <factory>, section_id: 'Optional[int]' = None, ties: 'Dict[int, List[Tuple[int, float]]]' = <factory>)
- **`parse_section_mesh_2d(path: str, materials: Mapping[int, SolidMaterial]) -> Dict[int, SectionMesh]`**
  Read ``section_mesh_2d.txt``: ``{section_id: SectionMesh}``.
- **`read_gmsh(path: str)`**
  Nodes ``(n, 2)`` and ``[(physical tag, physical name, kind, connectivity)]`` of a Gmsh 4.1 ASCII mesh's ``T6`` and ``Q8`` elements on physical surfaces.
- **class `SolidSectionAnalysis(mesh: SectionMesh) -> None`**
  Saint-Venant stiffness and properties of a section meshed in 2D.

### `pre_processing.section_integration.solid.cross_section_modes`

[pre_processing/section_integration/solid/cross_section_modes.py](../../pre_processing/section_integration/solid/cross_section_modes.py)

Cross-section deformation modes of sections meshed in 2D: GBT on a continuum section.

- **class `SolidCrossSectionModes(analysis, max_modes: int = 12) -> None`**
  GBT modes, stiffness and mass of a ``SolidSectionAnalysis`` (right-handed).

### `pre_processing.section_integration.solid.element_recovery`

[pre_processing/section_integration/solid/element_recovery.py](../../pre_processing/section_integration/solid/element_recovery.py)

3D stresses over a 2D-meshed section at a point of a beam element.

- **`element_section_state(element, U_e: npt.ArrayLike, xi: float, analysis) -> Tuple[Optional[np.ndarray], Optional[np.ndarray]]`**
  ``(forces, benscoter)`` at ``xi``: right-handed section forces (6,) for Euler-Bernoulli and Timoshenko elements, or right-handed Benscoter strains (8,) for Benscoter elements (the other entry ``None``).
- **`element_section_stresses(element, U_e: npt.ArrayLike, xi: float, analysis) -> List`**
  ``PointStress`` at every quadrature point of the section at ``xi`` of an element.

### `pre_processing.section_integration.solid.elements`

[pre_processing/section_integration/solid/elements.py](../../pre_processing/section_integration/solid/elements.py)

Quadratic 2D section elements: 6-node triangles (``T6``) and 8-node quadrilaterals (``Q8``).

- **`n_nodes(kind: str) -> int`**
- **`shape(kind: str, xi: float, eta: float) -> Tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]`**
  Shape functions ``N`` (n,) and their natural derivatives ``dN`` (2, n) at ``(xi, eta)``.
- **`quadrature(kind: str) -> Tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]`**
  Points ``(m, 2)`` and weights ``(m,)``: 3 x 3 Gauss on ``Q8``, 7-point (degree 5) on ``T6``.

### `pre_processing.section_integration.solid.from_geometry`

[pre_processing/section_integration/solid/from_geometry.py](../../pre_processing/section_integration/solid/from_geometry.py)

Conforming 2D section meshes of an implicit-geometry-kernel blade section.

- **`section_mesh_from_cut(cut, laminates: Mapping[str, object], library: Optional[Mapping[int, object]] = None, element_size: Optional[float] = None, rows_per_ply: int = 1, thickness_tol: float = 0.01, origin: str = 'pitch_axis', materials: Optional[Mapping[object, SolidMaterial]] = None, section_id: Optional[int] = None) -> SectionMesh`**
  Layered ``Q8`` mesh of an implicit-geometry-kernel ``SectionCut``.
- **`write_section_mesh_2d(path: str, meshes: Mapping[int, SectionMesh]) -> None`**
  Write ``section_mesh_2d.txt`` for ``{section_id: SectionMesh}``, with its materials.

### `pre_processing.section_integration.solid.material`

[pre_processing/section_integration/solid/material.py](../../pre_processing/section_integration/solid/material.py)

Three-dimensional orthotropic materials for 2D section meshes, and their rotation to section axes.

- **class `SolidMaterial(E1: float, E2: float, E3: float, G12: float, G13: float, G23: float, nu12: float, nu13: float, nu23: float, rho: float = 0.0) -> None`**
  Orthotropic material in its own axes ``1, 2, 3`` (engineering constants).
- **`material_axes(theta_deg: float, alpha_deg: float) -> npt.NDArray[np.float64]`**
  Columns ``e1, e2, e3`` in section axes ``(x, y, z)`` (see the module notes).
- **`bond_matrix(R: npt.ArrayLike) -> npt.NDArray[np.float64]`**
  Stress transformation ``sigma_section = T sigma_material`` for the rotation ``R`` (columns: material axes in section axes), Voigt order ``VOIGT``.

### `pre_processing.section_integration.solid.mesh`

[pre_processing/section_integration/solid/mesh.py](../../pre_processing/section_integration/solid/mesh.py)

2D finite-element meshes of beam cross-sections, with a material and orientation per element.

- **class `SectionMesh(nodes: npt.NDArray[np.float64], elements: List[npt.NDArray[np.int64]], kinds: List[str], material_ids: List[int], theta: List[float], alpha: List[float], materials: Dict[int, SolidMaterial] = <factory>, section_id: Optional[int] = None, ties: Dict[int, List[Tuple[int, float]]] = <factory>) -> None`**
  SectionMesh(nodes: 'npt.NDArray[np.float64]', elements: 'List[npt.NDArray[np.int64]]', kinds: 'List[str]', material_ids: 'List[int]', theta: 'List[float]', alpha: 'List[float]', materials: 'Dict[int, SolidMaterial]' = <factory>, section_id: 'Optional[int]' = None, ties: 'Dict[int, List[Tuple[int, float]]]' = <factory>)
- **class `QuadraturePoint(element: int, y: float, z: float, weight: float, N: np.ndarray, dN: np.ndarray, nodes: np.ndarray) -> None`**
  QuadraturePoint(element: 'int', y: 'float', z: 'float', weight: 'float', N: 'np.ndarray', dN: 'np.ndarray', nodes: 'np.ndarray')
- **`parse_section_mesh_2d(path: str, materials: Mapping[int, SolidMaterial]) -> Dict[int, SectionMesh]`**
  Read ``section_mesh_2d.txt``: ``{section_id: SectionMesh}``.
- **`read_gmsh(path: str)`**
  Nodes ``(n, 2)`` and ``[(physical tag, physical name, kind, connectivity)]`` of a Gmsh 4.1 ASCII mesh's ``T6`` and ``Q8`` elements on physical surfaces.

### `pre_processing.section_integration.solid.section_analysis`

[pre_processing/section_integration/solid/section_analysis.py](../../pre_processing/section_integration/solid/section_analysis.py)

Saint-Venant analysis of beam sections meshed in 2D: ``SolidSectionAnalysis``.

- **class `SolidSectionAnalysis(mesh: SectionMesh) -> None`**
  Saint-Venant stiffness and properties of a section meshed in 2D.
- **class `PointStress(element: int, y: float, z: float, strain: npt.NDArray[np.float64], stress: npt.NDArray[np.float64], stress_material: npt.NDArray[np.float64], stress_wall: npt.NDArray[np.float64]) -> None`**
  3D state at a point of the section.
- **`BENSCOTER_STRAINS`** (constant)

### `pre_processing.section_integration.thin_walled`

[pre_processing/section_integration/thin_walled/__init__.py](../../pre_processing/section_integration/thin_walled/__init__.py)

Thin-walled anisotropic sections: wall laminates, condensed wall laws, mid-line geometry, Saint-Venant analysis.

- **`DEFAULT_SHEAR_CORRECTION`** (constant)
- **class `Ply(material: PlyMaterial, angle_deg: float, thickness: float) -> None`**
  One ply of a wall laminate.
- **class `PlyMaterial(E1: float, E2: float, G12: float, nu12: float, G13: float, G23: float, rho: float = 0.0) -> None`**
  Orthotropic ply material in its own axes (1 = fibre, 2 = transverse, 3 = normal).
- **class `WallLaminate(plies: Sequence[Ply], shear_correction: float = 0.8333333333333334)`**
  A wall laminate: plies stacked from bottom to top about the mid-surface.
- **`FORCES`** (constant)
- **`STRAINS`** (constant)
- **class `ThinWalledSectionAnalysis(section: ThinWalledSection, laminates: Mapping[int, WallLaminate], hoop: str = 'free', n_sub: int = 1) -> None`**
  Saint-Venant stiffness and section properties of a thin-walled section.
- **`stiffness_about(K: npt.ArrayLike, y0: float, z0: float) -> npt.NDArray[np.float64]`**
  Stiffness (6, 6) about the point ``(y0, z0)`` from a stiffness about the origin.
- **class `ThinWalledSection(points: Mapping[int, Tuple[float, float]], walls: Iterable[Wall], laminate_ids: Optional[Iterable[int]] = None, section_id: int = 0) -> None`**
  Validated mid-line geometry of a thin-walled section.
- **class `Wall(wall_id: int, start: int, end: int, laminate_id: int) -> None`**
  A straight wall between two section points.
- **`WALL_RESULTANTS`** (constant)
- **`WALL_STRAINS`** (constant)
- **class `WallStiffness(matrix: npt.NDArray[np.float64], strains: Tuple[str, ...]) -> None`**
  A condensed wall stiffness and the names of its kept strains.
- **`beam_wall_stiffness(laminate: WallLaminate) -> WallStiffness`**
  Wall law for thin-walled beam theory: ``N_ss = M_ss = 0``, (4, 4).
- **`condense(K: npt.NDArray[np.float64], released: Iterable[int]) -> Tuple[npt.NDArray[np.float64], Tuple[int, ...]]`**
  Statically condense out the ``released`` rows/columns of a symmetric stiffness.
- **`condensed_wall_stiffness(laminate: WallLaminate, released: Iterable[str] = ('N_ss', 'M_ss')) -> WallStiffness`**
  Condensed in-plane wall stiffness with the named resultants released.
- **`gbt_wall_stiffness(laminate: WallLaminate) -> WallStiffness`**
  Wall law for cross-section distortion modes: ``N_ss = 0``, (5, 5).

### `pre_processing.section_integration.thin_walled.cross_section_modes`

[pre_processing/section_integration/thin_walled/cross_section_modes.py](../../pre_processing/section_integration/thin_walled/cross_section_modes.py)

Cross-section deformation modes of thin-walled sections (conventional GBT modes).

- **class `CrossSectionModes(analysis: ThinWalledSectionAnalysis) -> None`**
  Conventional (Vlasov) cross-section modes and the GBT section matrices.

### `pre_processing.section_integration.thin_walled.element_recovery`

[pre_processing/section_integration/thin_walled/element_recovery.py](../../pre_processing/section_integration/thin_walled/element_recovery.py)

Wall and ply stresses at a point of a beam element with a thin-walled section.

- **`element_wall_field(element, U_e: npt.ArrayLike, xi: float, analysis, wall_id: int) -> Tuple[KeptField, str]`**
  Kept-strain field ``s -> (strains, x rate)`` along wall ``wall_id`` at ``xi``.
- **`element_wall_stress(element, U_e: npt.ArrayLike, xi: float, analysis, wall_id: int, s: float, points_per_ply: int = 3) -> WallPointStress`**
  Recovered wall and ply state at ``xi`` of an element and ``s`` along a wall.

### `pre_processing.section_integration.thin_walled.from_geometry`

[pre_processing/section_integration/thin_walled/from_geometry.py](../../pre_processing/section_integration/thin_walled/from_geometry.py)

Thin-walled sections cut from an implicit-geometry-kernel blade.

- **class `GeometrySection(section: ThinWalledSection, laminates: Dict[int, WallLaminate], wall_names: Dict[int, str], wall_layers: Dict[int, Tuple[str, ...]], stack_names: Dict[int, str], origin: Tuple[float, float], z: float) -> None`**
  A thin-walled section built from a ``SectionCut``.
- **`thin_walled_from_cut(cut, laminates: LaminateMap, library: Optional[Mapping[int, WallLaminate]] = None, thickness_tol: float = 0.01, origin: str = 'pitch_axis', section_id: int = 0) -> GeometrySection`**
  ``GeometrySection`` of an implicit-geometry-kernel ``SectionCut``.
- **`laminate_lookup(wall, laminates: LaminateMap, library: Optional[Mapping[int, WallLaminate]] = None) -> WallLaminate`**
  The laminate of a ``SectionWall``: by label, then region, then kind.

### `pre_processing.section_integration.thin_walled.laminate`

[pre_processing/section_integration/thin_walled/laminate.py](../../pre_processing/section_integration/thin_walled/laminate.py)

Classical laminate theory for the walls of a thin-walled beam section.

- **class `PlyMaterial(E1: float, E2: float, G12: float, nu12: float, G13: float, G23: float, rho: float = 0.0) -> None`**
  Orthotropic ply material in its own axes (1 = fibre, 2 = transverse, 3 = normal).
- **`strain_rotation(theta_rad: float) -> npt.NDArray[np.float64]`**
  Engineering in-plane strain transform ``T_eps`` with ``eps_ply = T_eps @ eps_wall``.
- **`shear_rotation(theta_rad: float) -> npt.NDArray[np.float64]`**
  Transverse shear strain transform with ``[gamma_13, gamma_23] = R @ [gamma_xn, gamma_sn]``.
- **class `Ply(material: PlyMaterial, angle_deg: float, thickness: float) -> None`**
  One ply of a wall laminate.
- **class `WallLaminate(plies: Sequence[Ply], shear_correction: float = 0.8333333333333334)`**
  A wall laminate: plies stacked from bottom to top about the mid-surface.

### `pre_processing.section_integration.thin_walled.layup_sensitivity`

[pre_processing/section_integration/thin_walled/layup_sensitivity.py](../../pre_processing/section_integration/thin_walled/layup_sensitivity.py)

Exact derivatives of thin-walled section properties with respect to ply angles and thicknesses.

- **class `LayupParameter(laminate_id: int, ply: int, kind: str) -> None`**
  Ply ``ply`` (0-based, bottom first) of laminate ``laminate_id``: its ``angle`` (per degree) or ``thickness`` (per metre).
- **`laminate_derivative(lam: WallLaminate, param: LayupParameter) -> Tuple[np.ndarray, np.ndarray]`**
  ``d ABD`` (6, 6) and ``d [m0, m1, m2]`` (the mass moments ``int rho n^p dn``, p = 0, 1, 2).
- **`wall_law_derivative(lam: WallLaminate, dABD: np.ndarray, hoop: str = 'free') -> np.ndarray`**
  ``d wall_law(lam, hoop)`` (4, 4) for a change ``dABD`` of the laminate stiffness.
- **class `SectionDerivative(parameter: LayupParameter, flexibility: npt.NDArray[np.float64], stiffness: npt.NDArray[np.float64], mass: Dict[str, float], mass_matrix: npt.NDArray[np.float64], benscoter_stiffness: Optional[npt.NDArray[np.float64]] = None, warping_inertia: float = 0.0, geometric_terms: Optional[npt.NDArray[np.float64]] = None) -> None`**
  Derivatives of a thin-walled section for one ``LayupParameter``.
- **`section_derivative(an, param: LayupParameter, benscoter: bool = False, geometric: bool = False) -> SectionDerivative`**
  Exact derivatives of the section stiffness, flexibility and mass for ``param``.
- **`element_section_derivative(sa, d: SectionDerivative)`**
  A ``SectionAnalysis`` like ``sa`` whose stiffness and mass matrix are the derivatives ``d`` (element signs).

### `pre_processing.section_integration.thin_walled.section_analysis`

[pre_processing/section_integration/thin_walled/section_analysis.py](../../pre_processing/section_integration/thin_walled/section_analysis.py)

Saint-Venant analysis of anisotropic thin-walled sections.

- **`wall_law(laminate: WallLaminate, hoop: str = 'free') -> npt.NDArray[np.float64]`**
  Wall stiffness (4, 4) in order ``WALL_STRAINS``.
- **class `ThinWalledSectionAnalysis(section: ThinWalledSection, laminates: Mapping[int, WallLaminate], hoop: str = 'free', n_sub: int = 1) -> None`**
  Saint-Venant stiffness and section properties of a thin-walled section.
- **`translation_matrix(dy: float, dz: float) -> npt.NDArray[np.float64]`**
  Strain map ``psi_new = P @ psi`` when the reference point moves by ``(dy, dz)``.
- **`stiffness_about(K: npt.ArrayLike, y0: float, z0: float) -> npt.NDArray[np.float64]`**
  Stiffness (6, 6) about the point ``(y0, z0)`` from a stiffness about the origin.

### `pre_processing.section_integration.thin_walled.stress_recovery`

[pre_processing/section_integration/thin_walled/stress_recovery.py](../../pre_processing/section_integration/thin_walled/stress_recovery.py)

Wall and ply stresses of a thin-walled section from its beam solution.

- **`mid_surface_strains(laminate: WallLaminate, kept: npt.ArrayLike, hoop: str = 'free') -> Array`**
  All six mid-surface strains from the strains a beam theory carries.
- **`wall_resultants(laminate: WallLaminate, strains6: npt.ArrayLike) -> Array`**
  ``[N_xx, N_ss, N_xs, M_xx, M_ss, M_xs]`` for mid-surface strains (6,).
- **class `PlyPoint(ply: int, n: float, strain_wall: Array, stress_wall: Array, strain_ply: Array, stress_ply: Array, tau_xn: float = 0.0, tau_sn: float = 0.0, sigma_nn: float = 0.0, tau_13: float = 0.0, tau_23: float = 0.0) -> None`**
  Stresses at one height ``n`` through a wall (``n = 0`` the mid-surface).
- **`through_thickness(laminate: WallLaminate, strains6: npt.ArrayLike, rate_x: Optional[npt.ArrayLike] = None, rate_s: Optional[npt.ArrayLike] = None, rate_xs: Optional[npt.ArrayLike] = None, rate_ss: Optional[npt.ArrayLike] = None, points_per_ply: int = 3) -> List[PlyPoint]`**
  Ply stresses at ``points_per_ply`` heights per ply (bottom, ..., top of each ply).
- **class `WallPointStress(wall_id: int, s: float, position: Array, strains: Array, resultants: Array, transverse_shear: Array, plies: List[PlyPoint]) -> None`**
  Recovered state at one point of a wall.
- **`wall_point_stress(laminate: WallLaminate, field: KeptField, wall_id: int, s: float, position: npt.ArrayLike, element_span: Tuple[float, float], hoop: str = 'free', points_per_ply: int = 3) -> WallPointStress`**
  Recover the state at ``s`` from a kept-strain field along the wall.
- **class `SaintVenantRecovery(analysis) -> None`**
  Wall and ply stresses of a ``ThinWalledSectionAnalysis`` under section forces.

### `pre_processing.section_integration.thin_walled.taper`

[pre_processing/section_integration/thin_walled/taper.py](../../pre_processing/section_integration/thin_walled/taper.py)

Taper coupling of a thin-walled section: the Resal effect.

- **class `SectionTaper(coupling: npt.NDArray[np.float64], shear: npt.NDArray[np.float64]) -> None`**
  Right-handed taper coupling of one section.
- **`section_taper(start: ThinWalledSection, end: ThinWalledSection, laminates: Mapping[int, WallLaminate], distance: float, n_gauss: int = 3) -> SectionTaper`**
  Taper coupling of the section midway between ``start`` and ``end``.
- **`section_taper_from_geometry(geometry, z: float, laminates, step: float, library: Optional[Mapping[int, WallLaminate]] = None, thickness_tol: float = 0.01) -> SectionTaper`**
  Taper coupling at span ``z`` of an implicit-geometry-kernel blade, from cuts at ``z -+ step``.

### `pre_processing.section_integration.thin_walled.thin_walled_geometry`

[pre_processing/section_integration/thin_walled/thin_walled_geometry.py](../../pre_processing/section_integration/thin_walled/thin_walled_geometry.py)

Mid-line geometry of a thin-walled beam section.

- **class `Wall(wall_id: int, start: int, end: int, laminate_id: int) -> None`**
  A straight wall between two section points.
- **class `ThinWalledSection(points: Mapping[int, Tuple[float, float]], walls: Iterable[Wall], laminate_ids: Optional[Iterable[int]] = None, section_id: int = 0) -> None`**
  Validated mid-line geometry of a thin-walled section.

### `pre_processing.section_integration.thin_walled.wall_stiffness`

[pre_processing/section_integration/thin_walled/wall_stiffness.py](../../pre_processing/section_integration/thin_walled/wall_stiffness.py)

Condensed wall stiffness for beam and cross-section-mode theories.

- **class `WallStiffness(matrix: npt.NDArray[np.float64], strains: Tuple[str, ...]) -> None`**
  A condensed wall stiffness and the names of its kept strains.
- **`condense(K: npt.NDArray[np.float64], released: Iterable[int]) -> Tuple[npt.NDArray[np.float64], Tuple[int, ...]]`**
  Statically condense out the ``released`` rows/columns of a symmetric stiffness.
- **`condensed_wall_stiffness(laminate: WallLaminate, released: Iterable[str] = ('N_ss', 'M_ss')) -> WallStiffness`**
  Condensed in-plane wall stiffness with the named resultants released.
- **`beam_wall_stiffness(laminate: WallLaminate) -> WallStiffness`**
  Wall law for thin-walled beam theory: ``N_ss = M_ss = 0``, (4, 4).
- **`gbt_wall_stiffness(laminate: WallLaminate) -> WallStiffness`**
  Wall law for cross-section distortion modes: ``N_ss = 0``, (5, 5).

# Section analysis

The section provider seam and the section input and recovery that 1d-fem keeps. The cross-section analysis itself (thin-walled and 2D sections, GBT modes, stress recovery, ply failure, layup sensitivities) comes from the section provider, anisotropic-section-model by default, and so do the isotropic `section.txt` sections.

### `onedfem.sections`

[src/onedfem/sections.py](../../src/onedfem/sections.py)

The section provider: where 1d-fem gets its cross-section analyses.

- **`DEFAULT_SECTION_PROVIDER`** (constant)
- **`PROVIDER_MODULES`** (constant)
- **`SECTION_API`** (constant)
- **class `SectionProviderError`**
  No section provider is installed, or it lacks a module 1d-fem needs.
- **`check_api_version(module: ModuleType) -> Tuple[int, int]`**
  The contract version ``module`` declares, checked against ``SECTION_API``.
- **`check_section_provider(provider: Union[str, ModuleType, None] = None) -> List[str]`**
  Modules and names of ``PROVIDER_MODULES`` that ``provider`` (default: the current one) lacks.
- **`section_import(name: str, *names: str)`**
  ``from <provider>.<name> import <names>``: the objects, in order (one object for one name).
- **`section_module(name: str) -> ModuleType`**
  The provider's module ``name`` (for example ``"thin_walled.laminate"``).
- **`section_provider() -> ModuleType`**
  The section provider package, imported on first use and checked against ``SECTION_API``.
- **`set_section_provider(provider: Union[str, ModuleType, None]) -> None`**
  Use ``provider`` (a package name or module) for sections; ``None`` restores the default.

### `onedfem.pre_processing.element_library.isotropic_section`

[src/onedfem/pre_processing/element_library/isotropic_section.py](../../src/onedfem/pre_processing/element_library/isotropic_section.py)

Isotropic section properties, from the section provider.

- **`isotropic_section(*, E: float, G: float, A: float, I_y: float = 0.0, I_z: float = 0.0, J_t: float = 0.0, kappa: float = 0.8333333333333334, I_x: Optional[float] = None, y_sc: float = 0.0, z_sc: float = 0.0, Gamma: float = 0.0, rho: float = 0.0)`**
  The provider's ``SectionAnalysis.from_isotropic`` for these constants (cached, read-only).
- **`isotropic_section_of(element)`**
  ``isotropic_section`` from an element's ``section_array`` and ``material_array``.
- **`isotropic_stiffness(**constants) -> np.ndarray`**
  Section ``D`` (6, 6) on ``[eps, kappa_y, kappa_z, gamma_xy, gamma_xz, phi]`` (a copy).
- **`isotropic_mass(section) -> Dict[str, float]`**
  Mass per length of a section: ``m``, ``I_p`` (polar), ``I_yy``, ``I_zz`` (``rho`` times areas).
- **`nodal_mass_diagonal(section, n_dof: int = 12) -> np.ndarray`**
  Per-DOF mass densities ``mu`` of a two-node element from the section mass.

### `onedfem.post_processing.thin_walled.element_recovery`

[src/onedfem/post_processing/thin_walled/element_recovery.py](../../src/onedfem/post_processing/thin_walled/element_recovery.py)

Wall and ply stresses at a point of a beam element with a thin-walled section.

- **`element_wall_field(element, U_e: npt.ArrayLike, xi: float, analysis, wall_id: int) -> Tuple[KeptField, str]`**
  Kept-strain field ``s -> (strains, x rate)`` along wall ``wall_id`` at ``xi``.
- **`element_wall_stress(element, U_e: npt.ArrayLike, xi: float, analysis, wall_id: int, s: float, points_per_ply: int = 3) -> WallPointStress`**
  Recovered wall and ply state at ``xi`` of an element and ``s`` along a wall.

### `onedfem.post_processing.thin_walled.wall_stresses`

[src/onedfem/post_processing/thin_walled/wall_stresses.py](../../src/onedfem/post_processing/thin_walled/wall_stresses.py)

Write wall and ply stresses of a solved job with thin-walled sections.

- **`element_wall_rows(element, U_e, analysis, points_per_wall: int = 5, xi: float = 0.0, mode='', wall_points: Optional[Sequence] = None) -> List[list]`**
  CSV rows (``WALL_STRESS_COLUMNS``) for one element; ``x`` is the global ``X`` of the axis at ``xi``.
- **`write_wall_stresses(job_results_dir: str, elements: Sequence, analyses: Sequence, U_global: np.ndarray, points_per_wall: Optional[int] = None, stations: Optional[Sequence] = None, modes: Optional[np.ndarray] = None, wall_points: Optional[Sequence] = None) -> str`**
  Write ``wall_stresses.csv`` for the elements that have a thin-walled analysis; return its path.

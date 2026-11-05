# Section analysis

The section provider seam and the section input and recovery that 1d-fem keeps. The cross-section analysis itself (thin-walled and 2D sections, GBT modes, stress recovery, ply failure, layup sensitivities) comes from the section provider, anisotropic-section-model by default.

### `onedfem.sections`

[src/onedfem/sections.py](../../src/onedfem/sections.py)

The section provider: where 1d-fem gets its cross-section analyses.

- **`DEFAULT_SECTION_PROVIDER`** (constant)
- **`PROVIDER_MODULES`** (constant)
- **class `SectionProviderError`**
  No section provider is installed, or it lacks a module 1d-fem needs.
- **`check_section_provider(provider: Union[str, ModuleType, None] = None) -> List[str]`**
  Modules and names of ``PROVIDER_MODULES`` that ``provider`` (default: the current one) lacks.
- **`section_import(name: str, *names: str)`**
  ``from <provider>.<name> import <names>``: the objects, in order (one object for one name).
- **`section_module(name: str) -> ModuleType`**
  The provider's module ``name`` (for example ``"thin_walled.laminate"``).
- **`section_provider() -> ModuleType`**
  The section provider package, imported on first use.
- **`set_section_provider(provider: Union[str, ModuleType, None]) -> None`**
  Use ``provider`` (a package name or module) for sections; ``None`` restores the default.

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

- **`element_wall_rows(element, U_e, analysis, points_per_wall: int = 5, xi: float = 0.0, mode='') -> List[list]`**
  CSV rows (``WALL_STRESS_COLUMNS``) for one element; ``x`` is the global ``X`` of the axis at ``xi``.
- **`write_wall_stresses(job_results_dir: str, elements: Sequence, analyses: Sequence, U_global: np.ndarray, points_per_wall: Optional[int] = None, stations: Optional[Sequence] = None, modes: Optional[np.ndarray] = None) -> str`**
  Write ``wall_stresses.csv`` for the elements that have a thin-walled analysis; return its path.

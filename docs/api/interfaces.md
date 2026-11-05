# Interfaces

Adapters from other repositories' data (a blade's geometry, duck typed) to sections and jobs. The core does not depend on them.

### `onedfem.interfaces`

[src/onedfem/interfaces/__init__.py](../../src/onedfem/interfaces/__init__.py)

Adapters from other repositories' data to 1d-fem jobs.


### `onedfem.interfaces.blade`

[src/onedfem/interfaces/blade/__init__.py](../../src/onedfem/interfaces/blade/__init__.py)

Blade adapters: a blade's geometry (duck typed, e.g. from implicit-geometry-kernel) to thin-walled and 2D sections and to a whole beam job.


### `onedfem.interfaces.blade.blade_beam_job`

[src/onedfem/interfaces/blade/blade_beam_job.py](../../src/onedfem/interfaces/blade/blade_beam_job.py)

A 1d-fem beam job for a blade or wing lofted by implicit-geometry-kernel.

- **class `BladeBeamJob(job_dir: str, kind: str, z: np.ndarray, nodes: np.ndarray, z_cut: np.ndarray, twist: np.ndarray, pretwist: np.ndarray, sections: Dict[int, object], element_type: str, loads: Dict[int, np.ndarray] = <factory>, tapers: Dict[int, object] = <factory>) -> None`**
  What ``blade_beam_job`` wrote.
- **`blade_beam_job(job_dir: str, geometry, laminates: Mapping[str, object], n_elements: Optional[int] = None, kind: str = 'thin_walled', *, z_nodes: Optional[Sequence[float]] = None, library: Optional[Mapping[int, object]] = None, element_type: str = 'LinearTimoshenkoBeamElement3D', pretwist: bool = True, taper: bool = True, n_sub: int = 1, mesh_options: Optional[Mapping[str, object]] = None, loads: Optional[Mapping[float, Sequence[float]]] = None, line_loads: Optional[Mapping[float, Sequence[float]]] = None, thickness_tol: float = 0.01, settings: Optional[str] = None, gbt_modes: int = 0) -> BladeBeamJob`**
  Write a cantilever beam job for ``geometry`` into ``job_dir`` (replaced).
- **`blade_to_job(v: Sequence[float]) -> np.ndarray`**
  Blade-frame ``(x, y, z)`` components (or a force-moment 6-vector) in job axes ``(X, Y, Z)``.
- **`job_to_blade(v: Sequence[float]) -> np.ndarray`**
  Job-axis ``(X, Y, Z)`` components (or a 6-vector) in the blade frame ``(x, y, z)``.
- **`NONLINEAR_SETTINGS`** (constant)

### `onedfem.interfaces.blade.solid_sections`

[src/onedfem/interfaces/blade/solid_sections.py](../../src/onedfem/interfaces/blade/solid_sections.py)

Conforming 2D section meshes of an implicit-geometry-kernel blade section.

- **`section_mesh_from_cut(cut, laminates: Mapping[str, object], library: Optional[Mapping[int, object]] = None, element_size: Optional[float] = None, rows_per_ply: int = 1, thickness_tol: float = 0.01, origin: str = 'pitch_axis', materials: Optional[Mapping[object, SolidMaterial]] = None, section_id: Optional[int] = None) -> SectionMesh`**
  Layered ``Q8`` mesh of an implicit-geometry-kernel ``SectionCut``.
- **`write_section_mesh_2d(path: str, meshes: Mapping[int, SectionMesh]) -> None`**
  Write ``section_mesh_2d.txt`` for ``{section_id: SectionMesh}``, with its materials.

### `onedfem.interfaces.blade.thin_walled_sections`

[src/onedfem/interfaces/blade/thin_walled_sections.py](../../src/onedfem/interfaces/blade/thin_walled_sections.py)

Thin-walled sections cut from an implicit-geometry-kernel blade.

- **class `GeometrySection(section: ThinWalledSection, laminates: Dict[int, WallLaminate], wall_names: Dict[int, str], wall_layers: Dict[int, Tuple[str, ...]], stack_names: Dict[int, str], origin: Tuple[float, float], z: float) -> None`**
  A thin-walled section built from a ``SectionCut``.
- **`thin_walled_from_cut(cut, laminates: LaminateMap, library: Optional[Mapping[int, WallLaminate]] = None, thickness_tol: float = 0.01, origin: str = 'pitch_axis', section_id: int = 0) -> GeometrySection`**
  ``GeometrySection`` of an implicit-geometry-kernel ``SectionCut``.
- **`laminate_lookup(wall, laminates: LaminateMap, library: Optional[Mapping[int, WallLaminate]] = None) -> WallLaminate`**
  The laminate of a ``SectionWall``: by label, then region, then kind.
- **`section_taper_from_geometry(geometry, z: float, laminates, step: float, library: Optional[Mapping[int, WallLaminate]] = None, thickness_tol: float = 0.01) -> SectionTaper`**
  Taper coupling at span ``z`` of an implicit-geometry-kernel blade, from cuts at ``z -+ step``.

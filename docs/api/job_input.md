# Job input

Parsers for the job directory's text files, meshes and loads.

### `onedfem.pre_processing.parsing`

[src/onedfem/pre_processing/parsing/__init__.py](../../src/onedfem/pre_processing/parsing/__init__.py)


### `onedfem.pre_processing.parsing.distributed_load_parser`

[src/onedfem/pre_processing/parsing/distributed_load_parser.py](../../src/onedfem/pre_processing/parsing/distributed_load_parser.py)

- **`parse_distributed_load(file_path)`**
  Parses discrete distributed load vectors from a structured text file and returns a 2D NumPy array.

### `onedfem.pre_processing.parsing.element_parser`

[src/onedfem/pre_processing/parsing/element_parser.py](../../src/onedfem/pre_processing/parsing/element_parser.py)

- **class `ElementParser(filepath: str, job_results_dir: str) -> None`**
  Parses an [Element] section and returns a dictionary whose structure is identical to the other parsers in the suite:

### `onedfem.pre_processing.parsing.element_section_parser`

[src/onedfem/pre_processing/parsing/element_section_parser.py](../../src/onedfem/pre_processing/parsing/element_section_parser.py)

Parse ``element_section.txt``: which thin-walled section each element uses.

- **`parse_element_section_file(filepath: str) -> Dict[int, int]`**
  Return ``{element_id: section_id}``.
- **`parse_element_section_modes(filepath: str) -> Dict[int, Tuple[int, int, int]]`**
  Return ``{element_id: (section_id, modes, shear_modes)}`` (0 for absent columns).

### `onedfem.pre_processing.parsing.grid_parser`

[src/onedfem/pre_processing/parsing/grid_parser.py](../../src/onedfem/pre_processing/parsing/grid_parser.py)

- **class `GridParser(filepath: str, job_results_dir: str) -> None`**
  Parses a [Grid] section and returns::

### `onedfem.pre_processing.parsing.laminate_parser`

[src/onedfem/pre_processing/parsing/laminate_parser.py](../../src/onedfem/pre_processing/parsing/laminate_parser.py)

Parse ``laminate.txt``: ply materials and wall laminates for thin-walled sections.

- **`preprocess_lines(filepath: str) -> List[str]`**
  Stripped non-empty lines with ``#`` comments (whole-line or trailing) removed.
- **`split_blocks(lines: Sequence[str], headers: Sequence[str]) -> Dict[str, Tuple[str, List[str]]]`**
  Split lines into ``{header: (subheader_line, data_lines)}`` for the known block headers.
- **`tokens_match(line: str, expected: Sequence[str]) -> bool`**
- **`to_float(tok: str, what: str, line: str) -> float`**
- **`to_int(tok: str, what: str, line: str) -> int`**
- **`parse_ply_materials(sub: str, rows: Sequence[str]) -> Dict[int, PlyMaterial]`**
- **`parse_laminates(sub: str, rows: Sequence[str], materials: Dict[int, PlyMaterial]) -> Dict[int, WallLaminate]`**
- **`parse_laminate_file(filepath: str) -> Dict[int, WallLaminate]`**
  Read ``laminate.txt`` and return ``{laminate_id: WallLaminate}``.
- **`parse_ply_material_file(filepath: str) -> Dict[int, PlyMaterial]`**
  Read only the ``[Ply_Material]`` block of ``laminate.txt``: ``{material_id: PlyMaterial}``.
- **`format_laminate_file(laminates: Dict[int, WallLaminate]) -> str`**
  Text of a ``laminate.txt`` for ``{laminate_id: WallLaminate}``.

### `onedfem.pre_processing.parsing.layup_sensitivity_parser`

[src/onedfem/pre_processing/parsing/layup_sensitivity_parser.py](../../src/onedfem/pre_processing/parsing/layup_sensitivity_parser.py)

Parse ``layup_sensitivity.txt``: the ply parameters and responses of a sensitivity run.

- **class `Response(type: str, node_id: Optional[int] = None, dof: Optional[str] = None, mode: Optional[int] = None) -> None`**
  A response of a sensitivity run (see the module docstring).
- **`parse_layup_sensitivity_file(path: str) -> Tuple[List[LayupParameter], List[Response]]`**
  Parameters and responses from ``path``; two empty lists when the file is absent.
- **`format_layup_sensitivity_file(params: Sequence[LayupParameter], responses: Sequence[Response]) -> str`**
  The text of ``layup_sensitivity.txt`` for ``params`` and ``responses``.

### `onedfem.pre_processing.parsing.material_parser`

[src/onedfem/pre_processing/parsing/material_parser.py](../../src/onedfem/pre_processing/parsing/material_parser.py)

- **class `MaterialParser(filepath: str, job_results_dir: str) -> None`**
  Parses a [Material] section and returns

### `onedfem.pre_processing.parsing.ply_strength_parser`

[src/onedfem/pre_processing/parsing/ply_strength_parser.py](../../src/onedfem/pre_processing/parsing/ply_strength_parser.py)

Parse ``ply_strength.txt``: ply strengths per material id, for the failure criteria.

- **`parse_ply_strength_file(path: str) -> Dict[int, PlyStrength]`**
  ``{material_id: PlyStrength}`` from ``path``; an empty dict when the file is absent.
- **`format_ply_strength_file(strengths: Mapping[int, PlyStrength]) -> str`**
  Text of ``ply_strength.txt`` with every column.

### `onedfem.pre_processing.parsing.point_load_parser`

[src/onedfem/pre_processing/parsing/point_load_parser.py](../../src/onedfem/pre_processing/parsing/point_load_parser.py)

- **`parse_point_load(file_path)`**
  Parses point load vector/s from a structured text file and returns a 2D NumPy array.

### `onedfem.pre_processing.parsing.precurvature_parser`

[src/onedfem/pre_processing/parsing/precurvature_parser.py](../../src/onedfem/pre_processing/parsing/precurvature_parser.py)

Parse optional per-element reference curvature / twist for straight beam elements.

- **`parse_precurvature(file_path: Optional[str], element_ids: npt.NDArray[np.int64]) -> npt.NDArray[np.float64]`**
  Read ``precurvature.txt`` and return ``(N_e, 3)`` columns ``[k_x0, k_y0, k_z0]`` (1/m), row order matching *element_ids* (same order as ``element_dictionary['ids']``).
- **`element_reference_strain_voigt(element_dictionary: Dict[str, Any], element_id: int) -> npt.NDArray[np.float64]`**
  Return ``(6,)`` Voigt ``E_0`` for *element_id* from ``precurvature_per_element`` if present, else zeros.
- **`reference_strain_voigt(k_xyz: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]`**
  Build 6-vector ``E_0`` matching Voigt ``[eps_x, kappa_y, kappa_z, gamma_xy, gamma_xz, phi_x]``.
- **`voigt_standard_beam_to_third_order_beam(E_standard: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]`**
  Map standard beam Voigt ``[ε_x, κ_y, κ_z, …]`` to Levinson/Reddy row order ``[ε_x, κ_z, κ_y, …]``.

### `onedfem.pre_processing.parsing.prescribed_displacement_parser`

[src/onedfem/pre_processing/parsing/prescribed_displacement_parser.py](../../src/onedfem/pre_processing/parsing/prescribed_displacement_parser.py)

- **`parse_prescribed_displacement(file_path, dof_per_node: int = 6)`**
  Parses prescribed displacement conditions from a structured text file.

### `onedfem.pre_processing.parsing.region_laminate_parser`

[src/onedfem/pre_processing/parsing/region_laminate_parser.py](../../src/onedfem/pre_processing/parsing/region_laminate_parser.py)

Parse ``region_laminate.txt``: which ``laminate.txt`` laminate each geometry region uses.

- **`parse_region_laminate_file(filepath: str, laminate_ids: Optional[Iterable[int]] = None) -> Dict[str, int]`**
  Read ``region_laminate.txt`` and return ``{label: laminate_id}``.
- **`format_region_laminate(mapping: Dict[str, int]) -> str`**
  Text of a ``region_laminate.txt`` for ``{label: laminate_id}``.

### `onedfem.pre_processing.parsing.section_geometry_parser`

[src/onedfem/pre_processing/parsing/section_geometry_parser.py](../../src/onedfem/pre_processing/parsing/section_geometry_parser.py)

Parse ``section_geometry.txt``: mid-line points and walls of thin-walled sections.

- **`parse_section_mesh(filepath: str) -> Dict[int, int]`**
  Return ``{section_id: n_sub}`` from the optional ``[Section_Mesh]`` block (empty without it).
- **`parse_section_geometry_file(filepath: str, laminate_ids: Optional[Iterable[int]] = None) -> Dict[int, ThinWalledSection]`**
  Read ``section_geometry.txt`` and return ``{section_id: ThinWalledSection}``.
- **`format_section_geometry(sections: Dict[int, ThinWalledSection], n_sub: Optional[Dict[int, int]] = None) -> str`**
  Text of a ``section_geometry.txt`` for ``{section_id: ThinWalledSection}``.

### `onedfem.pre_processing.parsing.section_parser`

[src/onedfem/pre_processing/parsing/section_parser.py](../../src/onedfem/pre_processing/parsing/section_parser.py)

- **class `SectionParser(filepath: str, job_results_dir: str) -> None`**
  Parses a [Section] block and returns

### `onedfem.pre_processing.parsing.section_taper_parser`

[src/onedfem/pre_processing/parsing/section_taper_parser.py](../../src/onedfem/pre_processing/parsing/section_taper_parser.py)

Parse ``section_taper.txt``: the taper coupling of each element's section (the Resal effect).

- **`parse_section_taper_file(filepath: str) -> "Dict[int, 'SectionTaper']"`**
  Return ``{element_id: SectionTaper}``.
- **`format_section_taper_file(tapers: "Dict[int, 'SectionTaper']") -> str`**
  Text of ``section_taper.txt`` for ``{element_id: SectionTaper}``.

### `onedfem.pre_processing.parsing.simulation_settings_parser`

[src/onedfem/pre_processing/parsing/simulation_settings_parser.py](../../src/onedfem/pre_processing/parsing/simulation_settings_parser.py)

- **`parse_simulation_settings(file_path)`**
  Parses simulation settings from a structured simulation_settings.txt file.

### `onedfem.pre_processing.parsing.simulation_settings_resolution`

[src/onedfem/pre_processing/parsing/simulation_settings_resolution.py](../../src/onedfem/pre_processing/parsing/simulation_settings_resolution.py)

Canonical taxonomy §1–§5 resolution for simulation_settings after parsing simulation_settings.txt.

- **`effective_transient_config(settings: Dict[str, Any]) -> Dict[str, Any]`**
  Merge ``[Transient]`` over legacy ``[Dynamic]``; transient keys win.
- **`effective_eigen_config(settings: Dict[str, Any]) -> Dict[str, Any]`**
  Vibration eigenproblem settings (legacy ``modal`` + ``eigen`` section).
- **`effective_buckling_config(settings: Dict[str, Any]) -> Dict[str, Any]`**
  Linear buckling settings (legacy ``modal`` + ``buckling`` section).
- **`finalize_simulation_settings(settings: Dict[str, Any], *, type_line_explicit: bool) -> None`**
  Mutates ``settings`` in place: hydrate taxonomy dicts, validate ``enabled`` sections, and set canonical ``settings['type']``.

### `onedfem.pre_processing.parsing.wall_load_parser`

[src/onedfem/pre_processing/parsing/wall_load_parser.py](../../src/onedfem/pre_processing/parsing/wall_load_parser.py)

Parse ``wall_load.txt``: uniform surface loads on the walls of thin-walled sections.

- **`parse_wall_load_file(filepath: str) -> Dict[int, Dict[int, Tuple[float, float, float]]]`**
  Return ``{element_id: {wall_id: (p_x, p_s, p_n)}}``.

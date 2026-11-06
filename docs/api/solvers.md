# Solvers

Static, nonlinear, eigen, buckling, harmonic, transient and spectral solvers and their runners.

### `onedfem.processing.buckling`

[src/onedfem/processing/buckling/__init__.py](../../src/onedfem/processing/buckling/__init__.py)

Buckling and stability kernels (taxonomy §5): geometric stiffness and linear buckling eigenpairs.

- **`assemble_global_geometric_stiffness(elements: List[Any], U_global: np.ndarray, total_dof: int) -> csr_matrix`**
  Assemble :math:`\mathbf{K}_g` (global geometric / stress stiffness) from element :meth:`linear_geometric_stiffness_matrix` and a reference displacement field.
- **`apply_buckling_boundary_conditions(K_global: csr_matrix, Kg_global: csr_matrix, prescribed_displacements: Optional[dict] = None, fixed_dofs: Optional[Sequence[int]] = None) -> Tuple[csr_matrix, csr_matrix, np.ndarray]`**
  Apply the **same** constrained DOFs as modal vibration: penalty on ``K``, zero rows/columns on ``Kg`` at constrained DOFs (no artificial geometric stiffness at supports).
- **`solve_linear_buckling_eigenpairs(K_mod: csr_matrix, Kg_mod: csr_matrix, num_modes: int, constrained_dofs: Optional[np.ndarray] = None) -> Tuple[np.ndarray, np.ndarray]`**
  Solve the generalized eigenproblem :math:`\mathbf{K}\mathbf{x} = \mu (-\mathbf{K}_g)\mathbf{x}` associated with :math:`(\mathbf{K} + \lambda \mathbf{K}_g)\mathbf{x}=\mathbf{0}`.

### `onedfem.processing.buckling.p_delta`

[src/onedfem/processing/buckling/p_delta.py](../../src/onedfem/processing/buckling/p_delta.py)

Second-order (P-delta) static solve: ``(K + K_sigma(U)) U = F`` by fixed-point iteration.

- **class `PDeltaResult(U: np.ndarray, converged: bool, iterations: int, history: List[float] = <factory>) -> None`**
  Displacements and the iteration record of a P-delta solve.
- **`solve_p_delta(K, F: np.ndarray, elements: Sequence[Any], fixed_dofs: Sequence[int], *, tolerance: float = 1e-08, max_iterations: int = 50, U0: Optional[np.ndarray] = None, solver: Optional[Callable[[Any, np.ndarray], np.ndarray]] = None) -> PDeltaResult`**
  Solve ``(K + K_sigma(U)) U = F`` with homogeneous constraints on ``fixed_dofs``.

### `onedfem.processing.common`

[src/onedfem/processing/common/__init__.py](../../src/onedfem/processing/common/__init__.py)

Shared processing utilities.

- **`init_stage_logger(class_name: str, job_results_dir: Optional[Union[str, Path]]) -> logging.Logger`**
  Return a logger with optional ``logs/<ClassName>.log`` under the job tree.

### `onedfem.processing.dynamic`

[src/onedfem/processing/dynamic/__init__.py](../../src/onedfem/processing/dynamic/__init__.py)

Deprecated compatibility shim for transient processing.

- **`assemble_global_system(*args, **kwargs)`**
- **`apply_boundary_conditions(*args, **kwargs)`**
- **`newmark_integrate(*args, **kwargs)`**

### `onedfem.processing.dynamic.operations`

[src/onedfem/processing/dynamic/operations/__init__.py](../../src/onedfem/processing/dynamic/operations/__init__.py)

Staged transient dynamics pipeline.

- **class `AssembleDynamicGlobalSystem(elements: List, element_stiffness_matrices, element_mass_matrices, total_dof: int, job_results_dir: Optional[Union[str, Path]] = None, element_damping_matrices=None)`**
- **class `ModifyDynamicGlobalSystem(fixed_dofs: Optional[Sequence[int]] = None, prescribed_displacements: Optional[dict] = None, job_results_dir: Optional[Union[str, Path]] = None)`**
- **class `IntegrateTransientSystem(t_grid: np.ndarray, force_func: Callable[[float], np.ndarray], job_results_dir: Optional[Union[str, Path]] = None)`**
  Wrap :func:`onedfem.processing.dynamic.time_integration.newmark_integrate`.

### `onedfem.processing.dynamic.operations._logging`

[src/onedfem/processing/dynamic/operations/_logging.py](../../src/onedfem/processing/dynamic/operations/_logging.py)

Re-export shared stage logger for transient operations.

- **`init_stage_logger(class_name: str, job_results_dir: Optional[Union[str, Path]]) -> logging.Logger`**
  Return a logger with optional ``logs/<ClassName>.log`` under the job tree.

### `onedfem.processing.eigen`

[src/onedfem/processing/eigen/__init__.py](../../src/onedfem/processing/eigen/__init__.py)

Eigenproblems and spectral assembly helpers (taxonomy §2).

- **`PENALTY`** (constant)
- **`_compute_local_global_dof_map(elements: List, total_dof: int) -> List[numpy.ndarray]`**
  Compute global DOF indices for each element (same convention as static).
- **`_scatter_element_matrix(mat: Union[numpy.ndarray, scipy.sparse._coo.coo_matrix], dof_map: numpy.ndarray) -> Tuple[numpy.ndarray, numpy.ndarray, numpy.ndarray]`**
  Return (rows, cols, data) in global indexing for one element matrix.
- **`assemble_global_matrices(elements: List, element_stiffness_matrices: Optional[List] = None, element_mass_matrices: Optional[List] = None, total_dof: Optional[int] = None, job_results_dir: Optional[str] = None) -> Tuple[scipy.sparse._csr.csr_matrix, scipy.sparse._csr.csr_matrix, List[numpy.ndarray]]`**
  Assemble global K and M for eigen / modal vibration analysis.
- **`apply_boundary_conditions(K_global: scipy.sparse._csr.csr_matrix, M_global: scipy.sparse._csr.csr_matrix, fixed_dofs: Union[Sequence[int], numpy.ndarray, NoneType] = None, prescribed_displacements: Optional[dict] = None) -> Tuple[scipy.sparse._csr.csr_matrix, scipy.sparse._csr.csr_matrix, numpy.ndarray]`**
  Apply BCs to both K and M (penalty method): fixed DOFs get zero row/col, diagonal K = penalty, diagonal M = 1.0.

### `onedfem.processing.eigen.metrics`

[src/onedfem/processing/eigen/metrics/__init__.py](../../src/onedfem/processing/eigen/metrics/__init__.py)

Native eigen modal metrics (no onedfem.processing.static dependency).

- **`global_uz_unit_pattern(total_dof: int, dof_per_node: int) -> np.ndarray`**
  Return a unit-norm vector activating uz on every node (zero if unsupported layout).
- **`modal_effective_mass_fraction_z(M_mod: csr_matrix, mode_shapes: np.ndarray, *, dof_per_node: int) -> np.ndarray | None`**
  Per mode *j*, scalar ((φ_jᵀ M r)²) / ((rᵀ M r)(φ_jᵀ M φ_j))) with *r* = global uz unit pattern.

### `onedfem.processing.eigen.rotating_frame`

[src/onedfem/processing/eigen/rotating_frame.py](../../src/onedfem/processing/eigen/rotating_frame.py)

Free vibration of a beam structure spinning at a constant rate, in the rotating frame.

- **class `RotatingFrameResult(rotation: Rotation, frequencies_hz: npt.NDArray[np.float64], modes: npt.NDArray[np.complex128], U_centrifugal: npt.NDArray[np.float64], f_centrifugal: npt.NDArray[np.float64], element_objects: List[RotatingFrameObject] = <factory>) -> None`**
  Frequencies (Hz, rotating frame), modes (total_dof, n) and the centrifugal state.
- **`assemble_rotating_frame(elements: Sequence, coordinates: npt.ArrayLike, connectivity: Sequence, rotation: Rotation, total_dof: int)`**
  Global ``f_c`` (n,), ``K_soft``, ``K_rot``, ``G`` (sparse) and the element objects.
- **`rotating_frame_vibration(elements: Sequence, coordinates, connectivity, K: sp.spmatrix, M: sp.spmatrix, bc_dofs, rotation: Rotation, num_modes: int, gyroscopic: bool = True, prestress: bool = True) -> RotatingFrameResult`**
  Rotating-frame frequencies and modes (see the module docstring).
- **`campbell(elements, coordinates, connectivity, K, M, bc_dofs, rotation: Rotation, omegas: Sequence[float], num_modes: int, gyroscopic: bool = True, prestress: bool = True) -> npt.NDArray[np.float64]`**
  Rotating-frame frequencies (Hz) at each speed: rows ``[omega, f_1, ..., f_n]``.

### `onedfem.processing.harmonic`

[src/onedfem/processing/harmonic/__init__.py](../../src/onedfem/processing/harmonic/__init__.py)

Frequency-domain dynamics (taxonomy §4): harmonic / steady-state response.

- **`frequency_grid_hz(f_min_hz: float, f_max_hz: float, n_points: int) -> np.ndarray`**
  Uniformly spaced frequencies in Hz (inclusive endpoints).
- **`harmonic_damping_matrix(M_mod: csr_matrix, K_mod: csr_matrix, zeta: float, omega_ref_rad: float, rayleigh_alpha: float, rayleigh_beta: float) -> csr_matrix`**
  Combined damping: mass-proportional reference scaling plus optional Rayleigh terms.
- **`harmonic_displacement_modal_superposition(omega_natural_rad: np.ndarray, mode_shapes: np.ndarray, F: np.ndarray, frequencies_hz: np.ndarray, zeta: float, *, zeta_per_mode: Optional[np.ndarray] = None) -> np.ndarray`**
  Steady-state response using classical viscous modal damping.
- **`mass_proportional_damping(M_mod: csr_matrix, zeta: float, omega_ref_rad: float) -> csr_matrix`**
  Mass-proportional damping: C = 2 * zeta * omega_ref * M.
- **`rayleigh_damping(M_mod: csr_matrix, K_mod: csr_matrix, alpha: float, beta: float) -> csr_matrix`**
  Rayleigh damping contribution: alpha * M + beta * K (alpha [1/s], beta [s]).
- **`solve_one_frequency(K_mod: csr_matrix, M_mod: csr_matrix, C_mod: csr_matrix, F_mod: np.ndarray, omega_rad: float, *, linear_solver: str = 'spsolve', sparsity_signature_cache: Optional[list] = None, splu_csc_workspace: Optional[dict] = None) -> np.ndarray`**
  Solve (-omega^2 M + i omega C + K) u = F for complex u.
- **`sweep_displacements(K_mod: csr_matrix, M_mod: csr_matrix, C_mod: csr_matrix, F_mod: np.ndarray, frequencies_hz: np.ndarray, *, u_prescribed: Optional[np.ndarray] = None, dof_unknown: Optional[np.ndarray] = None, parallel: bool = False, max_workers: Optional[int] = None, mp_damping_reference: str = 'geometric_mean', zeta_mp: float = 0.0, rayleigh_alpha_mp: float = 0.0, rayleigh_beta_mp: float = 0.0, linear_solver: str = 'spsolve', zeta_per_frequency: Optional[np.ndarray] = None, geometric_omega_ref_rad: Optional[float] = None) -> np.ndarray`**
  Column k is u(2*pi*f_k) for each frequency in *frequencies_hz*.
- **`undamped_natural_modes(K_mod: csr_matrix, M_mod: csr_matrix, num_modes: int, *, dense_threshold: int = 512) -> tuple[np.ndarray, np.ndarray]`**
  Lowest undamped natural frequencies (rad/s) and mode shapes (columns), M-normalized.

### `onedfem.processing.harmonic.operations`

[src/onedfem/processing/harmonic/operations/__init__.py](../../src/onedfem/processing/harmonic/operations/__init__.py)

Staged harmonic pipeline (structural assembly, BCs, damping, frequency solve).

- **class `AssembleHarmonicStructuralMatrices(elements: List, element_stiffness_matrices: List, element_mass_matrices: List, total_dof: int, job_results_dir: Optional[Union[str, Path]] = None)`**
- **class `AssembleHarmonicLoadVector(elements: List, element_stiffness_matrices: List, element_force_vectors: List, total_dof: int, job_results_dir: Optional[Union[str, Path]] = None)`**
- **class `ModifyHarmonicStructuralMatrices(prescribed_displacements: Optional[dict] = None, job_results_dir: Optional[Union[str, Path]] = None, fixed_dofs: Optional[Sequence[int]] = None)`**
- **class `BuildHarmonicDampingMatrix(job_results_dir: Optional[Union[str, Path]] = None)`**
- **class `HarmonicSweepConfig(use_modal: bool, num_modal: int, basis_dir: Optional[str], basis_job: Optional[str], zeta: float, use_parallel: bool, mp_ref: str, rayleigh_alpha: float, rayleigh_beta: float, lin_solver: str, zpf: Optional[np.ndarray], ft_tab: Any, zt_tab: Any) -> None`**
  HarmonicSweepConfig(use_modal: 'bool', num_modal: 'int', basis_dir: 'Optional[str]', basis_job: 'Optional[str]', zeta: 'float', use_parallel: 'bool', mp_ref: 'str', rayleigh_alpha: 'float', rayleigh_beta: 'float', lin_solver: 'str', zpf: 'Optional[np.ndarray]', ft_tab: 'Any', zt_tab: 'Any')
- **class `SolveHarmonicFrequencySweep(job_results_dir: Optional[Union[str, Path]] = None, resolve_job_path: Optional[Callable[[str], str]] = None)`**

### `onedfem.processing.harmonic.operations._logging`

[src/onedfem/processing/harmonic/operations/_logging.py](../../src/onedfem/processing/harmonic/operations/_logging.py)

Re-export shared stage logger for harmonic operations.

- **`init_stage_logger(class_name: str, job_results_dir: Optional[Union[str, Path]]) -> logging.Logger`**
  Return a logger with optional ``logs/<ClassName>.log`` under the job tree.

### `onedfem.processing.spectral`

[src/onedfem/processing/spectral/__init__.py](../../src/onedfem/processing/spectral/__init__.py)

Spectral analysis helpers (Sections 2 and 5): diagnostics and staged operations.

- **`log_spectral_constrained_dofs(bc_dofs, total_dof: int | None = None, job_results_dir=None)`**
  Log penalty / fixed DOF count after spectral matrix modification.
- **`log_spectral_diagnostics(K, M, job_results_dir=None)`**
  Log basic diagnostics for global K and M (shape, nnz, Frobenius norm of data).

### `onedfem.processing.spectral.operations`

[src/onedfem/processing/spectral/operations/__init__.py](../../src/onedfem/processing/spectral/operations/__init__.py)

Staged spectral pipeline (eigen vibration + linear buckling matrix steps).

- **class `PrepareSpectralLocalMatrices(element_stiffness_matrices: Optional[Sequence] = None, element_mass_matrices: Optional[Sequence] = None, job_results_dir: Optional[Union[str, Path]] = None)`**
  Validate and convert element stiffness / mass matrices to COO float64 lists, mirroring :class:`onedfem.processing.static.operations.preparation.PrepareLocalSystem` for the stiffness path (mass has no force analogue).
- **class `AssembleSpectralGlobalSystem(elements: List, element_stiffness_matrices: List, element_mass_matrices: List, total_dof: int, job_results_dir: Optional[Union[str, Path]] = None)`**
  Scatter element matrices to global ``K`` and ``M`` (delegates to :func:`onedfem.processing.eigen.assembly.assemble_global_matrices`), then log diagnostics.
- **class `ModifySpectralGlobalSystem(job_results_dir: Optional[Union[str, Path]] = None, fixed_dofs: Optional[Sequence[int]] = None, prescribed_displacements: Optional[dict] = None)`**
  Wrap :func:`onedfem.processing.eigen.boundary_conditions.apply_boundary_conditions`.
- **class `SolveGeneralizedEigenproblem(num_modes: int, context: str = 'eigen vibration', dense_threshold: int = 512, job_results_dir: Optional[Union[str, Path]] = None)`**
  Solve :math:`K x = \lambda M x` for the smallest-magnitude eigenvalues (delegates to :func:`onedfem.processing.eigen.smallest_generalized_eigenpairs`).
- **class `AssembleBucklingGeometricStiffness(elements: List, U_global: np.ndarray, total_dof: int, job_results_dir: Optional[Union[str, Path]] = None)`**
  Build global geometric stiffness from prestress displacements.
- **class `ModifyBucklingGlobalMatrices(job_results_dir: Optional[Union[str, Path]] = None, prescribed_displacements: Optional[dict] = None, fixed_dofs: Optional[Sequence[int]] = None)`**
  Apply buckling BCs to elastic K and geometric K_sigma.
- **class `SolveLinearBucklingEigenpairs(num_modes: int, job_results_dir: Optional[Union[str, Path]] = None)`**
  Linear buckling eigenproblem on (K, K_sigma) pencil.

### `onedfem.processing.static`

[src/onedfem/processing/static/__init__.py](../../src/onedfem/processing/static/__init__.py)

Static FEM processing: assembly, boundary conditions, condensation, solve, results.

- **class `PrepareLocalSystem(Ke_raw: Union[List[Union[numpy.ndarray, scipy.sparse._coo.coo_matrix]], numpy.ndarray], Fe_raw: Union[List[Union[numpy.ndarray, list]], numpy.ndarray], job_results_dir: Union[str, pathlib.Path, NoneType] = None)`**
  Validates and formats local element stiffness matrices (Ke) and force vectors (Fe) into the formats required by AssembleGlobalSystem: - Ke: scipy.sparse.coo_matrix with float64 precision - Fe: 1D numpy arrays with float64 precision
- **class `AssembleGlobalSystem(elements: List[object], element_stiffness_matrices: Optional[List[scipy.sparse._coo.coo_matrix]] = None, element_force_vectors: Optional[List[numpy.ndarray]] = None, total_dof: Optional[int] = None, job_results_dir: Optional[str] = None, symmetry_tol: float = 1e-08)`**
  Assembles the global stiffness matrix and force vector for a finite element system.
- **class `ModifyGlobalSystem(K_global: scipy.sparse._csr.csr_matrix, F_global: numpy.ndarray, local_global_dof_map: Sequence[numpy.ndarray], job_results_dir: Optional[str] = None, fixed_dofs: Union[Sequence[int], numpy.ndarray, NoneType] = None, prescribed_displacements: Optional[dict] = None)`**
  High-performance boundary condition applier with numerical stabilization and diagnostics.
- **class `CondenseModifiedSystem(K_mod: scipy.sparse._csr.csr_matrix, F_mod: numpy.ndarray, fixed_dofs: Union[Sequence[int], numpy.ndarray], local_global_dof_map: list[numpy.ndarray], job_results_dir: Optional[str] = None, base_tol: float = 1e-12)`**
  Static condensation system with advanced validation and adaptive numerics.
- **class `SolveCondensedSystem(K_cond: sp.spmatrix, F_cond: np.ndarray, solver_name: str, job_results_dir: str, *, preconditioner: str | None = 'auto', max_mem_gb: float = 10.0, tolerance: float = 1e-06, max_iterations: int = 1000, restart: int = 20, ilu_drop_tol: float = 1e-06, ilu_fill_factor: float = 1.0, disable_scaling: bool = False)`**
  Solve `K_cond · U_cond = F_cond`.
- **class `ReconstructGlobalSystem(active_dofs: numpy.ndarray, U_cond: numpy.ndarray, total_dofs: int, job_results_dir: pathlib.Path, *, local_global_dof_map: Sequence[numpy.ndarray], fixed_dofs: Optional[numpy.ndarray] = None, inactive_dofs: Optional[numpy.ndarray] = None)`**
  High-performance displacement reconstruction system with validation and diagnostics.
- **class `DisassembleGlobalSystem(*, U_global: numpy.ndarray, R_global: numpy.ndarray, R_residual: numpy.ndarray, local_global_dof_map: Sequence[numpy.ndarray], job_results_dir: str | pathlib.Path | None = None, F_e: Optional[Sequence[numpy.ndarray]] = None) -> None`**
  Reverse-assembles global results to each element and writes:

### `onedfem.processing.static.follower_loads`

[src/onedfem/processing/static/follower_loads.py](../../src/onedfem/processing/static/follower_loads.py)

Nodal follower loads: forces and moments that turn with the node they act on.

- **`parse_follower_load(path: str) -> npt.NDArray[np.float64]`**
  Rows ``[x, y, z, F_x, F_y, F_z, M_x, M_y, M_z]`` under a ``[Follower load]`` header.
- **class `FollowerLoads(rows: npt.ArrayLike, coordinates: npt.ArrayLike, dof_per_node: int, total_dof: int)`**
  Follower loads on nodes; ``forces(U)`` and ``stiffness(U)`` on the global DOFs.
- **`FOLLOWER_LOAD_FILE`** (constant)

### `onedfem.processing.transient`

[src/onedfem/processing/transient/__init__.py](../../src/onedfem/processing/transient/__init__.py)

Transient analysis: canonical Section 3 package name.

- **`assemble_global_system(elements: List, element_stiffness_matrices: Optional[List] = None, element_mass_matrices: Optional[List] = None, element_damping_matrices: Optional[List] = None, total_dof: Optional[int] = None, job_results_dir: Optional[str] = None) -> Tuple[scipy.sparse._csr.csr_matrix, scipy.sparse._csr.csr_matrix, Optional[scipy.sparse._csr.csr_matrix], List[numpy.ndarray]]`**
  Assemble global K, M, and optionally C for dynamic analysis.
- **`apply_boundary_conditions(K_global: scipy.sparse._csr.csr_matrix, M_global: scipy.sparse._csr.csr_matrix, C_global: Optional[scipy.sparse._csr.csr_matrix] = None, fixed_dofs: Union[Sequence[int], numpy.ndarray, NoneType] = None, prescribed_displacements: Optional[dict] = None) -> Tuple[scipy.sparse._csr.csr_matrix, scipy.sparse._csr.csr_matrix, Optional[scipy.sparse._csr.csr_matrix], numpy.ndarray]`**
  Apply BCs to K, M, and optionally C (penalty method).
- **`newmark_integrate(K: scipy.sparse._csr.csr_matrix, M: scipy.sparse._csr.csr_matrix, C: Optional[scipy.sparse._csr.csr_matrix], u0: numpy.ndarray, v0: numpy.ndarray, t_grid: numpy.ndarray, F_func: Callable[[float], numpy.ndarray], beta: float = 0.25, gamma: float = 0.5) -> Tuple[numpy.ndarray, numpy.ndarray, numpy.ndarray]`**
  Newmark-beta time integration: M u'' + C u' + K u = F(t).

### `onedfem.processing.transient.assembly`

[src/onedfem/processing/transient/assembly.py](../../src/onedfem/processing/transient/assembly.py)

- **`assemble_global_system(elements: List, element_stiffness_matrices: Optional[List] = None, element_mass_matrices: Optional[List] = None, element_damping_matrices: Optional[List] = None, total_dof: Optional[int] = None, job_results_dir: Optional[str] = None) -> Tuple[scipy.sparse._csr.csr_matrix, scipy.sparse._csr.csr_matrix, Optional[scipy.sparse._csr.csr_matrix], List[numpy.ndarray]]`**
  Assemble global K, M, and optionally C for dynamic analysis.
- **`assemble_global_force_vector(elements: List, element_force_vectors: Optional[Sequence], total_dof: int, *, job_results_dir: Union[str, pathlib.Path, NoneType] = None, local_global_dof_map: Optional[List[numpy.ndarray]] = None) -> numpy.ndarray`**
  Scatter element nodal loads ``F_e`` into a global vector (same pattern as static force assembly).

### `onedfem.processing.transient.boundary_conditions`

[src/onedfem/processing/transient/boundary_conditions.py](../../src/onedfem/processing/transient/boundary_conditions.py)

- **`apply_boundary_conditions(K_global: scipy.sparse._csr.csr_matrix, M_global: scipy.sparse._csr.csr_matrix, C_global: Optional[scipy.sparse._csr.csr_matrix] = None, fixed_dofs: Union[Sequence[int], numpy.ndarray, NoneType] = None, prescribed_displacements: Optional[dict] = None) -> Tuple[scipy.sparse._csr.csr_matrix, scipy.sparse._csr.csr_matrix, Optional[scipy.sparse._csr.csr_matrix], numpy.ndarray]`**
  Apply BCs to K, M, and optionally C (penalty method).

### `onedfem.processing.transient.diagnostics.transient_run_diagnostic`

[src/onedfem/processing/transient/diagnostics/transient_run_diagnostic.py](../../src/onedfem/processing/transient/diagnostics/transient_run_diagnostic.py)

- **`log_transient_modified_system(K_mod: csr_matrix, M_mod: csr_matrix, C_mod: Optional[csr_matrix], *, n_bc_dofs: int, job_results_dir: Optional[Union[str, Path]] = None) -> None`**
  Log Frobenius norms and BC counts after dynamic BC application.

### `onedfem.processing.transient.operations`

[src/onedfem/processing/transient/operations/__init__.py](../../src/onedfem/processing/transient/operations/__init__.py)

Staged transient dynamics pipeline.

- **class `AssembleTransientGlobalSystem(elements: List, element_stiffness_matrices, element_mass_matrices, total_dof: int, job_results_dir: Optional[Union[str, Path]] = None, element_damping_matrices=None)`**
- **class `ModifyTransientGlobalSystem(fixed_dofs: Optional[Sequence[int]] = None, prescribed_displacements: Optional[dict] = None, job_results_dir: Optional[Union[str, Path]] = None)`**
- **class `IntegrateTransientSystem(t_grid: np.ndarray, force_func: Callable[[float], np.ndarray], job_results_dir: Optional[Union[str, Path]] = None)`**
  Wrap :func:`onedfem.processing.dynamic.time_integration.newmark_integrate`.
- **class `AssembleDynamicGlobalSystem(elements: List, element_stiffness_matrices, element_mass_matrices, total_dof: int, job_results_dir: Optional[Union[str, Path]] = None, element_damping_matrices=None)`**
- **class `ModifyDynamicGlobalSystem(fixed_dofs: Optional[Sequence[int]] = None, prescribed_displacements: Optional[dict] = None, job_results_dir: Optional[Union[str, Path]] = None)`**

### `onedfem.processing.transient.operations._logging`

[src/onedfem/processing/transient/operations/_logging.py](../../src/onedfem/processing/transient/operations/_logging.py)

- **`init_stage_logger(class_name: str, job_results_dir: Optional[Union[str, Path]]) -> logging.Logger`**
  Return a logger with optional ``logs/<ClassName>.log`` under the job tree.

### `onedfem.processing.transient.operations.integrate_transient_system`

[src/onedfem/processing/transient/operations/integrate_transient_system.py](../../src/onedfem/processing/transient/operations/integrate_transient_system.py)

- **class `IntegrateTransientSystem(t_grid: np.ndarray, force_func: Callable[[float], np.ndarray], job_results_dir: Optional[Union[str, Path]] = None)`**
  Wrap :func:`onedfem.processing.dynamic.time_integration.newmark_integrate`.

### `onedfem.processing.transient.time_integration`

[src/onedfem/processing/transient/time_integration.py](../../src/onedfem/processing/transient/time_integration.py)

- **`newmark_integrate(K: scipy.sparse._csr.csr_matrix, M: scipy.sparse._csr.csr_matrix, C: Optional[scipy.sparse._csr.csr_matrix], u0: numpy.ndarray, v0: numpy.ndarray, t_grid: numpy.ndarray, F_func: Callable[[float], numpy.ndarray], beta: float = 0.25, gamma: float = 0.5) -> Tuple[numpy.ndarray, numpy.ndarray, numpy.ndarray]`**
  Newmark-beta time integration: M u'' + C u' + K u = F(t).

### `onedfem.processing.transient.transient_forcing`

[src/onedfem/processing/transient/transient_forcing.py](../../src/onedfem/processing/transient/transient_forcing.py)

- **`build_transient_force_func(F_ref: np.ndarray, dyn_config: dict, *, total_dof: int, job_dir: Optional[str], end_time: float) -> Callable[[float], np.ndarray]`**
  Build ``F(t)`` from assembled reference load ``F_ref`` and merged transient settings.

### `onedfem.simulation_runner.buckling`

[src/onedfem/simulation_runner/buckling/__init__.py](../../src/onedfem/simulation_runner/buckling/__init__.py)

§5 Buckling / stability runners.

- **class `BucklingSimulationRunner(*args, **kwargs)`**
  Deprecated alias for :class:`LinearBucklingSimulationRunner`.
- **class `LinearBucklingSimulationRunner(settings, job_name)`**
  Linear buckling eigenproblem about a prestressed state (§5 linearized theory).
- **class `NonlinearBucklingSimulationRunner(settings: dict, job_name: str)`**
  Load-control nonlinear buckling / continuation MVP.

### `onedfem.simulation_runner.eigen`

[src/onedfem/simulation_runner/eigen/__init__.py](../../src/onedfem/simulation_runner/eigen/__init__.py)

§2 Eigen / spectral analysis runners.

- **class `EigenSimulationRunner(settings, job_name)`**
  Natural frequencies and mode shapes (same physics pipeline as legacy modal vibration).

### `onedfem.simulation_runner.harmonic`

[src/onedfem/simulation_runner/harmonic/__init__.py](../../src/onedfem/simulation_runner/harmonic/__init__.py)

§4 Frequency-domain harmonic runners (stub).

- **class `HarmonicSimulationRunner(settings: dict, job_name: str)`**
  Harmonic / frequency-response analysis using direct frequency sweep.

### `onedfem.simulation_runner.spectral`

[src/onedfem/simulation_runner/spectral/__init__.py](../../src/onedfem/simulation_runner/spectral/__init__.py)

Shared spectral analysis backend for eigen vibration and linear buckling (analysis sections 2 and 5).

- **class `VibrationBucklingBackend(settings, job_name)`**
  Shared implementation for Section 2 eigen vibration and Section 5 linear buckling (undamped :math:`K` / :math:`M` and buckling pencil). Subclass with a concrete ``run()`` in :class:`~onedfem.simulation_runner.eigen.eigen_simulation.EigenSimulationRunner` or :class:`~onedfem.simulation_runner.buckling.buckling_simulation.LinearBucklingSimulationRunner`.
- **`log_spectral_constrained_dofs(bc_dofs, total_dof: int | None = None, job_results_dir=None)`**
  Log penalty / fixed DOF count after spectral matrix modification.
- **`log_spectral_diagnostics(K, M, job_results_dir=None)`**
  Log basic diagnostics for global K and M (shape, nnz, Frobenius norm of data).

### `onedfem.simulation_runner.spectral.spectral_diagnostics`

[src/onedfem/simulation_runner/spectral/spectral_diagnostics.py](../../src/onedfem/simulation_runner/spectral/spectral_diagnostics.py)

Re-export spectral diagnostics from :mod:`onedfem.processing.spectral` (single source of truth).

- **`log_spectral_constrained_dofs(bc_dofs, total_dof: int | None = None, job_results_dir=None)`**
  Log penalty / fixed DOF count after spectral matrix modification.
- **`log_spectral_diagnostics(K, M, job_results_dir=None)`**
  Log basic diagnostics for global K and M (shape, nnz, Frobenius norm of data).

### `onedfem.simulation_runner.transient`

[src/onedfem/simulation_runner/transient/__init__.py](../../src/onedfem/simulation_runner/transient/__init__.py)

§3 Transient dynamics runners.

- **class `DynamicSimulationRunner(*args, **kwargs)`**
  Deprecated alias for :class:`TransientSimulationRunner`.
- **class `TransientSimulationRunner(settings: dict, job_name: str)`**
  Transient finite element analysis: assemble K/M/C, apply BCs, Newmark integrate.

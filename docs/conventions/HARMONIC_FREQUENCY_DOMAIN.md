# Harmonic analysis (§4) — frequency-domain design

**Status:** Implemented in [`simulation_runner/harmonic/harmonic_simulation.py`](../../simulation_runner/harmonic/harmonic_simulation.py) with kernels in [`processing/harmonic/frequency_response.py`](../../processing/harmonic/frequency_response.py). This note fixes notation, v1 scope, and extension points.

## Governing semi-discrete equation

For angular frequency \(\omega\) (rad/s) and consistent undamped structural matrices **K**, **M** assembled in the global frame:

\[
\left(-\omega^2 \mathbf{M} + i \omega \mathbf{C} + \mathbf{K}\right) \mathbf{u}(\omega) = \mathbf{f}(\omega).
\]

Loads are taken as **real** harmonic amplitudes \(\mathbf{f}\) (same pattern as linear static assembly); the unknown \(\mathbf{u}\) is **complex** to capture magnitude and phase.

Sign conventions match standard vibration textbooks and the undamped eigenproblem used in §2 (**K**, **M** from [`processing.eigen.assembly`](../../processing/eigen/assembly.py)).

## Boundary conditions

After global assembly, **K** and **M** receive the same **penalty** boundary treatment as §2 eigen ([`processing.eigen.boundary_conditions.apply_boundary_conditions`](../../processing/eigen/boundary_conditions.py)), including optional **`prescribed_displacement.txt`** wiring from **`run_job`**: fixed / prescribed-zero DOFs enforce large stiffness on **K** and identity scaling on **M**. The assembled load vector **F** has entries zeroed at constrained DOFs.

**Non-zero prescribed displacement amplitudes** (real values from the same parser used for static jobs) define **known harmonic motion** on those global DOFs; optional per-row phase (sixth column) and/or global **`[Harmonic] prescribed_motion_phase_rad`**. The solver partitions \(\mathbf{u} = \mathbf{u}^{(p)} + \mathbf{u}^{(f)}\) with \(\mathbf{u}^{(p)}\) fixed on prescribed-motion DOFs and solves the reduced system on the remaining free DOFs. **`point_load.txt`** supports an optional **10th column** **`phase_rad`** per row for complex equivalent nodal loads; **`[Harmonic] load_phase_rad`** still applies an extra **global** phase on the assembled load vector.

## Damping model

**Default:** mass-proportional damping

\[
\mathbf{C}_{\mathrm{mp}} = 2 \zeta \, \omega_{\mathrm{ref}} \, \mathbf{M}_{\mathrm{mod}},
\]

where \(\zeta\) is **`modal_damping_ratio`** from **`[Harmonic]`**, and \(\omega_{\mathrm{ref}} = 2\pi f_{\mathrm{ref}}\) with \(f_{\mathrm{ref}} = \sqrt{f_{\min} f_{\max}}\) (Hz).

**Optional Rayleigh** terms (after BCs, same layout as **K** / **M**):

\[
\mathbf{C} = \mathbf{C}_{\mathrm{mp}} + \alpha \, \mathbf{M}_{\mathrm{mod}} + \beta \, \mathbf{K}_{\mathrm{mod}},
\]

with **`rayleigh_alpha`** \(\alpha\) (1/s) and **`rayleigh_beta`** \(\beta\) (s). Defaults \(\alpha=\beta=0\) preserve legacy behaviour.

**Mass-proportional reference:** **`mp_damping_reference`** selects how \(\omega_{\mathrm{ref}}\) enters \(\mathbf{C}_{\mathrm{mp}}\):

- **`geometric_mean`** (default): \(\omega_{\mathrm{ref}} = 2\pi\sqrt{f_{\min} f_{\max}}\) — one **C** for the whole sweep (legacy behaviour).
- **`current_frequency`**: replace \(\omega_{\mathrm{ref}}\) with the **current** sample \(\omega\) when forming \(\mathbf{C}_{\mathrm{mp}}\) at each frequency (Rayleigh terms unchanged). Implemented by rebuilding **C** each sample inside [`sweep_displacements`](../../processing/harmonic/frequency_response.py).

**Optional modal superposition path:** **`use_modal_superposition`** — solve undamped modes \((\omega_r, \boldsymbol{\phi}_r)\) from **K**, **M** (same BCs), then build \(\mathbf{u}(\omega)\) with diagonal **viscous modal damping** using **`modal_damping_ratio`** \(\zeta\) on each modal oscillator ([`modal_superposition.py`](../../processing/harmonic/modal_superposition.py)), unless a **`damping_zeta_table_file`** is set (then \(\zeta_r\) is interpolated at each mode’s natural frequency in Hz). Incompatible with non-zero prescribed-motion partitioning in the current implementation (use the direct sweep instead). Rayleigh terms are **not** applied on this path (a warning is logged; use the direct sweep for full **C**).

**Prescribed motion phase:** optional **sixth numeric column** on each row of **`prescribed_displacement.txt`** after **`type`** gives per-DOF phase (rad), added to the global **`prescribed_motion_phase_rad`** on motion DOFs (when the column is omitted, only the global scalar applies).

**Frequency-dependent ζ:** **`[Harmonic] damping_zeta_table_file`** — two-column text (`frequency_hz`, `zeta`); values are interpolated for each sweep sample (direct path) or each retained mode (modal path). Job-relative paths resolve against the job folder (**`job_dir`** wiring from **`run_job`**).

**§2 mode reuse:** **`harmonic_modal_basis_dir`** points at a **`modal_results`** directory from a prior **`[Type] eigen`** job; optional **`harmonic_modal_basis_job_name`** selects **`{name}_frequencies.txt`** / **`{name}_mode_shapes.txt`**. When set, harmonic skips **`eigsh`** for the modal superposition basis.

**Sparse structure note:** SciPy’s high-level **`splu`** does not expose reuse of symbolic factorization when only matrix **values** change with \(\omega\). The runner can still reuse a **CSC data buffer** for **`splu`** (default on; set **`FEM_HARMONIC_SPLU_CSC_BUFFER=0`** to disable). Set **`FEM_HARMONIC_VERIFY_A_PATTERN=1`** to assert a fixed CSC pattern for **A(ω)** across the sweep (debug). Throughput: prefer **`parallel_frequency_sweep`**; optional script **`benchmarks/harmonic_sparse_microbench.py`** compares backends locally (**not** CI).

**Profiling gate (direct sweep):** Run the micro-benchmark on a mesh-sized surrogate or instrument a real job; if linear solves dominate wall time after enabling **`parallel_frequency_sweep`**, consider experimental solver backends only behind explicit settings and parity tests — high-level **`splu`** still performs a full numeric factorization each \(\omega\) unless a lower-level API or external package exposes symbolic reuse safely.

## Job inputs

Parser keys (see [`simulation_runner/harmonic/README.md`](../../simulation_runner/harmonic/README.md)) populate **`simulation_settings["harmonic"]`**:

| Key | Role |
|-----|------|
| `frequency_min_hz` | Lower sweep bound (Hz). |
| `frequency_max_hz` | Upper sweep bound (Hz). |
| `num_frequency_points` | Number of uniformly spaced samples (inclusive endpoints). |
| `modal_damping_ratio` | \(\zeta\) for mass-proportional **C**. |
| `rayleigh_alpha` | \(\alpha\) in \(\mathbf{C}\) (alias `rayleigh_m`). |
| `rayleigh_beta` | \(\beta\) in \(\mathbf{C}\) (alias `rayleigh_k`). |
| `load_phase_rad` | Uniform phase (rad) on the assembled load vector (alias `load_phase`). |
| `parallel_frequency_sweep` | If `true`, frequency samples use a concurrent thread pool. |
| `mp_damping_reference` | `geometric_mean` or `current_frequency` — see *Damping model*. |
| `use_modal_superposition` | If `true`, use undamped modes + modal harmonic expansion instead of assembling **C** once (see *Damping model*). |
| `modal_superposition_num_modes` | Number of modes retained for the modal path. |
| `prescribed_motion_phase_rad` | Phase (rad) on prescribed-motion DOFs (alias `prescribed_motion_phase`). |
| `harmonic_linear_solver` | `spsolve` (default) or `splu` — backend for each sparse complex solve (primarily diagnostic/benchmark). |
| `damping_zeta_table_file` | Optional path (job-relative or absolute) to ζ(**f**) table for mass-proportional damping per frequency / per mode. |
| `harmonic_modal_basis_dir` | Optional path to **`modal_results`** from a §2 eigen job (mode reuse). |
| `harmonic_modal_basis_job_name` | Optional eigen job name prefix for **`{name}_frequencies.txt`** / **`{name}_mode_shapes.txt`**. |

If any of these are omitted, the runner applies documented defaults (see `effective_harmonic_config` in code).

Loads: optional **10th column** on each **`point_load.txt`** row — **`phase_rad`** — phases that row’s equivalent nodal contribution \(\exp(i\phi)\) (complex **`F_e`** assembly).

**`[PostProcessing]`:** When **`run_secondary_tertiary_harmonic`** is **true**, secondary/tertiary snapshots use **`real(U[:, k])`**, **`imag(U[:, k])`**, or **both** per **`harmonic_secondary_tertiary_displacement_component`** (`real` \| `imag` \| `both`; default **`real`**). Column indices **`k`** come from **`harmonic_frequency_index`** (default `0`), **`harmonic_secondary_tertiary_frequency_indices`**, or **`harmonic_secondary_tertiary_all_frequencies`** — see **`RESULTS_DESIGN.md`** for precedence.

## Outputs

Under **`primary_results/harmonic_results/`**:

- **`{job_name}_frequencies_hz.txt`** — column of sampled frequencies.
- **`{job_name}_displacement_real.txt`** and **`..._imag.txt`** — DOF-by-frequency matrices (rows = global DOF index, columns = frequency index); complex response \(\mathbf{u}(\omega)\).
- **`{job_name}_displacement_abs.txt`** and **`{job_name}_displacement_phase_rad.txt`** — magnitude and phase (radians), same layout.

Optional secondary/tertiary outputs follow modal/dynamic when **`run_secondary_tertiary_harmonic`** is enabled (snapshot driven by **`harmonic_secondary_tertiary_displacement_component`**; **`harmonic_secondary_tertiary_all_frequencies`** exports per column).

### Primary output schema (v1)

- **`logs/primary_artifacts.json`** (when the job completes) lists the five primary text files above under **`artifacts`** with paths **relative to** the job results root. Field **`schema_version`** is **`1.0`**; **`family`** is **`harmonic`**. Downstream tools should treat column order as identical to the frequency vector in **`frequencies_hz`**.
- Helpers for a future **native complex stress** path live in [`processing/harmonic/complex_displacement_recovery.py`](../../processing/harmonic/complex_displacement_recovery.py) (`static_recovery_pair_from_complex_column`).

## Roadmap (next)

Status of the follow-on track (after the ζ tables, mode reuse, phased loads/BCs, selective post indices, **`FEM_HARMONIC_VERIFY_A_PATTERN`**, and modal dense-fallback work already in tree):

| Priority | Track | Notes |
|----------|--------|--------|
| 1 | **Low-level sparse reuse** | **Partial —** optional CSC shell reuse for **`splu`**: set **`FEM_HARMONIC_SPLU_CSC_BUFFER=0`** to disable. Each frequency still performs a full numeric factorization; this path cuts Python alloc / CSR→CSC churn when **A(ω)** keeps a fixed sparsity pattern. True SuperLU symbolic-only reuse is not exposed in SciPy’s high-level API. |
| 2 | **Stress harmonics in post** | **Partial —** **`[PostProcessing] harmonic_secondary_tertiary_displacement_component`** = **`real`** (default) \| **`imag`** \| **`both`** drives secondary/tertiary snapshots from **`real(U)`** / **`imag(U)`** (see **`RESULTS_DESIGN.md`**). |
| 3 | **Large-n modal eigen** | **Partial —** shared **`processing/eigen/smallest_generalized_eigenpairs.py`** (**SM** \(\rightarrow\) shift-invert \(\rightarrow\) dense fallback) powers harmonic **`undamped_natural_modes`** and **`VibrationBucklingBackend.solve_modal_vibration`**. |

## Relation to Section 2 and Section 3

- **Section 2 eigen:** Undamped modes diagnose resonance; harmonic sweep gives forced response at arbitrary \(\omega\).
- **Section 3 transient:** Time-domain Newmark; harmonic is **frequency-domain** steady-state sinusoidal response (different physics hook).

Use **direct frequency sweep** when full **M**, **C**, **K** matrices are required; use **`use_modal_superposition`** when a truncated modal basis is acceptable for performance (same BC metric as Section 2 assembly).

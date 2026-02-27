# Abaqus validation visualisers

This package **validates** the in-house FEM solver by running the same problem in **Abaqus** and comparing results. Unlike **verification_visualisers** (which compare FEM to analytical solutions like Roark), validation compares FEM to a second solver (Abaqus).

## Scope

- **Elements**: 2-node 3D Euler-Bernoulli (`EulerBernoulliBeamElement3D` → Abaqus B33) and 2-node 3D Timoshenko (`TimoshenkoBeamElement3D` → Abaqus B31).
- **Comparisons**: Deformation DOFs (u_y, θ_z and optionally all 6), GCI–Richardson (tip deflection/rotation vs mesh), and section forces (shear force diagram, bending moment diagram).

## Environment

- **abqpy**: Validation runs the generated script with the **project Python**; [abqpy](https://pypi.org/project/abqpy/) provides type hints and reimplements `mdb.saveAs()` to launch Abaqus. Install with: `pip install abqpy==2021.*` (match your Abaqus version). Optional dependency list: `post_processing/validation_visualisers/requirements-validation.txt`.
- **Abaqus**: Required only to *run* Abaqus and *extract* results from the ODB. You need an Abaqus installation with **CAE**; scripts are run with Abaqus’s bundled Python, e.g. `abaqus cae noGUI=<script>`. The run wrapper uses `C:\SIMULIA\CAE` by default (or env **ABAQUS_CAE_ROOT**); else `abaqus` must be on PATH. The runner sets **ABAQUS_BAT_PATH** from config so abqpy finds the launcher.
- **Comparison scripts** (deformation, GCI–Richardson, section forces) use only the project’s standard stack (numpy, matplotlib, existing parsers) and do **not** import Abaqus. They can run without Abaqus installed if Abaqus result CSVs are already present.

## Layout

- **abaqus/**: Input translation (job dir → Abaqus Python script), run wrapper, ODB → CSV extraction.
- **abaqus_results/**: Abaqus output per job (`job_XXXX_nX/`): the two CSVs are **U_global.csv**, then **section_forces.csv**; plus `rotation_source.txt` (contents **ODB** when rotation is read from the ODB and sign-flipped to match FEM, or **none** when rotation is not in the ODB and is written as zero), and optionally the run’s `.inp`, `.odb`, `.sta` (status), and `.msg` (messages) copied into the job dir so each result directory is self-contained and debuggable. Not committed (gitignore).
- **deformation/**: FEM vs Abaqus u_y and θ_z profile comparison. **Key outputs:** FEM n128 vs **Abaqus n500** (converged reference), one plot per base job. Writes to **deformation/deformation_plots/**.
- **section_forces/**: Section forces comparison (SFD/BMD: Vy, Mz). **Key outputs:** FEM n128 vs **Abaqus n500** (converged reference), one plot per base job. Writes to **section_forces/section_forces_plots/**.
- **grid_convergence_study/**: Convergence across meshes (GCI/Richardson, largest-mesh review, LaTeX table). Uses 3-grid GCI (n32, 64, 128) and **Abaqus n500** as reference; convergence behaviour across the 6 meshes (n4–n128) lives here. Writes to **grid_convergence_study/gci_tables/**.
- **output/**: Reserved for **review_abaqus_results.py** only (abaqus_results_review.csv, abaqus_performance_summary.md, errors log). Comparison plots and GCI/review CSVs are under the subtrees above.

## How Abaqus results are generated

Abaqus results in this project are **always generated using the Python Abaqus package** (abqpy). The pipeline is:

1. **Job input** (from `jobs/job_XXXX_nN/`) is translated into an Abaqus CAE Python script (`abaqus/generated/run_<job>.py`) by `job_to_abaqus_script.py`.
2. The script is **run with project Python**; [abqpy](https://pypi.org/project/abqpy/) provides the Abaqus scripting API and reimplements `mdb.saveAs()` so that it **launches Abaqus** (e.g. `abaqus cae noGUI=<script>`).
3. Abaqus executes the script: builds the model, runs the analysis, exports the ODB to CSV (U_global.csv, section_forces.csv, etc.) and writes them into `abaqus_results/job_XXXX_nN/`.

There is no separate .inp-only or command-line `abaqus job=... input=...` path; the **only** way to produce Abaqus results here is via the generated CAE script and the Python Abaqus package. Use `run_abaqus_cae.py` for one or more jobs, or `run_all_abaqus_jobs.py` (optionally with `--from-results`) to generate scripts and run Abaqus for many jobs at once.

## How to run

1. **Generate Abaqus script** for a job (e.g. `job_0000_n8`):
   ```bash
   python post_processing/validation_visualisers/abaqus/job_to_abaqus_script.py --job-dir jobs/job_0000_n8
   ```
   This writes a script under `abaqus/generated/`.

2. **Run Abaqus** (uses project Python + abqpy; **ABAQUS_BAT_PATH** is set from `C:\SIMULIA\CAE` or **ABAQUS_CAE_ROOT**):
   ```bash
   python post_processing/validation_visualisers/abaqus/run_abaqus_cae.py --job job_0000_n8
   ```
   Run several jobs by repeating `--job`:
   ```bash
   python post_processing/validation_visualisers/abaqus/run_abaqus_cae.py --job job_0000_n8 --job job_0005_n16
   ```
   This runs the generated script with project Python; abqpy's `saveAs()` launches Abaqus, which runs the script again to build the model, run the job, and export ODB → CSV into `abaqus_results/job_XXXX_nX/`.

3. **Run comparison** (after FEM n128 and Abaqus n500 results exist). Key plots compare **FEM n128 vs Abaqus n500** (one per base job):
   ```bash
   python post_processing/validation_visualisers/deformation/deformation_comparison.py
   python post_processing/validation_visualisers/section_forces/section_forces_comparison.py
   python post_processing/validation_visualisers/grid_convergence_study/gci_richardson_abaqus_report.py
   ```
   Or run all validation visualisers:
   ```bash
   python post_processing/validation_visualisers/run_all_validation_visualisers.py
   ```
   Outputs: `deformation/deformation_plots/`, `section_forces/section_forces_plots/`, `grid_convergence_study/gci_tables/`.

**Plots show only FEM (no Abaqus curves)?** Comparison scripts look for Abaqus **n500** reference CSVs in `abaqus_results/job_XXXX_n500/`. If those folders are empty or missing, plots will show "FEM only". Use the **n500 reference batch** (see below) to generate all 12 Abaqus n500 results, then re-run the comparison scripts.

**GCI–Richardson vs Abaqus:** To build the deflection/rotation table, run `grid_convergence_study/gci_richardson_abaqus_report.py`. It requires FEM results at n=32, 64, 128 and **Abaqus results at n=500** (converged reference) for jobs 0,1,2,5,6,7. Use the n500 reference batch (see below) to generate Abaqus n500 results, then run the GCI report and LaTeX table. Rotation comparison requires UR in the ODB; the generated script prints which field output was requested and whether the ODB contains UR (see Table caveats below).

**Table caveats:** (1) **Tip rotation:** Abaqus reference is zero unless the ODB contains rotational DOF (UR). The Abaqus script template now requests `U` and `UR`; re-generate run scripts, re-run Abaqus, and re-extract ODB→CSV to get non-zero rotation comparison. When you run the generated script, it prints **Field output requested: U, UR, SF** or **SF only (U/UR request failed)** and **ODB has UR: True/False** after the job completes. To confirm that U and UR were correctly extracted and written, check `run_log.txt` for the line **U_global.csv: U and UR read from ODB and written (rotation sign flipped to match FEM).** (or, if UR was missing from the ODB, **U_global.csv: U read from ODB; rotation not in ODB, written as zero.**). If UR is still missing, check the Abaqus job output (e.g. `.sta` or `.msg` in the Abaqus working directory) for field-output or UR-related errors; some Abaqus or beam element setups may not output UR. When UR is not in the ODB, rotation is written as zero and `rotation_source.txt` is **none**. (2) **Distributed loads (Triangular, Parabolic):** If Abaqus deflection appears zero for jobs 6 or 7, the reference is missing or tip node ordering differs—only point-load and UDL rows then have a meaningful error vs Abaqus.

**Abaqus n500 reference batch:** Deformation and section-forces comparisons use **Abaqus at n=500** as the converged reference (commercial code); FEM is compared at n128. To create the full n500 result set:

1. **Create job dirs** (from project root). The mesh variant scripts already include n500 in `VARIANT_NS`; run them so `jobs/job_0000_n500` … `jobs/job_0011_n500` exist:
   ```bash
   python pre_processing/mesh_library/create_point_load_mesh_variants.py
   python pre_processing/mesh_library/create_distributed_mesh_variants.py
   python pre_processing/mesh_library/create_timoshenko_mesh_variants.py
   ```
2. **Run Abaqus for the n500 batch** (generates scripts and runs Abaqus for all 12 jobs):
   ```bash
   python post_processing/validation_visualisers/run_all_abaqus_jobs.py --n500-only
   ```
   Use `--n500-only --dry-run` to list jobs without running; `--n500-only --script-only` to generate only CAE scripts (no Abaqus license required). Results go to `abaqus_results/job_XXXX_n500/`.

**Batch validation (optional):** For regression runs when Abaqus is available, use [run_batch_validation.py](post_processing/validation_visualisers/run_batch_validation.py): it runs Abaqus for a fixed set of jobs (default `job_0000_n128`, `job_0005_n128`), runs the comparison scripts, and checks that output files exist. Use `--n500-reference` to run the full Abaqus n500 batch (12 jobs) then comparisons and checks. Example: `python post_processing/validation_visualisers/run_batch_validation.py --n500-reference` or `--compare-only` to skip Abaqus and only run comparisons + checks.

**Re-run all Abaqus results:** To regenerate and run Abaqus for every validation job (all `job_XXXX_nN` under `jobs/`), see [PLAN_RERUN_ABAQUS_RESULTS.md](post_processing/validation_visualisers/PLAN_RERUN_ABAQUS_RESULTS.md). Use: `python post_processing/validation_visualisers/run_all_abaqus_jobs.py` (optional: `--dry-run`, `--jobs job_0000_n8 ...`, `--no-regenerate`).

**Generate Abaqus files from post_processing/results:** To generate Abaqus scripts (and optionally run Abaqus) only for jobs that have at least one FEM result directory under `post_processing/results`, use `--from-results`. Jobs are discovered from timestamped result dir names (`job_XXXX_nN_<timestamp>_pid...`); only jobs that also exist under `jobs/` are run (others are skipped with a message). **Run these commands from the project root** (the `fem_model` directory), not from inside `validation_visualisers`:
```bash
cd path\to\fem_model
python post_processing/validation_visualisers/run_all_abaqus_jobs.py --from-results
python post_processing/validation_visualisers/run_all_abaqus_jobs.py --from-results --dry-run
```
Use `--script-only` to generate only the Abaqus CAE scripts under `abaqus/generated/` without running Abaqus (no Abaqus license required):
```bash
python post_processing/validation_visualisers/run_all_abaqus_jobs.py --from-results --script-only
```

## Job coverage

| Job(s) | Element | Load | Use |
|--------|--------|------|-----|
| job_0000, 0001, 0002 | Euler-Bernoulli (B33) | Point (end, mid, quarter) | Deformation, SFD/BMD |
| job_0003 | Timoshenko (B31) | End point load | Deformation, SFD/BMD |
| job_0005, 0006, 0007 | Euler-Bernoulli (B33) | UDL, triangular, parabolic | Deformation, convergence, SFD/BMD |

**Distributed loads:** When `distributed_load.txt` is present, the script uses **LineLoad** (UDL) if Fy is constant (job_0005), or **equivalent nodal forces** (ConcentratedForce per node from lumped distributed values) for triangular/parabolic (job_0006, job_0007) so Abaqus matches the FEM load distribution.

Mesh variants (e.g. job_0000_n4, n8, n16, n32, n64, n128) are used for GCI–Richardson and comparisons. See `jobs/README_JOBS.md` for mesh variant generation.

## End-to-end validation checklist

To confirm the full pipeline on your machine:

1. Run: `python post_processing/validation_visualisers/abaqus/run_abaqus_cae.py --job job_0000_n8` (abqpy and Abaqus installed; `ABAQUS_BAT_PATH` from `C:\SIMULIA\CAE` or **ABAQUS_CAE_ROOT**).
2. If the generated script fails in Abaqus, apply version fixes from the "Abaqus version notes" section in [job_to_abaqus_script.py](post_processing/validation_visualisers/abaqus/job_to_abaqus_script.py) or document workarounds here.
3. Confirm outputs exist: `abaqus_results/job_0000_n8/U_global.csv` and `section_forces.csv`.
4. After FEM results exist for the same job, run `python post_processing/validation_visualisers/run_all_validation_visualisers.py` and check overlay plots in `deformation/deformation_plots/` and `section_forces/section_forces_plots/`, and GCI/review CSVs in `grid_convergence_study/gci_tables/`.
5. If you use a different Abaqus version, add a line under "Abaqus version notes" (e.g. *Tested with Abaqus 2024.*).

## Reviewing Abaqus result directories and performance

To audit result directories and solver performance (completion status, run time, .msg warnings/errors), run:

```bash
python post_processing/validation_visualisers/abaqus/review_abaqus_results.py
```

Options: `--expected` to compare against all jobs under `jobs/` and list jobs with no result dir; `--output <path>` for the CSV report; `--md <path>` for the Markdown summary; `--no-log-errors` to skip the errors/inconsistencies log. Default outputs: `validation_visualisers/output/abaqus_results_review.csv`, `validation_visualisers/output/abaqus_performance_summary.md`, and `validation_visualisers/output/abaqus_results_errors_and_inconsistencies.log` (run_log and .msg error/warning lines plus inconsistencies such as missing files, bad U_global header, or rotation_source vs ODB UR mismatch). The report covers file presence (U_global.csv, section_forces.csv, .inp, .odb, .sta, .msg, etc.), status from `.sta` (COMPLETED/ABORTED), total time from `.sta`, and error/warning counts from `.msg`. Validation performance (FEM vs Abaqus agreement) is produced by the comparison scripts; see `deformation/deformation_plots/`, `section_forces/section_forces_plots/`, and `grid_convergence_study/gci_tables/` for those results.

## Axis and section convention

- **Beam axis**: x along the wire; y, z are cross-section axes. Our section uses A, I_y, I_z, J_t; Abaqus Beam General Section uses I11, I22, I12, J. Mapping: I_y → I11, I_z → I22, I12 = 0, J_t → J. Document any sign differences in plots if needed.

## Abaqus version notes

The generated script uses the Abaqus Scripting API. If it fails on your Abaqus version, adjust as follows:

- **BeamSection**: Some versions use `BeamGeneralSection` or different parameter names (e.g. `crossSectionArea` vs `area`). Consult the Abaqus Scripting Reference for your release.
- **WirePolyLine**: Ensure `points` is a tuple of 3-tuples `((x,y,z), ...)`; `mergeType` and `meshable` may differ.
- **Mesh order**: Use seed → setElementType → generateMesh; the region must be the correct edge set (e.g. `part.Set(edges=part.edges, name="BeamSet")`).
- **DisplacementBC**: Encastre at node 1 may need to be created in the Initial step with all values set in one call (e.g. `u1=0, u2=0, u3=0, ur1=0, ur2=0, ur3=0`) instead of `setValuesInStep` in a later step.
- **ODB**: Newer Abaqus may use `import odb` instead of `odbAccess`; the script tries both when exporting results.
- **Field output request (U, UR, SF, SM):** The script requests **U**, **UR** (deflection/rotation), **SF** (section forces: N, Vy, Vz), and **SM** (section moments: T, My, Mz), analogous to requesting both U and UR. On Abaqus 2021 both `model.FieldOutputRequest` and `step.FieldOutputRequest` (StaticStep) may be missing. The script tries model first, then step; if both fail, Abaqus default field output often still writes U and UR to the ODB. For section output, fallbacks try "SF, SM" then "SF" only. The export uses UR when present (with sign flip to match FEM); when UR is not in the ODB, rotation is written as zero and `rotation_source.txt` is **none**. In `abaqus_results/job_XXX_nN/run_log.txt` check "ODB has UR: True/False" and the line "U_global.csv: U and UR read from ODB and written ..." (or "U read from ODB; rotation not in ODB, written as zero.") to confirm extraction into U_global.csv.

After validating on your installation, note here: e.g. *Tested with Abaqus 2022.*  
*Tested with Abaqus 2021 (abqpy 2021.*).*

## Expanding validation coverage

Run `run_abaqus_cae.py` for additional jobs to populate `abaqus_results` for deformation, convergence, and section-forces comparisons across point-load and distributed-load cases, e.g.:

```bash
python post_processing/validation_visualisers/abaqus/run_abaqus_cae.py --job job_0005_n16 --job job_0006_n8 --job job_0007_n8
```

For new commits, follow [docs/conventions/](docs/conventions/) and use conventional commits (e.g. `docs(validation): tested with Abaqus 2024`).

## Section forces (SF/SM) and section_forces.csv

The script **requests** both **SF** (section forces) and **SM** (section moments). When the ODB contains beam section output **SF** (and optionally **SM**), the export writes **section_forces.csv** (columns `x,N,Vy,Vz,T,My,Mz`, one row per element, x = element centroid) and **nodal_section_forces.csv** in the same format as the FEM: `# column_order=resultant`, header `N,Vy,Vz,T,My,Mz`, one row per node in node label order (nodal values are averaged from element section forces at adjacent elements). **N, Vy, Vz** come from SF; **T, My, Mz** come from SM when present in the ODB, else from SF if it has six components, else zero. In `run_log.txt` the line "Section forces: SF + SM written to ..." (or "SF (6 components)..." or "SF only...") confirms how T, My, Mz were obtained. The section-forces comparison uses `nodal_section_forces.csv` when present for direct nodal comparison.

**Our convention:** N = axial, Vy/Vz = shear, T = torsion, My/Mz = bending.

**Abaqus definitions (from Abaqus documentation):**

| Abaqus | Meaning | Note |
|--------|---------|------|
| SF1 | Axial force | → our N |
| SF2 | Transverse shear in local 2-direction | Not available for B23, B23H, **B33, B33H** |
| SF3 | Transverse shear in local 1-direction | Beams in space only; not available for **B33, B33H** |
| SM1 | Bending moment about local 1-axis | → one of My/Mz depending on local vs global |
| SM2 | Bending moment about local 2-axis | Beams in space only |
| SM3 | Twisting moment about beam axis | → our T |

For **B33/B33H** (Euler–Bernoulli 3D), SF2 and SF3 are not available; only SF1 (axial) and SM (SM1, SM2, SM3) are used. The export maps SF indices 0,1,2 → N, Vy, Vz and SM indices 0,1,2 → T, My, Mz. Abaqus component order (SM1/SM2/SM3 vs T/My/Mz) may differ by element or version; if the section-forces comparison shows sign or component mismatches, document the observed Abaqus order and any sign conventions here.

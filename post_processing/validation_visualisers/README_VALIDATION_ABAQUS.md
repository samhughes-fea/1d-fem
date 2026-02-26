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
- **deflection_tables/**: Deformation comparison and GCI–Richardson report (FEM vs Abaqus tip deflection/rotation). Scripts write to **output/**.
- **section_forces/**: Section forces comparison (SFD/BMD). Writes to **output/**.
- **output/**: Validation comparison results (overlay plots, error CSVs). The GCI–Richardson report produces `output/gci_richardson_abaqus_deflection_rotation.csv`, analogous to the verification table `gci_richardson_roark_deflection_rotation.csv` but with Abaqus (fine-mesh tip values) as reference instead of Roark.

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

3. **Run comparison** (after both FEM and Abaqus results exist):
   ```bash
   python post_processing/validation_visualisers/deflection_tables/deformation_comparison.py
   python post_processing/validation_visualisers/deflection_tables/gci_richardson_abaqus_report.py
   python post_processing/validation_visualisers/section_forces/section_forces_comparison.py
   ```
   Or run all validation visualisers:
   ```bash
   python post_processing/validation_visualisers/run_all_validation_visualisers.py
   ```

**Plots show only FEM (no Abaqus curves)?** Comparison scripts look for Abaqus CSVs in `abaqus_results/job_XXXX_nN/` (e.g. `U_global.csv`, `section_forces.csv`). If that folder is empty or missing, plots will show "FEM only". Run step 2 above for each job you want to compare so Abaqus writes results into `abaqus_results/`, then re-run the comparison scripts.

**GCI–Richardson vs Abaqus:** To build the deflection/rotation table analogous to the verification `gci_richardson_roark_deflection_rotation.csv`, run `deflection_tables/gci_richardson_abaqus_report.py`. It requires FEM results at n=32, 64, 128 and **Abaqus results at n=128** for jobs 0,1,2,5,6,7. To regenerate Abaqus scripts and run all six n128 jobs (overwrites existing ODB/CSV in `abaqus_results/`):

```bash
python post_processing/validation_visualisers/abaqus/run_abaqus_cae.py --job job_0000_n128 --job job_0001_n128 --job job_0002_n128 --job job_0005_n128 --job job_0006_n128 --job job_0007_n128
```

Then run the GCI report and LaTeX table as above. Rotation comparison requires UR in the ODB; the generated script prints which field output was requested and whether the ODB contains UR (see Table caveats below).

**Table caveats:** (1) **Tip rotation:** Abaqus reference is zero unless the ODB contains rotational DOF (UR). The Abaqus script template now requests `U` and `UR`; re-generate run scripts, re-run Abaqus, and re-extract ODB→CSV to get non-zero rotation comparison. When you run the generated script, it prints **Field output requested: U, UR, SF** or **SF only (U/UR request failed)** and **ODB has UR: True/False** after the job completes. To confirm that U and UR were correctly extracted and written, check `run_log.txt` for the line **U_global.csv: U and UR read from ODB and written (rotation sign flipped to match FEM).** (or, if UR was missing from the ODB, **U_global.csv: U read from ODB; rotation not in ODB, written as zero.**). If UR is still missing, check the Abaqus job output (e.g. `.sta` or `.msg` in the Abaqus working directory) for field-output or UR-related errors; some Abaqus or beam element setups may not output UR. When UR is not in the ODB, rotation is written as zero and `rotation_source.txt` is **none**. (2) **Distributed loads (Triangular, Parabolic):** If Abaqus deflection appears zero for jobs 6 or 7, the reference is missing or tip node ordering differs—only point-load and UDL rows then have a meaningful error vs Abaqus.

**Batch validation (optional):** For regression runs when Abaqus is available, use [run_batch_validation.py](post_processing/validation_visualisers/run_batch_validation.py): it runs Abaqus for a fixed set of jobs (e.g. `job_0000_n8`, `job_0005_n16`), runs the comparison scripts, and checks that output files exist. Example: `python post_processing/validation_visualisers/run_batch_validation.py` or `--compare-only` to skip Abaqus and only run comparisons + checks.

**Re-run all Abaqus results:** To regenerate and run Abaqus for every validation job (all `job_XXXX_nN` under `jobs/`), see [PLAN_RERUN_ABAQUS_RESULTS.md](post_processing/validation_visualisers/PLAN_RERUN_ABAQUS_RESULTS.md). Use: `python post_processing/validation_visualisers/run_all_abaqus_jobs.py` (optional: `--dry-run`, `--jobs job_0000_n8 ...`, `--no-regenerate`).

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
4. After FEM results exist for the same job, run `python post_processing/validation_visualisers/run_all_validation_visualisers.py` and check overlay plots and error CSVs in `output/`.
5. If you use a different Abaqus version, add a line under "Abaqus version notes" (e.g. *Tested with Abaqus 2024.*).

## Reviewing Abaqus result directories and performance

To audit result directories and solver performance (completion status, run time, .msg warnings/errors), run:

```bash
python post_processing/validation_visualisers/abaqus/review_abaqus_results.py
```

Options: `--expected` to compare against all jobs under `jobs/` and list jobs with no result dir; `--output <path>` for the CSV report; `--md <path>` for the Markdown summary; `--no-log-errors` to skip the errors/inconsistencies log. Default outputs: `validation_visualisers/output/abaqus_results_review.csv`, `validation_visualisers/output/abaqus_performance_summary.md`, and `validation_visualisers/output/abaqus_results_errors_and_inconsistencies.log` (run_log and .msg error/warning lines plus inconsistencies such as missing files, bad U_global header, or rotation_source vs ODB UR mismatch). The report covers file presence (U_global.csv, section_forces.csv, .inp, .odb, .sta, .msg, etc.), status from `.sta` (COMPLETED/ABORTED), total time from `.sta`, and error/warning counts from `.msg`. Validation performance (FEM vs Abaqus agreement) is produced by the comparison scripts (deformation_comparison, section_forces_comparison, gci_richardson_abaqus_report); see `output/` for those results.

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

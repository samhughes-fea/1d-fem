# Plan: Re-run all Abaqus results for each job

## Scope

Re-run Abaqus for every validation job that has a job directory under `jobs/` (pattern `job_XXXX_nN`). Each run will:

1. Regenerate the Abaqus CAE script from the job input (grid, elements, loads, etc.).
2. Run Abaqus (CAE noGUI) to build the model, run the analysis, and export results.
3. Write into `abaqus_results/job_XXXX_nN/`:
   - **U_global.csv** (then **section_forces.csv** if ODB has SF/SM)
   - **rotation_source.txt**
   - **run_log.txt**
   - Copied run artifacts: **.inp**, **.odb**, **.sta**, **.msg**

Existing result directories will be overwritten for the jobs you run.

## Job set

All validation job directories under `jobs/` that match `job_XXXX_nN` (XXXX = base job, nN = mesh size). As of this plan, that includes:

- **Point load (Euler–Bernoulli):** job_0000, job_0001, job_0002 — variants n4, n8, n16, n32, n64, n128.
- **Timoshenko / other:** job_0003_n10, job_0004_n10.
- **Distributed load (Euler–Bernoulli):** job_0005, job_0006, job_0007 — variants n4, n8, n16, n32, n64, n128.

To see the current list on your machine:

```powershell
Get-ChildItem -Path "jobs" -Directory | Where-Object { $_.Name -match '^job_\d{4}_n\d+$' } | ForEach-Object { $_.Name } | Sort-Object
```

## Prerequisites

- **Python** (project environment), **abqpy** (`pip install abqpy==2021.*` or match your Abaqus version).
- **Abaqus** with CAE; launcher set via **ABAQUS_CAE_ROOT** or **ABAQUS_BAT_PATH** (see `README_VALIDATION_ABAQUS.md`).
- Run from **project root**: `c:\Users\s1834431\Programs\fem_model`.

## Steps

### 1. Regenerate scripts and run Abaqus for all jobs

Use the helper script that discovers all `job_XXXX_nN` dirs and runs Abaqus for each:

```bash
python post_processing/validation_visualisers/run_all_abaqus_jobs.py
```

Options:

- `--dry-run` — print job list and exit (no Abaqus run).
- `--no-regenerate` — use existing generated scripts; only run Abaqus.
- `--jobs job_0000_n8 job_0005_n16 ...` — run only these jobs (default: all discovered).

Example (subset):

```bash
python post_processing/validation_visualisers/run_all_abaqus_jobs.py --jobs job_0000_n8 job_0001_n8 job_0002_n8 job_0005_n16
```

### 2. (Optional) Run results review and performance summary

After the Abaqus run, you can audit result directories and solver performance (file presence, completion status, run time, .msg counts):

```bash
python post_processing/validation_visualisers/abaqus/review_abaqus_results.py --expected
```

This writes `output/abaqus_results_review.csv` and `output/abaqus_performance_summary.md`. Use `--expected` to list jobs that exist under `jobs/` but have no `abaqus_results/` dir.

### 3. (Optional) Run comparisons and checks

After Abaqus results are written, run validation visualisers and check outputs:

```bash
python post_processing/validation_visualisers/run_batch_validation.py --compare-only
```

Or pass the same job list if you only ran a subset:

```bash
python post_processing/validation_visualisers/run_batch_validation.py --compare-only --jobs job_0000_n8 job_0005_n16
```

Note: `run_batch_validation.py` by default uses a small fixed set of jobs for its output checks; use `--jobs` to match the jobs you ran.

### 4. Verify result directories

For each run job, confirm under `post_processing/validation_visualisers/abaqus_results/job_XXXX_nN/`:

- `U_global.csv`
- `section_forces.csv` (when ODB has section forces)
- `rotation_source.txt`
- `run_log.txt`
- Optional: `job_XXXX_nN_abaqus.inp`, `.odb`, `.sta`, `.msg`

## Manual alternative (without run_all_abaqus_jobs.py)

If you prefer to call the runner directly with a fixed list:

```bash
python post_processing/validation_visualisers/abaqus/run_abaqus_cae.py --job job_0000_n4 --job job_0000_n8 --job job_0000_n16 --job job_0000_n32 --job job_0000_n64 --job job_0000_n128 --job job_0001_n4 ... 
```

(Add all 44 job names; use the discovery command above to generate the list.)

## Time and resources

- Each job runs Abaqus once (model build + analysis + ODB export). Larger meshes (e.g. n128) take longer.
- Running all jobs sequentially can take a long time; run in a batch or overnight if needed.
- Abaqus license must be available for the full duration.

## After re-run

- Deformation, GCI–Richardson, and section-force comparison scripts read from `abaqus_results/`; re-run them to refresh plots and tables.
- GCI–Richardson vs Abaqus requires Abaqus results at **n=100** for jobs 0, 1, 2, 5, 6, 7 (see `README_VALIDATION_ABAQUS.md`).

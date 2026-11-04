# Running jobs

Running a job directory end to end.

### `workflow_orchestrator.run_job`

[workflow_orchestrator/run_job.py](../../workflow_orchestrator/run_job.py)

- **`configure_logging(log_file_path)`**
  Configures logging for the application.
- **`configure_child_logging(job_results_dir)`**
  Configure logging so orchestrator + simulation_runner share ``logs/process_job.log``.
- **`get_machine_specs()`**
  Returns extended system specifications as a formatted string.
- **`track_usage()`**
  Returns current memory, disk, and CPU usage.
- **`setup_job_results_directory(case_name: str) -> str`**
  Creates the main job results directory with standard subdirectories. Returns the absolute path to the created directory.
- **`process_job(job_dir, job_results_dir, job_times, job_start_end_times, force_serial=False, max_processes_per_job=None)`**
  Processes a single FEM simulation job.
- **`main()`**
  Runs FEM simulation jobs (serially or in batches of N with a process pool).

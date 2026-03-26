# run_job.py
#
# Runs jobs in batches of 4 by default (--parallel-jobs 4), with 2 cores per job
# (--cores-per-job 2). The next batch starts only after all jobs in the current
# batch have completed and written results. Use --serial to run one job at a time
# and disable in-job parallelism. When running many jobs (e.g. 6×12), the machine
# may crash or freeze due to memory (parallel workers per job), disk, or CPU/thermal;
# use --serial or reduce --parallel-jobs / --cores-per-job if needed.

import os
import sys
import glob
import logging
import time
import concurrent.futures
import numpy as np
import psutil
import shutil
import platform
import datetime
import uuid
from tabulate import tabulate
#import cpuinfo
import subprocess
from typing import List

# Adjust Python Path to include project root
script_dir = os.path.dirname(os.path.abspath(__file__))
fem_model_root = os.path.abspath(os.path.join(script_dir, '..'))
if fem_model_root not in sys.path:
    sys.path.insert(0, fem_model_root)

# Parsing class imports

from pre_processing.parsing.element_parser import ElementParser
from pre_processing.parsing.grid_parser import GridParser
from pre_processing.parsing.material_parser import MaterialParser
from pre_processing.parsing.section_parser import SectionParser

from pre_processing.parsing.simulation_settings_parser import parse_simulation_settings
from pre_processing.parsing.point_load_parser import parse_point_load
from pre_processing.parsing.distributed_load_parser import parse_distributed_load
from pre_processing.parsing.prescribed_displacement_parser import parse_prescribed_displacement

#from pre_processing.element_library.element_factory import ElementFactory # lazy import in process job

# Configure logging for the main process
def configure_logging(log_file_path):
    """Configures logging for the application."""
    logger = logging.getLogger("run_job")
    logger.setLevel(logging.DEBUG)

    for handler in logger.handlers[:]:
        logger.removeHandler(handler)

    file_handler = logging.FileHandler(log_file_path, mode="w", encoding="utf-8")
    stream_handler = logging.StreamHandler(sys.stdout)

    formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
    file_handler.setFormatter(formatter)
    stream_handler.setFormatter(formatter)

    logger.addHandler(file_handler)
    logger.addHandler(stream_handler)

    logger.debug(f"🔧 Logging configured — output will be written to: {log_file_path}")
    return logger

def configure_child_logging(job_results_dir):
    """Configures logging for a child process."""
    log_file_path = os.path.join(job_results_dir, "logs", "process_job.log")
    logger = logging.getLogger(f"child_{os.getpid()}")
    logger.setLevel(logging.DEBUG)

    for handler in logger.handlers[:]:
        logger.removeHandler(handler)

    file_handler = logging.FileHandler(log_file_path, mode="w", encoding="utf-8")
    formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)

    return logger

def get_machine_specs():
    """Returns extended system specifications as a formatted string."""
    #cpu_info = cpuinfo.get_cpu_info()
    #cpu_name = cpu_info.get('brand_raw', platform.processor())
    logical_cores = psutil.cpu_count(logical=True)
    physical_cores = psutil.cpu_count(logical=False)
    total_ram = psutil.virtual_memory().total / (1024 ** 3)
    disk_total, _, disk_free = shutil.disk_usage("/")

    specs = (
        f"Machine Specifications:\n"
        f"   - OS: {platform.system()} {platform.release()} ({platform.version()})\n"
        #f"   - CPU: {cpu_name}\n"
        f"       • Logical cores: {logical_cores}\n"
        f"       • Physical cores: {physical_cores}\n"
        f"   - RAM: {total_ram:.2f} GB\n"
        f"   - Disk: {disk_total / (1024 ** 3):.2f} GB total, {disk_free / (1024 ** 3):.2f} GB free\n"
        f"   - Python Version: {platform.python_version()} ({sys.executable})\n"
    )

    try:
        result = subprocess.check_output("wmic path win32_VideoController get name", shell=True)
        gpus = result.decode().split('\n')[1:]
        gpus = [g.strip() for g in gpus if g.strip()]
        specs += f"   - GPU(s): {', '.join(gpus)}\n"
    except Exception:
        specs += "   - GPU(s): Unable to detect (requires Windows & WMIC)\n"

    return specs

def track_usage():
    """Returns current memory, disk, and CPU usage."""
    process = psutil.Process(os.getpid())
    return {
        "Memory (MB)": process.memory_info().rss / (1024 * 1024),
        "Disk (GB)": psutil.disk_usage('/').used / (1024 ** 3),
        "CPU (%)": process.cpu_percent(interval=0.1)
    }

def setup_job_results_directory(case_name: str) -> str:
    """
    Creates the main job results directory with standard subdirectories.
    Returns the absolute path to the created directory.

    Raises:
        ValueError: If case_name is invalid
        OSError: If directories cannot be created
    """
    if not isinstance(case_name, str) or not case_name.strip():
        raise ValueError("case_name must be a non-empty string")

    timestamp = datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S-%f')
    pid = os.getpid()
    uid = uuid.uuid4().hex[:8]

    # Protect against injection or filesystem exploits
    sanitized_case = case_name.replace(os.sep, "_").replace(" ", "_")

    results_base = os.path.join(fem_model_root, "post_processing", "results")
    job_results_dir = os.path.join(results_base, f"{sanitized_case}_{timestamp}_pid{pid}_{uid}")

    try:
        os.makedirs(job_results_dir, exist_ok=False)
    except FileExistsError:
        raise RuntimeError(f"Job directory already exists: {job_results_dir}")
    except OSError as e:
        raise RuntimeError(f"Failed to create main job directory: {e}") from e

    subdirs: List[str] = [
        "element_stiffness_matrices",
        "element_force_vectors",
        "primary_results",
        "secondary_results",
        "logs",
        "maps",
        "diagnostics"
    ]

    for subdir in subdirs:
        full_path = os.path.join(job_results_dir, subdir)
        try:
            os.makedirs(full_path, exist_ok=False)
        except FileExistsError:
            raise RuntimeError(f"Subdirectory already exists unexpectedly: {full_path}")
        except OSError as e:
            raise RuntimeError(f"Failed to create subdirectory: {full_path} -> {e}") from e

    # Optional: Confirm permissions and structure
    for subdir in subdirs:
        full_path = os.path.join(job_results_dir, subdir)
        if not os.access(full_path, os.W_OK):
            raise PermissionError(f"Subdirectory not writable: {full_path}")

    return job_results_dir


def process_job(job_dir, job_results_dir, job_times, job_start_end_times, force_serial=False, max_processes_per_job=None):
    """Processes a single FEM simulation job.

    When force_serial is True, disables per-job parallelism (element instantiation
    and stiffness/force computation) to reduce memory use and avoid machine
    instability when running many jobs (e.g. 6×12). Use --serial when the machine
    crashes or freezes on large batch runs.
    When max_processes_per_job is set (e.g. 2), caps in-job parallelism to that
    many processes per job when running multiple jobs in parallel.
    """

    case_name = os.path.basename(job_dir)
    logger = configure_child_logging(job_results_dir)
    logger.info(f"🟢 Starting job: {case_name}")

    start_time = time.time()
    usage_start = track_usage()
    performance_log_path = os.path.join(job_results_dir, "logs", "job_performance.log")
    performance_data = [["Step", "Time (s)", "Memory (MB)", "Disk (GB)", "CPU (%)"]]

    try:
        # --- PARSING INPUT FILES ---
        step_start = time.time()


        element_dictionary = ElementParser(
            filepath=os.path.join(job_dir, "element.txt"),
            job_results_dir=job_results_dir
        ).parse()["element_dictionary"]

        grid_dictionary = GridParser(
            filepath=os.path.join(job_dir, "grid.txt"),
            job_results_dir=job_results_dir
        ).parse()["grid_dictionary"]

        material_dictionary = MaterialParser(
            filepath=os.path.join(job_dir, "material.txt"),
            job_results_dir=job_results_dir
        ).parse()["material_dictionary"]

        section_dictionary = SectionParser(
            filepath=os.path.join(job_dir, "section.txt"),
            job_results_dir=job_results_dir
        ).parse()["section_dictionary"]

        simulation_settings = parse_simulation_settings(os.path.join(job_dir, "simulation_settings.txt"))
        if force_serial:
            if "parallel" not in simulation_settings:
                simulation_settings["parallel"] = {}
            simulation_settings["parallel"]["enable_parallel_instantiation"] = False
            simulation_settings["parallel"]["enable_parallel_computation"] = False
            logger.info("Serial mode: per-job parallelism disabled (single process).")
        point_load_array = np.array([])
        distributed_load_array = np.array([])
        prescribed_displacement_dict = None

        point_load_path = os.path.join(job_dir, "point_load.txt")
        if os.path.exists(point_load_path):
            point_load_array = parse_point_load(point_load_path)

        distributed_load_path = os.path.join(job_dir, "distributed_load.txt")
        if os.path.exists(distributed_load_path):
            distributed_load_array = parse_distributed_load(distributed_load_path)

        prescribed_displacement_path = os.path.join(job_dir, "prescribed_displacement.txt")
        if os.path.exists(prescribed_displacement_path):
            prescribed_displacement_dict = parse_prescribed_displacement(prescribed_displacement_path)
            # Convert to format expected by ModifyGlobalSystem
            if len(prescribed_displacement_dict['global_dof']) > 0:
                prescribed_displacement_dict = {
                    'global_dof': prescribed_displacement_dict['global_dof'],
                    'value': prescribed_displacement_dict['value']
                }
                logger.info(f"Loaded {len(prescribed_displacement_dict['global_dof'])} prescribed displacement conditions")
            else:
                prescribed_displacement_dict = None

        parsing_time = time.time() - step_start
        performance_data.append(["Parsing", parsing_time, *track_usage().values()])

        # --- ELEMENT INSTANTIATION -----------------------------------------------
        step_start = time.time()

        element_ids = element_dictionary["ids"]          # (Ne,) NumPy array

        from pre_processing.element_library.element_factory import ElementFactory
        factory = ElementFactory(job_results_dir=job_results_dir)

        # Get parallel configuration for element instantiation
        parallel_config = simulation_settings.get("parallel", {})
        enable_parallel_instantiation = parallel_config.get("enable_parallel_instantiation", True)
        num_processes_instantiation = parallel_config.get("num_processes", "auto")
        
        # Handle "auto" for num_processes
        if num_processes_instantiation == "auto":
            num_processes_instantiation = os.cpu_count() or 1
        elif isinstance(num_processes_instantiation, str):
            try:
                num_processes_instantiation = int(num_processes_instantiation)
            except ValueError:
                logger.warning(f"Invalid num_processes value '{num_processes_instantiation}', using 'auto'")
                num_processes_instantiation = os.cpu_count() or 1
        if max_processes_per_job is not None:
            num_processes_instantiation = min(num_processes_instantiation, max_processes_per_job)

        all_elements = factory.create_elements_batch(
            element_ids            = element_ids,
            element_dictionary     = element_dictionary,
            grid_dictionary        = grid_dictionary,
            material_dictionary    = material_dictionary,
            section_dictionary     = section_dictionary,
            point_load_array       = point_load_array,
            distributed_load_array = distributed_load_array,
            enable_parallel        = enable_parallel_instantiation,
            num_processes          = num_processes_instantiation,
        )

        # -------------------------------------------------------------------------
        if any(elem is None for elem in all_elements):
            logger.error(f"❌ Error: Some elements failed to instantiate in {case_name}.")
            raise ValueError(f"❌ Invalid elements detected in {case_name}.")

        required_dirs = [
            os.path.join(job_results_dir, "element_stiffness_matrices"),
            os.path.join(job_results_dir, "element_force_vectors"),
        ]
        for d in required_dirs:
            if not os.path.exists(d):
                logger.error(f"Missing required directory: {d}")
                raise FileNotFoundError(f"Directory not created: {d}")

        element_creation_time = time.time() - step_start
        performance_data.append(
            ["Element Instantiation", element_creation_time, *track_usage().values()]
        )

        # --- ELEMENT COMPUTATIONS ---
        # Get parallel configuration
        parallel_config = simulation_settings.get("parallel", {})
        enable_parallel_computation = parallel_config.get("enable_parallel_computation", True)
        num_processes = parallel_config.get("num_processes", "auto")
        
        # Handle "auto" for num_processes
        if num_processes == "auto":
            num_processes = os.cpu_count() or 1
        elif isinstance(num_processes, str):
            try:
                num_processes = int(num_processes)
            except ValueError:
                logger.warning(f"Invalid num_processes value '{num_processes}', using 'auto'")
                num_processes = os.cpu_count() or 1
        if max_processes_per_job is not None:
            num_processes = min(num_processes, max_processes_per_job)

        # Compute element stiffness matrices
        step_start = time.time()
        if enable_parallel_computation:
            try:
                from pre_processing.element_library.parallel_compute import compute_element_stiffness_parallel
                element_objects = compute_element_stiffness_parallel(
                    all_elements,
                    num_processes=num_processes
                )
            except Exception as e:
                logger.warning(f"Parallel stiffness computation failed: {e}, falling back to sequential")
                vectorized_stiffness = np.vectorize(lambda elem: elem.element_stiffness_matrix() if elem else None, otypes=[object])
                element_objects = vectorized_stiffness(all_elements)
        else:
            vectorized_stiffness = np.vectorize(lambda elem: elem.element_stiffness_matrix() if elem else None, otypes=[object])
            element_objects = vectorized_stiffness(all_elements)
        stiffness_time = time.time() - step_start
        performance_data.append(["Element Stiffness Computation", stiffness_time, *track_usage().values()])

        # Compute element force vectors
        step_start = time.time()
        if enable_parallel_computation:
            try:
                from pre_processing.element_library.parallel_compute import compute_element_force_parallel
                force_objects = compute_element_force_parallel(
                    all_elements,
                    num_processes=num_processes
                )
            except Exception as e:
                logger.warning(f"Parallel force computation failed: {e}, falling back to sequential")
                vectorized_force = np.vectorize(lambda elem: elem.element_force_vector() if elem else None, otypes=[object])
                force_objects = vectorized_force(all_elements)
        else:
            vectorized_force = np.vectorize(lambda elem: elem.element_force_vector() if elem else None, otypes=[object])
            force_objects = vectorized_force(all_elements)
        force_time = time.time() - step_start
        performance_data.append(["Element Force Computation", force_time, *track_usage().values()])

        # Validate no None in element/force objects (e.g. from failed parallel fallback or disk errors)
        elem_none = sum(1 for o in element_objects if o is None)
        force_none = sum(1 for o in force_objects if o is None)
        if elem_none or force_none:
            raise RuntimeError(
                f"Cannot run simulation: {elem_none} element object(s) and {force_none} force object(s) are None. "
                "Check disk space and logs (e.g. [Errno 28] No space left on device)."
            )

        # Compute element mass matrices when needed for modal/dynamic
        solver_type = simulation_settings.get("type", "").lower()
        element_mass_matrices = None
        if solver_type in ("modal", "dynamic"):
            step_start = time.time()
            if enable_parallel_computation:
                try:
                    from pre_processing.element_library.parallel_compute import compute_element_mass_parallel
                    mass_objects = compute_element_mass_parallel(
                        all_elements,
                        num_processes=num_processes
                    )
                except Exception as e:
                    logger.warning(f"Parallel mass computation failed: {e}, falling back to sequential")
                    mass_objects = np.array(
                        [elem.element_mass_matrix() if elem else None for elem in all_elements],
                        dtype=object
                    )
            else:
                mass_objects = np.array(
                    [elem.element_mass_matrix() if elem else None for elem in all_elements],
                    dtype=object
                )
            mass_none = sum(1 for o in mass_objects if o is None)
            if mass_none:
                raise RuntimeError(
                    f"Cannot run {solver_type} simulation: {mass_none} element(s) do not implement element_mass_matrix(). "
                    "Only element types with mass implemented (e.g. Bar-3D) are supported for modal/dynamic."
                )
            element_mass_matrices = np.array([obj.M_e for obj in mass_objects], dtype=object)
            performance_data.append(["Element Mass Computation", time.time() - step_start, *track_usage().values()])

        # --- SIMULATION EXECUTION ---
        step_start = time.time()

        if solver_type == "static":
            from simulation_runner.static.static_simulation import StaticSimulationRunner
            runner = StaticSimulationRunner(
                elements                   = all_elements,
                grid_dictionary            = grid_dictionary,
                element_dictionary         = element_dictionary,
                material_dictionary        = material_dictionary,
                section_dictionary         = section_dictionary,
                point_load_array           = point_load_array,
                distributed_load_array     = distributed_load_array,
                element_objects            = element_objects,      # NEW: ElementObject[]
                force_objects              = force_objects,        # NEW: ForceObject[]
                job_name                   = case_name,
                job_results_dir            = job_results_dir,
                simulation_settings        = simulation_settings   # NEW: Pass simulation settings
            )
            # Set prescribed displacements if provided
            if prescribed_displacement_dict is not None:
                runner.prescribed_displacements = prescribed_displacement_dict

        elif solver_type == "static_nonlinear":
            from simulation_runner.static.nonlinear_static_simulation import NonlinearStaticSimulationRunner
            runner = NonlinearStaticSimulationRunner(
                elements=all_elements,
                grid_dictionary=grid_dictionary,
                element_dictionary=element_dictionary,
                material_dictionary=material_dictionary,
                section_dictionary=section_dictionary,
                point_load_array=point_load_array,
                distributed_load_array=distributed_load_array,
                element_objects=element_objects,
                force_objects=force_objects,
                job_name=case_name,
                job_results_dir=job_results_dir,
                simulation_settings=simulation_settings,
            )
            if prescribed_displacement_dict is not None:
                runner.prescribed_displacements = prescribed_displacement_dict

        elif solver_type == "modal":
            from simulation_runner.modal.modal_simulation import ModalSimulationRunner
            # Extract stiffness and mass from element/mass objects
            element_stiffness_matrices = np.array([obj.K_e for obj in element_objects], dtype=object)
            modal_settings = {
                "elements": all_elements,
                "mesh_dictionary": {
                    "node_ids": grid_dictionary.get("ids", []),
                    "coordinates": grid_dictionary.get("coordinates", [])
                },
                "element_stiffness_matrices": element_stiffness_matrices,
                "element_mass_matrices": element_mass_matrices,
                "element_dictionary": element_dictionary,
                "grid_dictionary": grid_dictionary,
                "material_dictionary": material_dictionary,
                "section_dictionary": section_dictionary,
                "point_load_array": point_load_array,
                "distributed_load_array": distributed_load_array,
                "job_results_dir": job_results_dir,
                "simulation_settings": simulation_settings,
            }
            runner = ModalSimulationRunner(
                settings=modal_settings,
                job_name=case_name
            )

        elif solver_type == "dynamic":
            from simulation_runner.dynamic.dynamic_simulation import DynamicSimulationRunner
            element_stiffness_matrices_dyn = np.array([obj.K_e for obj in element_objects], dtype=object)
            dynamic_settings = {
                "elements": all_elements,
                "mesh_dictionary": {
                    "node_ids": grid_dictionary.get("ids", []),
                    "coordinates": grid_dictionary.get("coordinates", []),
                },
                "element_stiffness_matrices": element_stiffness_matrices_dyn,
                "element_mass_matrices": element_mass_matrices,
                "job_results_dir": job_results_dir,
                "simulation_settings": simulation_settings,
            }
            runner = DynamicSimulationRunner(settings=dynamic_settings, job_name=case_name)

        else:
            logger.error(f"❌ Unknown simulation type: '{solver_type}'")
            raise ValueError(f"Unknown simulation type: '{solver_type}'")

        runner.run()
        simulation_time = time.time() - step_start
        performance_data.append(["Full Simulation", simulation_time, *track_usage().values()])

        # --- Final Tracking ---
        end_time = time.time()
        usage_end = track_usage()
        job_times[case_name] = {"total_time": end_time - start_time}
        job_start_end_times[case_name] = (start_time, end_time)

        parallel_jobs = [
            job for job, (s, e) in job_start_end_times.items()
            if (s < end_time and e > start_time) and job != case_name
        ]

        with open(performance_log_path, "w") as f:
            f.write(get_machine_specs() + "\n")
            f.write(f"Job: {case_name}\n")
            f.write(f"Timestamp (job start): {datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S-%f')}\n")
            f.write(f"Total Time: {end_time - start_time:.2f} sec\n\n")
            f.write(tabulate(performance_data, headers="firstrow", tablefmt="grid") + "\n\n")
            f.write(f"Parallel Jobs: {', '.join(parallel_jobs) if parallel_jobs else 'None'}\n")
            f.write(f"Start Memory: {usage_start['Memory (MB)']:.2f} MB | End Memory: {usage_end['Memory (MB)']:.2f} MB\n")
            f.write(f"Start Disk Usage: {usage_start['Disk (GB)']:.2f} GB | End Disk Usage: {usage_end['Disk (GB)']:.2f} GB\n")
            f.write(f"Start CPU: {usage_start['CPU (%)']:.2f}% | End CPU: {usage_end['CPU (%)']:.2f}%\n")

    except Exception as e:
        logger.error(f"❌ Error in job {case_name}: {e}", exc_info=True)
        traceback_path = os.path.join(job_results_dir, "logs", "traceback.log")
        try:
            with open(traceback_path, "w") as f:
                import traceback
                traceback.print_exc(file=f)
        except Exception as trace_err:
            logger.error(f"⚠️ Failed to write traceback file: {trace_err}", exc_info=True)


def _run_one_job(job_dir, job_results_dir, force_serial, max_processes_per_job):
    """Worker for process pool: run one job and return timing data for main to merge."""
    job_times = {}
    job_start_end_times = {}
    process_job(
        job_dir,
        job_results_dir,
        job_times,
        job_start_end_times,
        force_serial=force_serial,
        max_processes_per_job=max_processes_per_job,
    )
    case_name = os.path.basename(job_dir)
    return (case_name, job_times.get(case_name), job_start_end_times.get(case_name))


def main():
    """Runs FEM simulation jobs (serially or in batches of N with a process pool)."""
    import argparse
    parser = argparse.ArgumentParser(
        description="Run FEM simulation jobs. Use --serial to disable in-job parallelism (recommended for large batches to avoid memory/CPU overload). "
        "Otherwise runs in batches of --parallel-jobs with --cores-per-job per job."
    )
    parser.add_argument(
        "--serial",
        action="store_true",
        help="Disable per-job parallelism (element instantiation and stiffness/force computation). Use when running many jobs (e.g. 6×12) to avoid machine crashes from memory or CPU load.",
    )
    parser.add_argument(
        "--parallel-jobs",
        type=int,
        default=4,
        metavar="N",
        help="Number of jobs to run in parallel per batch (default: 4). Ignored if --serial.",
    )
    parser.add_argument(
        "--cores-per-job",
        type=int,
        default=2,
        metavar="C",
        help="Max processes per job for in-job parallelism when using parallel batches (default: 2). Ignored if --serial.",
    )
    args = parser.parse_args()

    force_serial = args.serial
    parallel_jobs = 1 if force_serial else args.parallel_jobs
    cores_per_job = args.cores_per_job

    log_file_path = os.path.join(script_dir, "run_job.log")
    logger = configure_logging(log_file_path)
    logger.info("🚀 Starting FEM Simulation Workflow")

    jobs_dir = os.path.join(fem_model_root, 'jobs')
    job_dirs = [d for d in glob.glob(os.path.join(jobs_dir, 'job_*')) if os.path.isdir(d)]

    if not job_dirs:
        logger.warning("⚠️ No job directories found.")
        return

    if force_serial:
        logger.info("🟢 Serial mode: in-job parallelism disabled (single process per job).")
    job_times = {}
    job_start_end_times = {}

    if parallel_jobs == 1:
        logger.info(f"🟢 Running {len(job_dirs)} jobs in serial.")
        for job_dir in job_dirs:
            case_name = os.path.basename(job_dir)
            try:
                job_results_dir = setup_job_results_directory(case_name)
            except Exception as e:
                logger.error(f"❌ Failed to create job results dir for {case_name}: {e}", exc_info=True)
                continue
            process_job(job_dir, job_results_dir, job_times, job_start_end_times, force_serial=force_serial)
    else:
        logger.info(f"🟢 Running {len(job_dirs)} jobs in batches of {parallel_jobs} (parallel_jobs={parallel_jobs}, cores_per_job={cores_per_job}).")
        batches = [job_dirs[i:i + parallel_jobs] for i in range(0, len(job_dirs), parallel_jobs)]
        with concurrent.futures.ProcessPoolExecutor(max_workers=parallel_jobs) as executor:
            for batch in batches:
                batch_results_dirs = []
                for job_dir in batch:
                    case_name = os.path.basename(job_dir)
                    try:
                        job_results_dir = setup_job_results_directory(case_name)
                        batch_results_dirs.append((job_dir, job_results_dir))
                    except Exception as e:
                        logger.error(f"❌ Failed to create job results dir for {case_name}: {e}", exc_info=True)
                if not batch_results_dirs:
                    continue
                futures = [
                    executor.submit(_run_one_job, job_dir, job_results_dir, force_serial, cores_per_job)
                    for job_dir, job_results_dir in batch_results_dirs
                ]
                for fut in concurrent.futures.as_completed(futures):
                    try:
                        case_name, times, start_end = fut.result()
                        if case_name and times is not None:
                            job_times[case_name] = times
                        if case_name and start_end is not None:
                            job_start_end_times[case_name] = start_end
                    except Exception as e:
                        logger.error(f"❌ Job failed: {e}", exc_info=True)

if __name__ == "__main__":
    main()
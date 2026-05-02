"""One-off: write 4×25-element cantilever jobs and run them (linear/nl × EB/Timoshenko)."""
from __future__ import annotations

import os
import sys

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

NELEM = 25
NNODES = NELEM + 1
L_BEAM = 2.0
DX = L_BEAM / NELEM

HEADER_ELE = (
    "[Element]\n"
    "[element_id]  [node1]  [node2]  [element_type]                [axial_order]     "
    "[bending_y_order] [bending_z_order] [shear_y_order]   [shear_z_order]   [torsion_order]   [load_order]\n"
)

MATERIAL_ROW = "0             2.1e+11     8.1000e+10   0.3    7850\n"
SECTION_ROW = (
    "0             0.00131       0.0           3.23400e-07       2.08769e-06       2.60673e-08\n"
)


def _write_common(job_dir: str) -> None:
    lines_grid = ["[Grid]", "[node_id]   [x]         [y]       [z]"]
    for i in range(NNODES):
        x = i * DX
        lines_grid.append(f"{i}           {x:.6f}    0.0       0.0")
    with open(os.path.join(job_dir, "grid.txt"), "w", encoding="utf-8") as f:
        f.write("\n".join(lines_grid) + "\n")

    lines_mat = ["[Material]", "[element_id]  [E]         [G]          [nu]   [rho]"]
    lines_sec = ["[Section]", "[element_id]  [A]           [I_x]         [I_y]             [I_z]             [J_t]"]
    for e in range(NELEM):
        lines_mat.append(MATERIAL_ROW.replace("0 ", f"{e} ", 1))
        lines_sec.append(SECTION_ROW.replace("0 ", f"{e} ", 1))
    with open(os.path.join(job_dir, "material.txt"), "w", encoding="utf-8") as f:
        f.write("\n".join(lines_mat) + "\n")
    with open(os.path.join(job_dir, "section.txt"), "w", encoding="utf-8") as f:
        f.write("\n".join(lines_sec) + "\n")

    pre = ["# test", "[Precurvature]", "[element_id] [k_x0] [k_y0] [k_z0]"]
    for e in range(NELEM):
        pre.append(f"{e} 0.0 0.0 0.0")
    with open(os.path.join(job_dir, "precurvature.txt"), "w", encoding="utf-8") as f:
        f.write("\n".join(pre) + "\n")

    pd = """[Prescribed Displacement]
[id]     [node_id]  [dof]   [value]     [type]          [comment]
0        0          UX      0.0         displacement     # Fixed
1        0          UY      0.0         displacement
2        0          UZ      0.0         displacement
3        0          RX      0.0         displacement
4        0          RY      0.0         displacement
5        0          RZ      0.0         displacement
"""
    with open(os.path.join(job_dir, "prescribed_displacement.txt"), "w", encoding="utf-8") as f:
        f.write(pd)

    pl = f"""# End load
[Point load]
         [x]          [y]          [z]        [F_x]        [F_y]        [F_z]        [M_x]        [M_y]        [M_z]
    {L_BEAM:.6f}     0.000000     0.000000     0.000000  -500.000000     0.000000     0.000000     0.000000     0.000000
"""
    with open(os.path.join(job_dir, "point_load.txt"), "w", encoding="utf-8") as f:
        f.write(pl)


def _write_element_file(job_dir: str, element_type: str, *, timoshenko: bool) -> None:
    rows = [HEADER_ELE]
    for e in range(NELEM):
        n1, n2 = e, e + 1
        if timoshenko:
            row = (
                f"{e}             {n1}        {n2}        {element_type}   3                 3                 3                 "
                f"2                 2                 3                 2                 \n"
            )
        else:
            row = (
                f"{e}             {n1}        {n2}        {element_type}   3                 3                 3                 "
                f"0                 0                 3                 2                 \n"
            )
        rows.append(row)
    with open(os.path.join(job_dir, "element.txt"), "w", encoding="utf-8") as f:
        f.write("".join(rows))


def _write_sim_settings(job_dir: str, *, nonlinear: bool) -> None:
    if not nonlinear:
        text = "[Simulation]\n[Type]\nStatic\n"
    else:
        text = (
            "[Simulation]\n[Type]\nstatic_nonlinear\n\n[Newton]\n"
            "tolerance = 1e-4\n"
            "relative_tolerance = 1e-6\n"
            "relative_reference = first_residual\n"
            "max_iterations = 40\n"
            "tolerance_delta_u = 1e-9\n\n"
            "[Solver]\n"
            "tolerance = 1e-10\n"
            "max_iterations = 2000\n\n"
            "[Nonlinear]\n"
            "num_increments = 1\n"
            "line_search = false\n"
        )
    with open(os.path.join(job_dir, "simulation_settings.txt"), "w", encoding="utf-8") as f:
        f.write(text)


def write_all_jobs() -> list[str]:
    jobs_dir = os.path.join(ROOT, "jobs")
    os.makedirs(jobs_dir, exist_ok=True)
    cases: list[tuple[str, str, bool, bool]] = [
        ("job_test_n25_linear_eb", "LinearEulerBernoulliBeamElement3D", False, False),
        ("job_test_n25_linear_timoshenko", "LinearTimoshenkoBeamElement3D", True, False),
        ("job_test_n25_nonlinear_eb", "NonlinearEulerBernoulliBeamElement3D", False, True),
        ("job_test_n25_nonlinear_timoshenko", "NonlinearTimoshenkoBeamElement3D", True, True),
    ]
    names: list[str] = []
    for folder, etype, is_timo, is_nl in cases:
        d = os.path.join(jobs_dir, folder)
        os.makedirs(d, exist_ok=True)
        _write_common(d)
        _write_element_file(d, etype, timoshenko=is_timo)
        _write_sim_settings(d, nonlinear=is_nl)
        names.append(folder)
    return names


def run_jobs(job_names: list[str]) -> None:
    from workflow_orchestrator.run_job import process_job, setup_job_results_directory

    job_times: dict = {}
    job_start_end: dict = {}
    for name in job_names:
        job_dir = os.path.join(ROOT, "jobs", name)
        res_dir = setup_job_results_directory(name)
        print(f"Running {name} -> {res_dir}")
        process_job(
            job_dir,
            res_dir,
            job_times,
            job_start_end,
            force_serial=True,
            max_processes_per_job=1,
        )
        print(f"  done. times: {job_times.get(name)}")


if __name__ == "__main__":
    n = write_all_jobs()
    run_jobs(n)

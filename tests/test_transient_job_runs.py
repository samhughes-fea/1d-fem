"""A transient job runs end to end through ``process_job``.

``process_job`` passed the elements as a list and the element matrices stacked into one
object array, and left out the loads; the runner failed on ``list.size`` and the failure
went only to ``logs/traceback.log``. The tests that called the runner directly missed it.
"""

from pathlib import Path

import numpy as np

from onedfem.workflow_orchestrator.run_job import process_job

JOB = Path(__file__).resolve().parents[1] / "jobs" / "job_benchmark_transient_cantilever_multidof"
SUBDIRS = ("element_stiffness_matrices", "element_force_vectors", "primary_results", "secondary_results",
           "logs", "maps", "diagnostics")


def test_transient_benchmark_job_runs(tmp_path):
    res = tmp_path / "res"
    for sub in SUBDIRS:
        (res / sub).mkdir(parents=True)
    process_job(str(JOB), str(res), {}, {}, force_serial=True, max_processes_per_job=1)
    tb = res / "logs" / "traceback.log"
    assert not (tb.exists() and tb.read_text().strip()), tb.read_text()
    out = res / "primary_results" / "dynamic_results"
    u = np.loadtxt(out / f"{JOB.name}_displacements.txt", comments="#", ndmin=2)
    t = np.loadtxt(out / f"{JOB.name}_time.txt", comments="#", ndmin=1)
    assert u.shape[0] == t.size > 1
    assert np.all(u[0] == 0.0)            # starts at rest
    assert np.abs(u[1:]).max() > 0.0      # the tip load drives the response

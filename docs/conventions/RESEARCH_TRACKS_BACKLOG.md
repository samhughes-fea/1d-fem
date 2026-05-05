# Research tracks (out of band from §4 harmonic releases)

Track these as separate issues or milestone epics; they are **not** blocking harmonic frequency-domain work.

| Track | Design / reference | Notes |
|-------|-------------------|--------|
| **Eigen / static assembly deduplication** | [`processing/eigen/assembly.py`](../../processing/eigen/assembly.py), [`processing/static/operations/assembly.py`](../../processing/static/operations/assembly.py) | Extract shared scatter/BC/validation helpers in a **dedicated refactor PR** (behavior-preserving); regress eigen/buckling/harmonic smokes and linear/nonlinear static jobs. |
| **New primary analysis domains** | [`simulation_runner/README.md`](../../simulation_runner/README.md) (section *Template for new analysis domains*) | Add `processing/<domain>/{operations,results,diagnostics}` and a thin `simulation_runner/<domain>/` runner; extend parser taxonomy and [`workflow_orchestrator/run_job.py`](../../workflow_orchestrator/run_job.py) dispatch when the analysis type is implemented. |
| **LTB** lateral-torsional buckling validation | [`MODAL_BUCKLING_LTB_VALIDATION.md`](MODAL_BUCKLING_LTB_VALIDATION.md) | Closed-form or benchmark **P_cr** / load factor once reference BCs and loads are fixed. |
| **GESDB** strain recovery milestones | [`geometrically_exact_shear_deformable_beam_formulation.md`](../element_library/geometrically_exact_shear_deformable_beam_formulation.md) | Native large-rotation beam roadmap items. |
| **Nonlinear prestress → linear buckling** | [`MODAL_BUCKLING_NONLINEAR_PRESTRESS_DESIGN.md`](MODAL_BUCKLING_NONLINEAR_PRESTRESS_DESIGN.md) | Prestressed **K** for §5 buckling vs transient nonlinear paths. |
| **Nonlinear buckling continuation benchmarks** | [`NONLINEAR_BUCKLING_CONTINUATION.md`](NONLINEAR_BUCKLING_CONTINUATION.md), [`NONLINEAR_BUCKLING_BENCHMARKS.md`](NONLINEAR_BUCKLING_BENCHMARKS.md) | Progress from smoke/regression to pinned imperfect-column and post-critical reference cases. |
| **Non-static production-readiness parity** | [`NONSTATIC_PRODUCTION_READINESS_CHECKLIST.md`](NONSTATIC_PRODUCTION_READINESS_CHECKLIST.md) | Bring eigen, transient, harmonic, and buckling runners/packages to static-grade artifact, diagnostics, and validation standards. |
| **Eigen + linear buckling acceptance benchmarks** | [`EIGEN_AND_LINEAR_BUCKLING_BENCHMARKS.md`](EIGEN_AND_LINEAR_BUCKLING_BENCHMARKS.md) | Pinned repository benchmark jobs and artifact/positivity acceptance checks before literature-calibrated closures. |

When opening a GitHub issue, link the doc above and the simulation type (**eigen**, **buckling**, **transient**, **harmonic**) affected.

## Branch WIP (snapshot)

Large batches of `jobs/job_*_*/simulation_settings.txt` updates that remain unstaged are intentionally split from taxonomy PRs: commit or revert them in a dedicated **jobs** or **chore** change when ready.

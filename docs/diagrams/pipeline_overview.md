# Pipeline overview

End-to-end FEM pipeline from job discovery through parsing, element creation, simulation execution, and results. The logic is implemented in `workflow_orchestrator/run_job.py` (`main()` and `process_job()`) and the simulation runners under `simulation_runner/`.

```mermaid
flowchart TB
    subgraph jobs [Job input]
        JobDirs["jobs/job_*/"]
    end

    subgraph orchestration [Workflow orchestrator]
        Main["run_job.main()"]
        ProcessJob["process_job()"]
    end

    subgraph parsing [Parsing]
        Parse["Parse input files\n(element, grid, material, section,\nsimulation_settings, loads, prescribed_disp)"]
    end

    subgraph elements [Element phase]
        Factory["ElementFactory.create_elements_batch()"]
        KFe["Compute K_e and F_e\n(parallel or sequential)"]
    end

    subgraph runner [Simulation type]
        Static["StaticSimulationRunner.run()"]
        Modal["ModalSimulationRunner.run()"]
    end

    subgraph results [Results]
        ResultDir["post_processing/results/\n{case}_{timestamp}_..."]
    end

    subgraph post [Post-processing]
        Graphical["graphical_visualisers"]
        Verification["verification_visualisers"]
        Tensor["tensor_visualisers"]
    end

    JobDirs --> Main
    Main --> ProcessJob
    ProcessJob --> Parse
    Parse --> Factory
    Factory --> KFe
    KFe --> Static
    KFe --> Modal
    Static --> ResultDir
    Modal --> ResultDir
```

No change needed — both runners can produce results. Continuing with the static simulation flow.
<｜tool▁calls▁begin｜><｜tool▁call▁begin｜>
TodoWrite
    ResultDir --> Graphical
    ResultDir --> Verification
    ResultDir --> Tensor
```

- **Job discovery:** `main()` finds all `job_*` directories under `jobs/`, creates a result directory per job, and runs `process_job()` (optionally in parallel).
- **Parsing:** Element, grid, material, section, simulation settings, and optional point/distributed loads and prescribed displacements are read from the job directory.
- **Element phase:** `ElementFactory` builds elements; then element stiffness matrices and force vectors are computed (parallel or sequential per `simulation_settings.parallel`).
- **Simulation:** Branch by `simulation_settings.type` — `"static"` runs `StaticSimulationRunner.run()`, `"modal"` runs `ModalSimulationRunner.run()`.
- **Results:** Written under `post_processing/results/{case}_{timestamp}_pid{pid}_{uid}/` (primary_results, secondary_results, tertiary_results, maps, logs, etc.).
- **Post-processing:** Separate scripts under `post_processing/` read these results (graphical, verification, tensor visualisers).

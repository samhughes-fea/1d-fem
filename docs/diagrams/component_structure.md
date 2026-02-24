# Component structure

Top-level modules and their main subcomponents. Dependencies flow from job input through pre_processing and processing into the simulation runners; post_processing consumes results.

```mermaid
flowchart TB
    subgraph workflow [workflow_orchestrator]
        RunJob["run_job.py\nmain(), process_job()"]
    end

    subgraph pre [pre_processing]
        Parsing["parsing/\nElementParser, GridParser, MaterialParser,\nSectionParser, simulation_settings,\npoint_load, distributed_load, prescribed_displacement"]
        ElementLib["element_library/\nElementFactory, element_1D_base,\ngauss_point_data, parallel_compute"]
        Bar["element_library/bar/"]
        Truss["element_library/truss/"]
        EulerBernoulli["element_library/euler_bernoulli/"]
        Timoshenko["element_library/timoshenko/"]
        Levinson["element_library/levinson/"]
        MeshLib["mesh_library/\ncreate_*_mesh_variants"]
        LoadLib["load_library/\nschemes (distributed_loads, etc.)"]
    end

    subgraph processing [processing]
        StaticOps["static/operations/\nPrepareLocalSystem, AssembleGlobalSystem,\nModifyGlobalSystem, CondenseModifiedSystem,\nSolveCondensedSystem, ReconstructGlobalSystem,\nDisassembleGlobalSystem"]
        StaticResults["static/results/\ncompute_primary, compute_secondary, compute_tertiary,\ncontainers, save_*_container"]
        StaticDiag["static/diagnostics/"]
        ModalProc["modal/\nassembly, boundary_conditions"]
    end

    subgraph runner [simulation_runner]
        StaticRunner["static/StaticSimulationRunner"]
        ModalRunner["modal/ModalSimulationRunner"]
    end

    subgraph post [post_processing]
        Graphical["graphical_visualisers/\nprimary, secondary, tertiary"]
        Verification["verification_visualisers/\nroarks_formulas, deflection_tables"]
        Tensor["tensor_visualisers/"]
    end

    RunJob --> Parsing
    RunJob --> ElementLib
    RunJob --> StaticRunner
    RunJob --> ModalRunner

    ElementLib --> Bar
    ElementLib --> Truss
    ElementLib --> EulerBernoulli
    ElementLib --> Timoshenko
    ElementLib --> Levinson

    StaticRunner --> StaticOps
    StaticRunner --> StaticResults
    StaticRunner --> StaticDiag
    ModalRunner --> ModalProc

    StaticResults --> post
```

## Module roles

| Module | Role |
|--------|------|
| **workflow_orchestrator** | Discovers jobs under `jobs/`, creates result dirs, invokes `process_job()` (parsing, element creation, K_e/F_e computation, runner selection and execution). |
| **pre_processing/parsing** | Reads job input files into dictionaries/arrays (element, grid, material, section, settings, loads, prescribed displacements). |
| **pre_processing/element_library** | ElementFactory builds elements; bar, truss, euler_bernoulli, timoshenko, levinson implement K_e and F_e; gauss_point_data and parallel_compute support formulation cache and parallel assembly. |
| **pre_processing/mesh_library** | Mesh generation variants (e.g. create_point_load_mesh_variants, create_distributed_mesh_variants). |
| **pre_processing/load_library** | Load schemes (distributed loads, equivalent line loads, etc.). |
| **processing/static** | Operations (prepare, assemble, modify, condense, solve, reconstruct, disassemble); results (primary, secondary, tertiary); containers; diagnostics. |
| **processing/modal** | Global assembly and boundary conditions for modal analysis. |
| **simulation_runner** | StaticSimulationRunner runs the full linear-static workflow; ModalSimulationRunner runs modal (K/M assembly, BCs, eigenvalue solve). |
| **post_processing** | Scripts that read result directories: graphical (deformation, load, stress, strain, section forces, etc.), verification (Roark, deflection convergence, GCI), tensor visualisers. |

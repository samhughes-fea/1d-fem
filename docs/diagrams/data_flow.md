# Data flow

Job input files and result output layout. Each job directory under `jobs/` (e.g. `jobs/job_0000_n8/`) is read by `process_job()`; results are written to a new directory under `post_processing/results/`.

```mermaid
flowchart LR
    subgraph inputs [Job input directory]
        ElementTxt["element.txt"]
        GridTxt["grid.txt"]
        MaterialTxt["material.txt"]
        SectionTxt["section.txt"]
        SettingsTxt["simulation_settings.txt"]
        PointLoad["point_load.txt\n(optional)"]
        DistLoad["distributed_load.txt\n(optional)"]
        PrescribedDisp["prescribed_displacement.txt\n(optional)"]
        PrecurvatureTxt["precurvature.txt\n(optional)"]
    end

    subgraph pipeline [Pipeline]
        ProcessJob["process_job()"]
    end

    subgraph outputs [Job results directory]
        PrimaryDir["primary_results/"]
        SecondaryDir["secondary_results/"]
        TertiaryDir["tertiary_results/"]
        MapsDir["maps/"]
        LogsDir["logs/"]
        DiagDir["diagnostics/"]
        KeDir["element_stiffness_matrices/"]
        FeDir["element_force_vectors/"]
    end

    subgraph artifacts [Key artifacts]
        UGlobal["U_global.csv"]
        RGlobal["R_global.csv"]
        IndexMaps["Index maps\n(assembly, condensation, etc.)"]
    end

    ElementTxt --> ProcessJob
    GridTxt --> ProcessJob
    MaterialTxt --> ProcessJob
    SectionTxt --> ProcessJob
    SettingsTxt --> ProcessJob
    PointLoad --> ProcessJob
    DistLoad --> ProcessJob
    PrescribedDisp --> ProcessJob
    PrecurvatureTxt --> ProcessJob

    ProcessJob --> PrimaryDir
    ProcessJob --> SecondaryDir
    ProcessJob --> TertiaryDir
    ProcessJob --> MapsDir
    ProcessJob --> LogsDir
    ProcessJob --> DiagDir
    ProcessJob --> KeDir
    ProcessJob --> FeDir

    PrimaryDir --> UGlobal
    PrimaryDir --> RGlobal
    MapsDir --> IndexMaps
```

## Input files (job directory)

| File | Required | Parser / usage |
|------|----------|----------------|
| element.txt | Yes | ElementParser — connectivity, element type |
| grid.txt | Yes | GridParser — node IDs, coordinates |
| material.txt | Yes | MaterialParser — E, nu, rho, etc. |
| section.txt | Yes | SectionParser — section properties |
| simulation_settings.txt | Yes | parse_simulation_settings — solver type (static/modal/dynamic), solver config, parallel options; **[Modal]** section keys include `analysis` (`vibration` \| `buckling`), `buckling_prestress` (`linear_static`; `none` is rejected), `buckling_load_factor` (scales reference loads for prestress used to build geometric stiffness) |
| point_load.txt | No | parse_point_load — concentrated loads |
| distributed_load.txt | No | parse_distributed_load — line loads |
| prescribed_displacement.txt | No | parse_prescribed_displacement — fixed/prescribed DOFs |
| precurvature.txt | No | parse_precurvature — per-element reference ``(k_x0, k_y0, k_z0)`` (1/m) for straight EB/Timoshenko; omitted or all-zero preserves prior behaviour |

**Beam warping and section tiers:** optional `[warping]` / `[curvature]` columns in `element.txt` and 6/8/10/11-column layouts for `section.txt` (including `Gamma` and shear centre) are summarized in [`docs/conventions/JOB_INPUT_BEAM_WARPING.md`](../conventions/JOB_INPUT_BEAM_WARPING.md).

## Output directory layout

Result root: `post_processing/results/{case_name}_{timestamp}_pid{pid}_{uid}/`

| Subdirectory | Contents |
|--------------|----------|
| primary_results | U_global, R_global, F_global, K_global (or summaries); CSVs per resolution (global, elemental); primary results summary |
| secondary_results | Strain, stress, energy density (Gauss and nodal); secondary summary |
| tertiary_results | Section forces, principal stress, etc.; tertiary summary |
| maps | Assembly, modification, condensation, reconstruction index maps |
| logs | process_job.log, job_performance.log, traceback on failure |
| diagnostics | Linear-static diagnostics, runtime telemetry |
| element_stiffness_matrices | Per-element K_e (from element phase) |
| element_force_vectors | Per-element F_e (from element phase) |

Post-processing scripts (graphical, verification, tensor visualisers) discover results by scanning `post_processing/results/` and parsing job result directory names.

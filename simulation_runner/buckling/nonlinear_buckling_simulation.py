from __future__ import annotations

import logging
import os

import numpy as np

from processing.common.nonlinear_buckling_continuation import (
    ArcLengthStepResult,
    ContinuationConfig,
    ImperfectionConfig,
    append_continuation_history_row,
    predictor_load_factor_for_increment,
    seed_initial_imperfection,
    write_continuation_history_header,
    write_continuation_summary,
)
from processing.common.nonlinear_equilibrium import solve_arc_length_corrector
from simulation_runner.static.nonlinear_static_simulation import NonlinearStaticSimulationRunner

logger = logging.getLogger(__name__)


class NonlinearBucklingSimulationRunner:
    """Load-control nonlinear buckling / continuation MVP."""

    def __init__(self, settings: dict, job_name: str):
        self.settings = settings
        self.job_name = job_name
        self.simulation_settings = settings.get("simulation_settings") or {}
        self.job_results_dir = settings.get("job_results_dir")
        self.primary_results_dir = os.path.join(self.job_results_dir, "primary_results")
        self.diagnostics_dir = os.path.join(self.job_results_dir, "diagnostics")
        self.nlb_results_dir = os.path.join(self.primary_results_dir, "nonlinear_buckling_results")

    def _build_nonlinear_settings(self) -> dict:
        sim = dict(self.simulation_settings)
        buckling = dict(sim.get("buckling") or {})
        nonlinear = dict(sim.get("nonlinear") or {})
        if "num_increments" not in nonlinear and "num_increments" in buckling:
            nonlinear["num_increments"] = buckling.get("num_increments")
        if "load_factors" not in nonlinear and buckling.get("load_factors") is not None:
            nonlinear["load_factors"] = buckling.get("load_factors")
        for key in ("line_search", "line_search_max_backtracks", "line_search_shrink"):
            if key not in nonlinear and key in buckling:
                nonlinear[key] = buckling.get(key)
        sim["nonlinear"] = nonlinear
        return sim

    def _continuation_method(self) -> str:
        buckling = dict(self.simulation_settings.get("buckling") or {})
        method = str(buckling.get("continuation_method", "load_control")).strip().lower()
        if method not in {"load_control", "arc_length"}:
            raise ValueError("buckling.continuation_method must be 'load_control' or 'arc_length'")
        return method

    def _arc_length_settings(self) -> tuple[float, float]:
        buckling = dict(self.simulation_settings.get("buckling") or {})
        radius = float(buckling.get("arc_length_radius", 1.0))
        alpha_scale = float(buckling.get("arc_length_alpha_scale", 1.0))
        return radius, alpha_scale

    def _imperfection_settings(self) -> tuple[str | None, int, float]:
        buckling = dict(self.simulation_settings.get("buckling") or {})
        source = buckling.get("imperfection_source")
        if source is None:
            return None, 0, 0.0
        source_s = str(source).strip().lower()
        mode_index = int(buckling.get("imperfection_mode_index", 0))
        scale = float(buckling.get("imperfection_scale", 0.0))
        return source_s, mode_index, scale

    def run(self) -> None:
        if not self.job_results_dir:
            raise ValueError("NonlinearBucklingSimulationRunner requires job_results_dir in settings")
        os.makedirs(self.diagnostics_dir, exist_ok=True)
        os.makedirs(self.nlb_results_dir, exist_ok=True)

        runner = NonlinearStaticSimulationRunner(
            elements=self.settings["elements"],
            grid_dictionary=self.settings["grid_dictionary"],
            element_dictionary=self.settings["element_dictionary"],
            material_dictionary=self.settings["material_dictionary"],
            section_dictionary=self.settings["section_dictionary"],
            point_load_array=self.settings["point_load_array"],
            distributed_load_array=self.settings["distributed_load_array"],
            element_objects=self.settings["element_objects"],
            force_objects=self.settings["force_objects"],
            job_name=self.job_name,
            job_results_dir=self.job_results_dir,
            simulation_settings=self._build_nonlinear_settings(),
        )
        runner.setup_simulation()
        runner.prepare_local_system(job_results_dir=runner.primary_results_dir)
        runner.assemble_global_system(runner.primary_results_dir)
        F_ext_full = np.asarray(runner.F_global, dtype=np.float64).ravel()

        U_global = np.zeros(runner.total_dof, dtype=np.float64)
        if runner.prescribed_displacements is not None:
            gd = np.asarray(runner.prescribed_displacements["global_dof"], dtype=np.int32)
            val = np.asarray(runner.prescribed_displacements["value"], dtype=np.float64)
            U_global[gd] = val
        source, mode_index, scale = self._imperfection_settings()
        U_global, imperfection_meta = seed_initial_imperfection(
            job_name=self.job_name,
            primary_results_dir=self.primary_results_dir,
            job_results_dir=self.job_results_dir,
            U_global=U_global,
            config=ImperfectionConfig(source=source, mode_index=mode_index, scale=scale),
        )

        load_factors = runner._compute_load_factors()
        continuation_method = self._continuation_method()
        arc_length_radius, arc_length_alpha_scale = self._arc_length_settings()
        continuation_config = ContinuationConfig(
            continuation_method=continuation_method,
            load_factors=np.asarray(load_factors, dtype=np.float64),
            arc_length_radius=arc_length_radius,
            arc_length_alpha_scale=arc_length_alpha_scale,
        )
        history_path = os.path.join(self.nlb_results_dir, "continuation_history.csv")
        completed = 0
        last_load_factor = 0.0
        write_continuation_history_header(history_path)
        tip_dof = int(max(0, runner.total_dof - runner.dof_per_node))
        for zero_i, lam in enumerate(load_factors):
            inc_i = zero_i + 1
            runner.load_increment_index = inc_i
            runner.last_load_factor = float(lam)
            predictor_load_factor = predictor_load_factor_for_increment(
                config=continuation_config,
                increment_index=zero_i,
                current_U=U_global,
                tip_dof=tip_dof,
                reference_load_vector=F_ext_full,
            )
            step_result = self._solve_continuation_step(
                runner=runner,
                continuation_config=continuation_config,
                increment_index=zero_i,
                requested_load_factor=float(lam),
                predictor_load_factor=float(predictor_load_factor),
                current_U=U_global,
                reference_load_vector=F_ext_full,
                tip_dof=tip_dof,
            )
            U_global = step_result.U_global if hasattr(step_result, "U_global") else U_global
            effective_load_factor = float(step_result.load_factor)
            converged = bool(step_result.converged)
            iterations_used = int(step_result.iterations_used)
            residual = step_result.last_norm_F_cond if hasattr(step_result, "last_norm_F_cond") else step_result.residual_norm
            append_continuation_history_row(
                path=history_path,
                increment_index=inc_i,
                load_factor=float(lam),
                converged=bool(converged),
                iterations_used=int(iterations_used),
                residual_norm=residual,
                tip_dof=tip_dof,
                tip_displacement=float(U_global[tip_dof]) if runner.total_dof else 0.0,
                continuation_method=continuation_method,
                predictor_load_factor=float(step_result.predictor_load_factor if hasattr(step_result, "predictor_load_factor") else predictor_load_factor),
            )
            last_load_factor = float(effective_load_factor)
            if converged:
                completed = inc_i
            else:
                break

        summary = {
            "job_name": self.job_name,
            "continuation_method": continuation_method,
            "increments_requested": int(len(load_factors)),
            "increments_completed": int(completed),
            "last_load_factor": float(last_load_factor),
            "history_csv": os.path.relpath(history_path, self.job_results_dir).replace("\\", "/"),
        }
        if imperfection_meta is not None:
            summary["imperfection"] = imperfection_meta
        summary_path = os.path.join(self.diagnostics_dir, "nonlinear_buckling_summary.json")
        write_continuation_summary(summary_path, summary)
        diag_log = os.path.join(self.diagnostics_dir, "nonlinear_buckling_diagnostic.log")
        with open(diag_log, "w", encoding="utf-8") as f:
            f.write(
                "Nonlinear buckling continuation MVP executed.\n"
                f"increments_completed={completed}\n"
                f"history_csv={summary['history_csv']}\n"
                "See docs/conventions/NONLINEAR_BUCKLING_CONTINUATION.md\n"
            )
        logger.info("Nonlinear buckling continuation history written to %s", history_path)

    def _solve_continuation_step(
        self,
        *,
        runner: NonlinearStaticSimulationRunner,
        continuation_config: ContinuationConfig,
        increment_index: int,
        requested_load_factor: float,
        predictor_load_factor: float,
        current_U: np.ndarray,
        reference_load_vector: np.ndarray,
        tip_dof: int,
    ):
        if continuation_config.continuation_method != "arc_length":
            F_ext_global = float(requested_load_factor) * reference_load_vector
            before = runner.newton_iterations_total
            U_global, _, converged = runner._newton_raphson_solve(
                F_ext_global,
                current_U,
                float(requested_load_factor),
            )
            return ArcLengthStepResult(
                load_factor=float(requested_load_factor),
                predictor_load_factor=float(predictor_load_factor),
                converged=bool(converged),
                iterations_used=int(runner.newton_iterations_total - before),
                residual_norm=runner.last_norm_F_cond,
                tip_displacement=float(U_global[tip_dof]) if runner.total_dof else 0.0,
            )

        def build_system_from_state(U_state: np.ndarray, load_factor_state: float):
            F_state = float(load_factor_state) * reference_load_vector
            with runner.monitor.stage("ArcLengthBuild"):
                K_T_list, F_int_list = runner._build_K_T_and_F_int(U_state)
                runner.assemble_global_system(
                    runner.primary_results_dir,
                    element_stiffness_matrices=K_T_list,
                    element_force_vectors=F_int_list,
                )
                F_int_global = np.asarray(runner.F_global, dtype=np.float64).ravel()
            R_global = F_state - F_int_global
            with runner.monitor.stage("ArcLengthModify"):
                runner.modify_global_system(
                    runner.K_global,
                    R_global,
                    runner.local_global_dof_map,
                    runner.primary_results_dir,
                    prescribed_displacements=getattr(runner, "prescribed_displacements", None),
                )
                R_mod = runner.F_mod
                runner.condense_modified_system(
                    runner.K_mod,
                    np.asarray(R_mod).ravel(),
                    runner.fixed_dofs,
                    runner.local_global_dof_map,
                    runner.primary_results_dir,
                    base_tol=runner.condensation_config.get("base_tol", 1e-12),
                )
            F_cond_arr = np.asarray(runner.F_cond, dtype=np.float64).ravel()
            return float(np.linalg.norm(R_global)), float(np.linalg.norm(F_cond_arr)), R_global

        def solve_condensed_step_from_state(iteration: int, load_factor_state: float) -> np.ndarray:
            return runner.solve_condensed_system(
                runner.K_cond,
                runner.F_cond,
                runner.primary_results_dir,
                solver_name=runner.solver_config.get("type", "cg"),
                preconditioner="auto",
                tolerance=runner.solver_config.get("tolerance", 1e-6),
                max_iterations=runner.solver_config.get("max_iterations", 1000),
                restart=runner.solver_config.get("restart", 20),
                ilu_drop_tol=runner.solver_config.get("ilu_drop_tol", 1e-6),
                ilu_fill_factor=runner.solver_config.get("ilu_fill_factor", 1.0),
                disable_scaling=runner.solver_config.get("disable_scaling", False),
                load_increment_index=increment_index + 1,
                newton_iter=iteration,
                load_factor=float(load_factor_state),
            )

        def reconstruct_delta(delta_u_cond: np.ndarray) -> np.ndarray:
            return runner.reconstruct_global_system(
                runner.condensed_dofs,
                delta_u_cond,
                runner.total_dof,
                runner.primary_results_dir,
                fixed_dofs=runner.fixed_dofs,
                inactive_dofs=runner.inactive_dofs,
                local_global_dof_map=runner.local_global_dof_map,
            )

        def on_iteration(record):
            runner.newton_iterations_total += 1
            runner.last_norm_F_cond = record.norm_F_cond

        predictor_du = np.zeros_like(current_U)
        if current_U.size:
            predictor_du[tip_dof] = float(continuation_config.arc_length_radius)
        result = solve_arc_length_corrector(
            U_prev=current_U,
            load_factor_prev=0.0 if increment_index == 0 else float(continuation_config.load_factors[increment_index - 1]),
            predictor_displacement=predictor_du,
            reference_load_vector=reference_load_vector,
            arc_length_radius=float(continuation_config.arc_length_radius),
            alpha_scale=float(continuation_config.arc_length_alpha_scale),
            newton_tol=float(runner.newton_tol),
            newton_max_iter=int(runner.newton_max_iter),
            build_system_from_state=build_system_from_state,
            solve_condensed_step_from_state=solve_condensed_step_from_state,
            reconstruct_delta=reconstruct_delta,
            iteration_callback=on_iteration,
            load_increment_index=increment_index + 1,
        )
        runner.last_norm_F_cond = result.last_norm_F_cond
        return result

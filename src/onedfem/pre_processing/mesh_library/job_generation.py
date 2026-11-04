from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Iterable, Mapping, Sequence

import numpy as np


NONLINEAR_STATIC_SETTINGS = """[Simulation]\n[Type]\nstatic_nonlinear\n\n[Newton]\ntolerance = 1e-4\nrelative_tolerance = 1e-6\nrelative_reference = first_residual\nmax_iterations = 20\ntolerance_delta_u = 1e-9\n\n[Solver]\ntolerance = 1e-10\nmax_iterations = 2000\n\n[Nonlinear]\nnum_increments = 2\nline_search = false\n"""

LINEAR_STATIC_SETTINGS = """[Simulation]\n[Type]      \nStatic             \n#Dynamic          \n#Modal\n"""


@dataclass(frozen=True)
class MeshSpec:
    length: float
    num_elements: int
    growth_factor: float = 0.0


@dataclass(frozen=True)
class ElementSpec:
    element_type: str
    axial_order: int = 3
    bending_y_order: int = 3
    bending_z_order: int = 3
    shear_y_order: int = 0
    shear_z_order: int = 0
    torsion_order: int = 3
    load_order: int = 2
    warping: int | None = None


@dataclass(frozen=True)
class MaterialSpec:
    E: float = 2.10e11
    G: float = 8.10e10
    nu: float = 0.3
    rho: float = 7850.0


@dataclass(frozen=True)
class SectionSpec:
    A: float = 1.3075519589902397e-3
    I_x: float = 0.0
    I_y: float = 3.2340029e-7
    I_z: float = 2.0876865e-6
    J_t: float = 2.606727e-8


@dataclass(frozen=True)
class PrescribedDisplacementSpec:
    id: int
    node_id: int
    dof: str
    value: float
    constraint_type: str = "displacement"
    comment: str = "# Fixed support"


@dataclass(frozen=True)
class PointLoadSpec:
    x: float
    y: float = 0.0
    z: float = 0.0
    F_x: float = 0.0
    F_y: float = 0.0
    F_z: float = 0.0
    M_x: float = 0.0
    M_y: float = 0.0
    M_z: float = 0.0

    def components(self) -> tuple[float, ...]:
        return (self.x, self.y, self.z, self.F_x, self.F_y, self.F_z, self.M_x, self.M_y, self.M_z)


@dataclass(frozen=True)
class DistributedLoadSpec:
    samples: tuple[PointLoadSpec, ...]


@dataclass(frozen=True)
class JobSpec:
    job_name: str
    mesh: MeshSpec
    element: ElementSpec
    material: MaterialSpec = field(default_factory=MaterialSpec)
    section: SectionSpec = field(default_factory=SectionSpec)
    prescribed_displacements: tuple[PrescribedDisplacementSpec, ...] = field(default_factory=tuple)
    point_loads: tuple[PointLoadSpec, ...] = field(default_factory=tuple)
    distributed_load: DistributedLoadSpec | None = None
    simulation_settings_text: str = LINEAR_STATIC_SETTINGS
    metadata_comments: Mapping[str, str] = field(default_factory=dict)
    readme_reference_text: str | None = None


def fixed_cantilever_support() -> tuple[PrescribedDisplacementSpec, ...]:
    return tuple(
        PrescribedDisplacementSpec(id=i, node_id=0, dof=dof, value=0.0)
        for i, dof in enumerate(("UX", "UY", "UZ", "RX", "RY", "RZ"))
    )


def generate_node_positions(mesh: MeshSpec) -> np.ndarray:
    num_nodes = int(mesh.num_elements) + 1
    if mesh.growth_factor == 0:
        return np.linspace(0.0, mesh.length, num_nodes)
    i = np.linspace(0.0, 1.0, num_nodes)
    norm = (np.exp(mesh.growth_factor * i) - 1.0) / (np.exp(mesh.growth_factor) - 1.0)
    return (1.0 - norm) * mesh.length


def generate_elements(node_positions: np.ndarray) -> list[tuple[int, int]]:
    return [(idx, idx + 1) for idx in range(len(node_positions) - 1)]


def distributed_samples(num_elements: int, length: float, fy_fn: Callable[[float], float]) -> tuple[PointLoadSpec, ...]:
    return tuple(
        PointLoadSpec(x=length * i / num_elements, F_y=float(fy_fn(length * i / num_elements)))
        for i in range(num_elements + 1)
    )


def write_job_from_spec(base_dir: str | Path, spec: JobSpec) -> Path:
    job_dir = Path(base_dir) / spec.job_name
    job_dir.mkdir(parents=True, exist_ok=True)

    node_positions = generate_node_positions(spec.mesh)
    elements = generate_elements(node_positions)

    _write_grid_file(job_dir, node_positions)
    _write_element_file(job_dir, elements, spec.element)
    _write_material_file(job_dir, len(elements), spec.material)
    _write_section_file(job_dir, len(elements), spec.section)
    _write_prescribed_displacement_file(job_dir, spec.prescribed_displacements)
    _write_simulation_settings_file(job_dir, spec.simulation_settings_text)

    if spec.point_loads:
        _write_load_file(
            job_dir / "point_load.txt",
            "[Point load]",
            spec.point_loads,
            spec.metadata_comments,
        )
    if spec.distributed_load is not None:
        _write_load_file(
            job_dir / "distributed_load.txt",
            "[Distributed load]",
            spec.distributed_load.samples,
            spec.metadata_comments,
        )
    if spec.readme_reference_text is not None:
        (job_dir / "README_REFERENCE.md").write_text(spec.readme_reference_text, encoding="utf-8")

    return job_dir


def _write_grid_file(job_dir: Path, node_positions: np.ndarray) -> None:
    with open(job_dir / "grid.txt", "w", encoding="utf-8") as f:
        f.write("[Grid]\n[node_id]   [x]         [y]       [z]\n")
        for idx, x in enumerate(node_positions):
            f.write(f"{idx:<12}{x:<12.6f}{0.0:<10.1f}{0.0:<10.1f}\n")


def _write_element_file(job_dir: Path, elements: Sequence[tuple[int, int]], spec: ElementSpec) -> None:
    include_warping = spec.warping is not None
    headers = [
        "[element_id]", "[node1]", "[node2]", "[element_type]",
        "[axial_order]", "[bending_y_order]", "[bending_z_order]",
        "[shear_y_order]", "[shear_z_order]", "[torsion_order]", "[load_order]",
    ]
    if include_warping:
        headers.append("[warping]")
    with open(job_dir / "element.txt", "w", encoding="utf-8") as f:
        f.write("[Element]\n")
        f.write("  ".join(headers) + "\n")
        for eid, (n1, n2) in enumerate(elements):
            values: list[object] = [
                eid, n1, n2, spec.element_type, spec.axial_order, spec.bending_y_order,
                spec.bending_z_order, spec.shear_y_order, spec.shear_z_order,
                spec.torsion_order, spec.load_order,
            ]
            if include_warping:
                values.append(int(spec.warping))
            f.write("  ".join(str(v) for v in values) + "\n")


def _write_material_file(job_dir: Path, num_elements: int, spec: MaterialSpec) -> None:
    with open(job_dir / "material.txt", "w", encoding="utf-8") as f:
        f.write("[Material]\n[element_id]  [E]         [G]          [nu]   [rho]\n")
        for idx in range(num_elements):
            f.write(f"{idx:<14}{spec.E:<12.1e}{spec.G:<13.4e}{spec.nu:<7.1f}{int(spec.rho)}\n")


def _write_section_file(job_dir: Path, num_elements: int, spec: SectionSpec) -> None:
    with open(job_dir / "section.txt", "w", encoding="utf-8") as f:
        f.write("[Section]\n[element_id]  [A]  [I_x]  [I_y]  [I_z]  [J_t]\n")
        for idx in range(num_elements):
            f.write(f"{idx}  {spec.A:.5f}  {spec.I_x:.1f}  {spec.I_y:.5e}  {spec.I_z:.5e}  {spec.J_t:.5e}\n")


def _write_prescribed_displacement_file(job_dir: Path, rows: Iterable[PrescribedDisplacementSpec]) -> None:
    with open(job_dir / "prescribed_displacement.txt", "w", encoding="utf-8") as f:
        f.write("[Prescribed Displacement]\n")
        f.write("[id]     [node_id]  [dof]   [value]     [type]          [comment]\n")
        for row in rows:
            f.write(
                f"{row.id}        {row.node_id}          {row.dof}      {row.value:.1f}         "
                f"{row.constraint_type}     {row.comment}\n"
            )


def _write_simulation_settings_file(job_dir: Path, text: str) -> None:
    (job_dir / "simulation_settings.txt").write_text(text, encoding="utf-8")


def _write_load_file(path: Path, heading: str, rows: Sequence[PointLoadSpec], metadata_comments: Mapping[str, str]) -> None:
    with open(path, "w", encoding="utf-8") as f:
        for key, value in metadata_comments.items():
            f.write(f"# {key}: {value}\n")
        f.write(f"{heading}\n")
        f.write("         [x]          [y]          [z]        [F_x]        [F_y]        [F_z]        [M_x]        [M_y]        [M_z]\n")
        for row in rows:
            f.write(
                "    " + "     ".join(f"{value:.6f}" for value in row.components()) + "\n"
            )

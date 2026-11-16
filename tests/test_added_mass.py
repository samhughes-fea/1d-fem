"""Labelled strip-heave added mass: formula, attachment, default mass unchanged.

No still-water frequency is asserted. The plate coefficient is Newman's, not fitted.
"""

from __future__ import annotations

import os
import shutil
import tempfile

import numpy as np
import pytest
import scipy.sparse as sp

from onedfem.pre_processing.element_library.added_mass import (
    STRIP_HEAVE_FORMULA,
    STRIP_HEAVE_SOURCE,
    AddedMass,
    added_mass_from_settings,
    element_added_mass,
    strip_heave_added_mass_per_length,
)
from onedfem.processing.eigen.added_mass import add_labelled_added_mass


ELEMENT = "LinearEulerBernoulliBeamElement3D"


def _cantilever(n_elem: int = 4, L: float = 2.0, rho: float = 7850.0):
    from onedfem.pre_processing.element_library.linear.beam.zosdt.euler_bernoulli.linear_euler_bernoulli_3D import (
        LinearEulerBernoulliBeamElement3D,
    )

    n_nodes = n_elem + 1
    dx = L / n_elem
    coordinates = np.array([[i * dx, 0.0, 0.0] for i in range(n_nodes)], dtype=np.float64)
    ids = np.arange(n_nodes, dtype=np.int32)
    section = {
        "A": np.full(n_elem, 0.01),
        "I_x": np.full(n_elem, 1.0e-9),
        "I_y": np.full(n_elem, 1.0e-6),
        "I_z": np.full(n_elem, 1.0e-6),
        "J_t": np.full(n_elem, 1.0e-9),
    }
    material = {
        "E": np.full(n_elem, 210.0e9),
        "G": np.full(n_elem, 210.0e9 / 2.6),
        "nu": np.full(n_elem, 0.3),
        "rho": np.full(n_elem, rho),
    }
    element_dictionary = {
        "ids": np.arange(n_elem, dtype=np.int32),
        "connectivity": np.array([[i, i + 1] for i in range(n_elem)], dtype=np.int32),
        "types": np.array([ELEMENT] * n_elem),
        "integration_orders": {
            name: np.full(n_elem, 4, dtype=np.int64)
            for name in ("axial", "bending_y", "bending_z", "shear_y", "shear_z", "torsion", "load")
        },
    }
    grid = {"coordinates": coordinates, "ids": ids}
    point_load = np.empty((0, 9))
    distributed = np.empty((0, 9))
    temp = tempfile.mkdtemp()
    job = os.path.join(temp, "job")
    os.makedirs(os.path.join(job, "element_stiffness_matrices"), exist_ok=True)
    os.makedirs(os.path.join(job, "element_force_vectors"), exist_ok=True)
    elements, mass = [], []
    for eid in range(n_elem):
        element = LinearEulerBernoulliBeamElement3D(
            element_id=eid,
            element_dictionary=element_dictionary,
            grid_dictionary=grid,
            material_dictionary=material,
            section_dictionary=section,
            point_load_array=point_load,
            distributed_load_array=distributed,
            job_results_dir=job,
        )
        elements.append(element)
        mass.append(np.asarray(element.element_mass_matrix().M_e, dtype=np.float64))
    return {
        "elements": elements,
        "coordinates": coordinates,
        "connectivity": element_dictionary["connectivity"],
        "mass": mass,
        "temp": temp,
        "L": L,
        "n_elem": n_elem,
        "n_nodes": n_nodes,
    }


def test_strip_heave_formula_is_newman_theodorsen_plate():
    rho, c = 1025.0, 1.5
    expected = rho * np.pi * (c / 2.0) ** 2
    assert strip_heave_added_mass_per_length(rho, c) == pytest.approx(expected)
    model = AddedMass(rho_fluid=rho, chord=c)
    assert model.mass_per_length == pytest.approx(expected)
    assert model.formula == STRIP_HEAVE_FORMULA
    assert "Newman" in STRIP_HEAVE_SOURCE and "Theodorsen" in STRIP_HEAVE_SOURCE


def test_default_settings_leave_the_model_off():
    assert added_mass_from_settings({}) is None
    assert added_mass_from_settings({"added_mass": {}}) is None
    assert added_mass_from_settings({
        "added_mass": {"model": None, "rho_fluid": None, "chord": None, "heave_local_axis": "z"},
    }) is None


def test_parser_reads_added_mass_block(tmp_path):
    from onedfem.pre_processing.parsing.simulation_settings_parser import parse_simulation_settings

    path = tmp_path / "simulation_settings.txt"
    path.write_text(
        "[Simulation]\n[Type]\neigen\n\n[AddedMass]\nmodel = strip_heave\n"
        "rho_fluid = 1025\nchord = 0.4\nheave_local_axis = z\n",
        encoding="utf-8",
    )
    settings = parse_simulation_settings(str(path))
    model = added_mass_from_settings(settings)
    assert model is not None
    assert model.rho_fluid == pytest.approx(1025.0)
    assert model.chord == pytest.approx(0.4)
    assert model.heave_local_axis == "z"


def test_strip_heave_requires_rho_and_chord():
    with pytest.raises(ValueError, match="rho_fluid and chord"):
        added_mass_from_settings({"added_mass": {"model": "strip_heave", "rho_fluid": 1025.0}})
    with pytest.raises(ValueError, match="strip_heave"):
        added_mass_from_settings({"added_mass": {"model": "fitted", "rho_fluid": 1.0, "chord": 1.0}})


def test_element_added_mass_rigid_heave_is_ma_times_length():
    built = _cantilever(n_elem=2, L=2.0)
    try:
        model = AddedMass(rho_fluid=1000.0, chord=0.4, heave_local_axis="z")
        elem = built["elements"][0]
        p0, p1 = built["coordinates"][0], built["coordinates"][1]
        obj = element_added_mass(elem, p0, p1, model)
        L_e = float(np.linalg.norm(p1 - p0))
        ones_w = np.zeros(obj.M_a.shape[0])
        ones_w[2] = 1.0
        ones_w[8] = 1.0
        assert ones_w @ obj.M_a @ ones_w == pytest.approx(model.mass_per_length * L_e, rel=1e-12)
        ones_v = np.zeros_like(ones_w)
        ones_v[1] = 1.0
        ones_v[7] = 1.0
        assert abs(ones_v @ obj.M_a @ ones_v) < 1e-12 * model.mass_per_length * L_e
    finally:
        shutil.rmtree(built["temp"], ignore_errors=True)


def test_default_add_labelled_added_mass_leaves_structural_mass():
    built = _cantilever()
    try:
        ndof = built["n_nodes"] * 6
        M = sp.csr_matrix((ndof, ndof))
        for element, Me in zip(built["elements"], built["mass"]):
            dof = np.asarray(element.assemble_global_dof_indices(), dtype=int)
            M = M + sp.csr_matrix(
                (Me.ravel(), (np.repeat(dof, dof.size), np.tile(dof, dof.size))),
                shape=(ndof, ndof),
            )
        M_out, objs = add_labelled_added_mass(
            M,
            elements=built["elements"],
            coordinates=built["coordinates"],
            connectivity=built["connectivity"],
            simulation_settings={},
            total_dof=ndof,
        )
        assert objs is None
        np.testing.assert_allclose(M_out.toarray(), M.toarray(), atol=0.0, rtol=0.0)
    finally:
        shutil.rmtree(built["temp"], ignore_errors=True)


def test_labelled_model_attaches_and_increases_heave_mass():
    built = _cantilever(n_elem=4, L=2.0)
    try:
        ndof = built["n_nodes"] * 6
        M = sp.csr_matrix((ndof, ndof))
        for element, Me in zip(built["elements"], built["mass"]):
            dof = np.asarray(element.assemble_global_dof_indices(), dtype=int)
            M = M + sp.csr_matrix(
                (Me.ravel(), (np.repeat(dof, dof.size), np.tile(dof, dof.size))),
                shape=(ndof, ndof),
            )
        settings = {"added_mass": {"model": "strip_heave", "rho_fluid": 1025.0, "chord": 0.5}}
        model = added_mass_from_settings(settings)
        M_out, objs = add_labelled_added_mass(
            M,
            elements=built["elements"],
            coordinates=built["coordinates"],
            connectivity=built["connectivity"],
            simulation_settings=settings,
            total_dof=ndof,
            results_dir=os.path.join(built["temp"], "primary_results"),
        )
        assert objs is not None and len(objs) == built["n_elem"]
        ones_w = np.zeros(ndof)
        ones_w[2::6] = 1.0
        added = float(ones_w @ (M_out - M) @ ones_w)
        assert added == pytest.approx(model.mass_per_length * built["L"], rel=1e-12)
        assert (M_out - M).nnz > 0
        model_path = os.path.join(built["temp"], "primary_results", "added_mass", "added_mass_model.txt")
        text = open(model_path, encoding="utf-8").read()
        assert STRIP_HEAVE_FORMULA in text
        assert "Newman" in text
        assert "not a measured frequency" in text
    finally:
        shutil.rmtree(built["temp"], ignore_errors=True)


def test_eigen_job_added_mass_block_attaches_and_default_job_is_unchanged(tmp_path):
    """``[AddedMass]`` writes the labelled model; a sibling job without the block matches air mass."""
    import logging

    from onedfem.workflow_orchestrator.run_job import process_job

    root = os.path.join(os.path.dirname(__file__), "..")
    src = os.path.join(root, "jobs", "job_smoke_eigen")

    def _run(job_dir, extra=""):
        shutil.copytree(src, job_dir)
        if extra:
            with open(os.path.join(job_dir, "simulation_settings.txt"), "a", encoding="utf-8") as handle:
                handle.write(extra)
        res = tmp_path / (os.path.basename(job_dir) + "_results")
        for sub in ("element_stiffness_matrices", "element_force_vectors", "primary_results",
                    "secondary_results", "logs", "maps", "diagnostics"):
            os.makedirs(res / sub)
        level = logging.getLogger().level
        try:
            process_job(str(job_dir), str(res), {}, {}, force_serial=True, max_processes_per_job=1)
        finally:
            logging.getLogger().setLevel(level)
        return res

    res_default = _run(tmp_path / "job_default")
    res_labelled = _run(
        tmp_path / "job_added",
        "\n[AddedMass]\nmodel = strip_heave\nrho_fluid = 1025\nchord = 0.4\nheave_local_axis = z\n",
    )
    def _freq(res):
        files = list((res / "primary_results" / "modal_results").glob("*_frequencies.txt"))
        if not files:
            files = list((res / "primary_results").rglob("*frequencies*"))
        assert files, sorted(p.name for p in (res / "primary_results").rglob("*"))
        return np.loadtxt(files[0])

    f0 = np.sort(np.atleast_1d(_freq(res_default)).astype(float))
    f1 = np.sort(np.atleast_1d(_freq(res_labelled)).astype(float))
    n = min(f0.size, f1.size)
    assert np.all(f0[:n] > 0.0) and np.all(f1[:n] > 0.0)
    # M_a is positive semi-definite: sorted frequencies cannot rise. Heave modes fall.
    # This is not a measured frequency in water.
    np.testing.assert_array_less(f1[:n] - 1e-12, f0[:n])
    assert np.any(f1[:n] < f0[:n] - 1e-9)
    assert (res_labelled / "primary_results" / "added_mass" / "added_mass_model.txt").is_file()
    assert not (res_default / "primary_results" / "added_mass").exists()

"""Mesh-wide nodal DOF layout: stride, names, and element DOF maps on any stride."""

from __future__ import annotations

import os
import shutil
import tempfile

import numpy as np
import pytest

from pre_processing.element_library.dof_layout import (
    dof_name_map,
    dof_names,
    mesh_dof_per_node,
    mesh_extra_dofs,
    row_extra_dofs,
)
from pre_processing.parsing.prescribed_displacement_parser import parse_prescribed_displacement

TS = "LinearTimoshenkoBeamElement3D"
ORDERS = {k: 2 for k in ("axial", "bending_y", "bending_z", "shear_y", "shear_z", "torsion", "load")}


def _ed(n_el, types=TS, warping=None, extra=None):
    ed = {
        "ids": np.arange(n_el),
        "connectivity": np.column_stack([np.arange(n_el), np.arange(1, n_el + 1)]),
        "types": np.array([types] * n_el) if isinstance(types, str) else np.array(types),
        "integration_orders": {k: np.full(n_el, v) for k, v in ORDERS.items()},
    }
    if warping is not None:
        ed["warping"] = np.asarray(warping, dtype=np.int8)
    if extra is not None:
        ed["extra_dofs"] = np.asarray(extra, dtype=int)
    return ed


def test_stride_from_rows():
    assert mesh_dof_per_node(None) == 6
    assert mesh_dof_per_node(_ed(2)) == 6
    assert mesh_dof_per_node(_ed(2, warping=[0, 1])) == 7
    ed = _ed(3, warping=[1, 0, 0], extra=[0, 4, 2])
    assert mesh_extra_dofs(ed) == 4 and mesh_dof_per_node(ed) == 10
    assert [row_extra_dofs(ed, i) for i in range(3)] == [1, 4, 2]


def test_names():
    assert dof_names(6) == ["UX", "UY", "UZ", "RX", "RY", "RZ"]
    assert dof_names(7)[6] == "W"
    assert dof_names(9)[6:] == ["Q1", "Q2", "Q3"]
    m = dof_name_map(9)
    assert m["W"] == m["CHI"] == m["Q1"] == 6 and m["Q3"] == 8
    assert "W" not in dof_name_map(6)
    with pytest.raises(ValueError):
        dof_name_map(5)


def test_prescribed_displacement_on_any_stride(tmp_path):
    p = tmp_path / "prescribed_displacement.txt"
    p.write_text(
        "[Prescribed Displacement]\n[id] [node_id] [dof] [value] [type] [comment]\n"
        "0 2 UZ 0.0 displacement\n1 2 W 0.0 displacement\n2 3 Q3 0.5 displacement\n",
        encoding="utf-8",
    )
    out = parse_prescribed_displacement(str(p), dof_per_node=9)
    np.testing.assert_array_equal(out["global_dof"], [2 * 9 + 2, 2 * 9 + 6, 3 * 9 + 8])
    out7 = parse_prescribed_displacement(str(p), dof_per_node=7)
    # Q3 does not exist on a 7-DOF mesh and is skipped
    np.testing.assert_array_equal(out7["global_dof"], [2 * 7 + 2, 2 * 7 + 6])


def _elements(ed, coords):
    from pre_processing.element_library.linear.beam.first_order_shear_deformation_theory.timoshenko.linear_timoshenko_3D import (
        LinearTimoshenkoBeamElement3D,
    )

    n_el = len(ed["ids"])
    jr = os.path.join(tempfile.mkdtemp(), "job")
    for sub in ("element_stiffness_matrices", "element_force_vectors", "logs"):
        os.makedirs(os.path.join(jr, sub), exist_ok=True)
    mat = {k: np.full(n_el, v) for k, v in dict(E=2.1e11, G=8.1e10, nu=0.3, rho=0.0).items()}
    sec = {k: np.full(n_el, v) for k, v in dict(A=1e-3, I_x=0.0, I_y=1e-6, I_z=2e-6, J_t=1e-7, kappa=5 / 6,
                                                alpha=0.0, y_sc=0.0, z_sc=0.0, Gamma=1e-9).items()}
    els = []
    for e in range(n_el):
        el = LinearTimoshenkoBeamElement3D(
            element_id=e, element_dictionary=ed, grid_dictionary={"coordinates": coords},
            material_dictionary=mat, section_dictionary=sec,
            point_load_array=np.empty((0, 9)), distributed_load_array=np.empty((0, 9)),
            job_results_dir=jr,
        )
        el.logger_operator = None
        els.append(el)
    return els, jr


@pytest.mark.parametrize("warping, extra", [(None, None), ([0, 0, 0, 0], [0, 3, 0, 0]), ([1, 1, 1, 1], [2, 0, 0, 0])])
def test_cantilever_is_the_same_on_any_stride(warping, extra):
    """A plain or warping Timoshenko cantilever gives the same tip response on a wider stride."""
    n_el, L = 4, 2.0
    x = np.linspace(0, L, n_el + 1)
    coords = np.column_stack([x, 0 * x, 0 * x])

    def solve(ed):
        els, jr = _elements(ed, coords)
        try:
            dpn = mesh_dof_per_node(ed)
            n = dpn * (n_el + 1)
            K = np.zeros((n, n))
            touched = np.zeros(n, bool)
            for el in els:
                assert el.mesh_dof_per_node == dpn
                g = el.assemble_global_dof_indices()
                K[np.ix_(g, g)] += el.element_stiffness_matrix().K_e
                touched[g] = True
            F = np.zeros(n)
            F[dpn * n_el:dpn * n_el + 6] = [1e3, 50.0, 80.0, 10.0, 0.0, 0.0]
            fixed = set(range(dpn)) | set(np.flatnonzero(~touched))
            if warping is not None and dpn > 6:
                fixed |= {dpn * i + 6 for i in range(n_el + 1)} if not any(warping) else set()
            free = np.array(sorted(set(range(n)) - fixed))
            U = np.zeros(n)
            U[free] = np.linalg.solve(K[np.ix_(free, free)], F[free])
            return U.reshape(n_el + 1, dpn)
        finally:
            shutil.rmtree(os.path.dirname(jr), ignore_errors=True)

    base_ed = _ed(n_el, warping=warping)
    wide_ed = _ed(n_el, warping=warping, extra=extra)
    U_base, U_wide = solve(base_ed), solve(wide_ed)
    k = U_base.shape[1]
    np.testing.assert_allclose(U_wide[:, :k], U_base, rtol=1e-12, atol=1e-15)
    assert np.all(U_wide[:, k:] == 0.0)

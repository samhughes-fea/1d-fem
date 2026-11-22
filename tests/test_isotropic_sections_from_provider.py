"""Isotropic sections come from the section provider (ASM), not from 1D-FEM's own formulas.

Every isotropic element (bar, truss, Euler-Bernoulli, Timoshenko, their nonlinear and
corotational variants) takes its section ``D`` and its section mass from
``SectionAnalysis.from_isotropic``, through ``isotropic_section``. The checks:

- each element's ``D`` is the provider's (sliced to the element's strains);
- the consistent mass uses the provider's mass per length, with the polar moment
  ``I_p`` on the twist DOFs (``I_x <= 0`` in ``section.txt`` means ``I_y + I_z``);
- swapping the provider for one that doubles every stiffness doubles every element ``D``,
  so no isotropic element keeps a private copy of the formulas.
"""

from __future__ import annotations

import os
import sys
import types

import numpy as np
import pytest

from anisotropic_section_model.section_analysis import SectionAnalysis

from onedfem import sections
from onedfem.pre_processing.element_library.isotropic_section import (
    _cached,
    isotropic_section,
    nodal_mass_diagonal,
)

E, G, RHO = 2.1e11, 8.1e10, 7850.0
SEC = dict(A=1e-3, I_x=0.0, I_y=1e-8, I_z=2e-8, J_t=1e-9, kappa=0.8, alpha=0.0)


def _kwargs(job_results_dir, etype, **section):
    for sub in ("element_stiffness_matrices", "element_force_vectors", "logs"):
        os.makedirs(os.path.join(job_results_dir, sub), exist_ok=True)
    ed = {
        "ids": np.array([0]),
        "connectivity": np.array([[0, 1]]),
        "types": np.array([etype]),
        "integration_orders": {k: np.array([2]) for k in
                               ("axial", "bending_y", "bending_z", "shear_y", "shear_z", "torsion", "load")},
    }
    return dict(
        element_id=0,
        element_dictionary=ed,
        grid_dictionary={"coordinates": np.array([[0.0, 0, 0], [1.0, 0, 0]])},
        material_dictionary={k: np.array([v]) for k, v in dict(E=E, G=G, nu=0.3, rho=RHO).items()},
        section_dictionary={k: np.array([v]) for k, v in {**SEC, **section}.items()},
        point_load_array=np.empty((0, 9)),
        distributed_load_array=np.empty((0, 9)),
        job_results_dir=job_results_dir,
    )


def _timoshenko():
    from onedfem.pre_processing.element_library.linear.beam.fosdt.timoshenko.linear_timoshenko_3D import (
        LinearTimoshenkoBeamElement3D,
    )
    return LinearTimoshenkoBeamElement3D


def _euler_bernoulli():
    from onedfem.pre_processing.element_library.linear.beam.zosdt.euler_bernoulli.linear_euler_bernoulli_3D import (
        LinearEulerBernoulliBeamElement3D,
    )
    return LinearEulerBernoulliBeamElement3D


def _bar():
    from onedfem.pre_processing.element_library.linear.bar.linear_bar_3D import LinearBarElement3D
    return LinearBarElement3D


def _truss():
    from onedfem.pre_processing.element_library.linear.truss.linear_truss_3D import LinearTrussElement3D
    return LinearTrussElement3D


def _reference(kappa=SEC["kappa"]):
    return SectionAnalysis.from_isotropic(E=E, G=G, A=SEC["A"], I_y=SEC["I_y"], I_z=SEC["I_z"],
                                          J_t=SEC["J_t"], kappa=kappa, rho=RHO)


def test_beam_elements_take_D_from_the_provider(tmp_path):
    ref = _reference()
    tim = _timoshenko()(**_kwargs(str(tmp_path / "t"), "LinearTimoshenkoBeamElement3D"))
    np.testing.assert_array_equal(tim.D_matrix(), ref.stiffness)
    np.testing.assert_array_equal(tim.material_stiffness_operator.assembly_form(), ref.stiffness)
    eb = _euler_bernoulli()(**_kwargs(str(tmp_path / "e"), "LinearEulerBernoulliBeamElement3D"))
    np.testing.assert_array_equal(eb.D_matrix(), ref.euler_bernoulli_D())
    np.testing.assert_array_equal(eb.material_stiffness_operator.assembly_form(), ref.euler_bernoulli_D())


def test_bar_and_truss_take_D_from_the_provider(tmp_path):
    ref = _reference()
    bar = _bar()(**_kwargs(str(tmp_path / "b"), "LinearBarElement3D"))
    np.testing.assert_array_equal(bar.material_stiffness_operator.assembly_form(),
                                  ref.stiffness[np.ix_([0, 5], [0, 5])])
    truss = _truss()(**_kwargs(str(tmp_path / "t"), "LinearTrussElement3D"))
    np.testing.assert_array_equal(truss.material_stiffness_operator.assembly_form(),
                                  ref.stiffness[np.ix_([0, 3, 5], [0, 3, 5])])


def test_shear_centre_offset_is_recorded_not_coupled(tmp_path):
    """``y_sc``, ``z_sc`` reach the provider as the shear centre; ``D`` stays diagonal."""
    kw = _kwargs(str(tmp_path / "t"), "LinearTimoshenkoBeamElement3D", y_sc=0.01, z_sc=-0.02)
    el = _timoshenko()(**kw)
    D = el.D_matrix()
    np.testing.assert_array_equal(D, np.diag(np.diag(D)))
    assert el.section_analysis.shear_centre == pytest.approx((0.01, -0.02))


@pytest.mark.parametrize("I_x, I_p", [(0.0, SEC["I_y"] + SEC["I_z"]), (5e-8, 5e-8)])
def test_twist_mass_is_the_polar_moment(tmp_path, I_x, I_p):
    for cls, etype in ((_timoshenko(), "LinearTimoshenkoBeamElement3D"),
                       (_euler_bernoulli(), "LinearEulerBernoulliBeamElement3D")):
        el = cls(**_kwargs(str(tmp_path / etype), etype, I_x=I_x))
        assert el.isotropic_section_analysis.mass["I_p"] == pytest.approx(RHO * I_p)
        assert el.I_p == pytest.approx(I_p)
        assert el.J_t == SEC["J_t"]                     # the torsion constant is a separate input
        mu = nodal_mass_diagonal(el.isotropic_section_analysis)
        np.testing.assert_allclose(mu[[0, 3, 4, 5]], RHO * np.array([SEC["A"], I_p, SEC["I_y"], SEC["I_z"]]))
        M = el.element_mass_matrix().M_e
        # twist is linear on both elements: the consistent twist mass is rho I_p L (2 1; 1 2) / 6
        np.testing.assert_allclose(M[np.ix_([3, 9], [3, 9])], RHO * I_p / 6 * np.array([[2, 1], [1, 2]]),
                                   rtol=1e-12)


def test_bar_mass_is_the_provider_mass_per_length(tmp_path):
    bar = _bar()(**_kwargs(str(tmp_path / "b"), "LinearBarElement3D"))
    M = bar.element_mass_matrix().M_e
    # axial translation is linear: rho A L (2 1; 1 2) / 6
    np.testing.assert_allclose(M[np.ix_([0, 6], [0, 6])], RHO * SEC["A"] / 6 * np.array([[2, 1], [1, 2]]),
                               rtol=1e-12)


def test_isotropic_section_caches_and_reads_I_x():
    a = isotropic_section(E=E, G=G, A=1e-3, I_y=1e-8, I_z=2e-8, J_t=1e-9, I_x=0.0, rho=1.0)
    b = isotropic_section(E=E, G=G, A=1e-3, I_y=1e-8, I_z=2e-8, J_t=1e-9, I_x=None, rho=1.0)
    assert a is b
    assert a.mass["I_p"] == pytest.approx(3e-8)


@pytest.fixture
def doubled_provider():
    """A provider that is ASM with every isotropic stiffness doubled."""
    import anisotropic_section_model as asm
    import anisotropic_section_model.section_analysis as asm_sa

    root = types.ModuleType("doubled_asm")
    root.__path__ = []                                   # a package, for section_module
    root.SECTION_API_VERSION = getattr(asm, "SECTION_API_VERSION", None)
    mod = types.ModuleType("doubled_asm.section_analysis")

    class Doubled(asm_sa.SectionAnalysis):
        @classmethod
        def from_isotropic(cls, **kw):
            base = asm_sa.SectionAnalysis.from_isotropic(**kw)
            return asm_sa.SectionAnalysis(2.0 * base.stiffness, 2.0 * base.classical_stiffness,
                                          shear_centre=base.shear_centre, mass=base.mass, source="isotropic")

    mod.SectionAnalysis = Doubled
    sys.modules["doubled_asm"] = root
    sys.modules["doubled_asm.section_analysis"] = mod
    previous = sections._provider
    sections.set_section_provider(root)
    try:
        yield
    finally:
        sections.set_section_provider(previous)
        sys.modules.pop("doubled_asm", None)
        sys.modules.pop("doubled_asm.section_analysis", None)
        _cached.cache_clear()


def test_every_isotropic_element_follows_the_provider(tmp_path, doubled_provider):
    ref = _reference()
    tim = _timoshenko()(**_kwargs(str(tmp_path / "t"), "LinearTimoshenkoBeamElement3D"))
    np.testing.assert_allclose(tim.D_matrix(), 2.0 * ref.stiffness)
    np.testing.assert_allclose(tim.material_stiffness_operator.assembly_form(), 2.0 * ref.stiffness)
    eb = _euler_bernoulli()(**_kwargs(str(tmp_path / "e"), "LinearEulerBernoulliBeamElement3D"))
    np.testing.assert_allclose(eb.D_matrix(), 2.0 * ref.euler_bernoulli_D())
    bar = _bar()(**_kwargs(str(tmp_path / "b"), "LinearBarElement3D"))
    np.testing.assert_allclose(bar.material_stiffness_operator.assembly_form(),
                               2.0 * ref.stiffness[np.ix_([0, 5], [0, 5])])
    truss = _truss()(**_kwargs(str(tmp_path / "r"), "LinearTrussElement3D"))
    np.testing.assert_allclose(truss.material_stiffness_operator.assembly_form(),
                               2.0 * ref.stiffness[np.ix_([0, 3, 5], [0, 3, 5])])

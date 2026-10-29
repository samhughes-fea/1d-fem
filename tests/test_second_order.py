"""Second-order geometric stiffness: buckling and P-delta against closed forms.

References (closed forms with fork supports unless noted):

- Euler and Engesser columns.
- Flexural-torsional buckling of a monosymmetric section in compression (Timoshenko and
  Gere, *Theory of Elastic Stability*, 1961, sec. 5.4), with the shear centre at ``y_0``
  from the centroid on the symmetry axis::

      (P_w - P)(P_phi - P) - P^2 y_0^2 / r_0^2 = 0,
      P_w = pi^2 E I_y / L^2,  P_phi = (GJ + pi^2 E Gamma / L^2) / r_0^2,
      r_0^2 = (I_y + I_z) / A + y_0^2.

- Lateral-torsional buckling under uniform moment (TG sec. 6.3)::

      M_cr = (pi / L) sqrt(E I_z GJ (1 + pi^2 E Gamma / (GJ L^2)))

  and of a monosymmetric I-beam with the Wagner coefficient
  ``beta = (1/I_y) int z (y^2 + z^2) dA - 2 z_0`` (Trahair, *Flexural-Torsional Buckling of
  Structures*, 1993, eq. 6.3)::

      M_cr = P_z (beta/2 +- sqrt(beta^2/4 + (Gamma/I_z)(1 + GJ L^2/(pi^2 E Gamma)))),
      P_z = pi^2 E I_z / L^2.

- Narrow rectangular beams with ``Gamma = 0`` (TG sec. 6.4 and 6.5): a cantilever with a
  tip load at the centroid buckles at ``4.013 sqrt(E I_z GJ) / L^2``; a simply supported
  beam with a central load at the centroid at ``16.94 sqrt(E I_z GJ) / L^2``.
"""

from __future__ import annotations

import importlib
import importlib.util
import os
import shutil
import tempfile

import numpy as np
import pytest
from scipy import linalg

from pre_processing.element_library.dof_layout import mesh_dof_per_node
from pre_processing.element_library.element_factory import ElementFactory
from pre_processing.element_library.linear.beam.second_order import geometric_section_matrix
from pre_processing.section_integration.section_analysis import SectionAnalysis
from pre_processing.section_integration.thin_walled import (
    Ply,
    PlyMaterial,
    ThinWalledSection,
    ThinWalledSectionAnalysis,
    Wall,
    WallLaminate,
)

E, NU = 2.0e11, 0.3
G = E / (2 * (1 + NU))
EB, TS, BS, REDDY = (
    "LinearEulerBernoulliBeamElement3D",
    "LinearTimoshenkoBeamElement3D",
    "LinearBenscoterBeamElement3D",
    "LinearReddyBeamElement3D",
)
CFRP = PlyMaterial(E1=140e9, E2=10e9, G12=5e9, nu12=0.3, G13=5e9, G23=3.5e9, rho=1600.0)


# --------------------------------------------------------------------- #
# Sections
# --------------------------------------------------------------------- #
def iso(t):
    return WallLaminate.isotropic(E, NU, t)


def channel(h=0.2, b=0.08, t=4e-3, shift=(0.0, 0.0), lam=None):
    """Channel, web on ``y = 0`` along ``z``, flanges towards ``+y``; origin moved by ``-shift``."""
    dy, dz = shift
    pts = {1: (b, -h / 2), 2: (0.0, -h / 2), 3: (0.0, h / 2), 4: (b, h / 2)}
    pts = {k: (y - dy, z - dz) for k, (y, z) in pts.items()}
    walls = [Wall(1, 1, 2, 1), Wall(2, 2, 3, 1), Wall(3, 3, 4, 1)]
    return ThinWalledSectionAnalysis(ThinWalledSection(pts, walls), {1: lam or iso(t)}, hoop="free", n_sub=2)


def channel_properties(h=0.2, b=0.08, t=4e-3):
    """Mid-line formulas: centroid, shear centre (from the web), I_y, I_z, J, Gamma."""
    A = t * (h + 2 * b)
    y_c = b * b / (h + 2 * b)
    e = 3 * b * b / (h + 6 * b)              # shear centre at y = -e
    I_y = t * h**3 / 12 + 2 * b * t * (h / 2) ** 2
    I_z = 2 * t * b**3 / 3 - A * y_c**2
    J = t**3 * (h + 2 * b) / 3
    Gamma = t * b**3 * h**2 / 12 * (3 * b + 2 * h) / (6 * b + h)
    return dict(A=A, y_c=y_c, y_s=-e, I_y=I_y, I_z=I_z, J=J, Gamma=Gamma)


def i_beam(h=0.3, b_top=0.15, b_bot=0.15, tf=6e-3, tw=4e-3, shift=0.0):
    pts = {1: (-b_top / 2, h / 2), 2: (0.0, h / 2), 3: (b_top / 2, h / 2),
           4: (-b_bot / 2, -h / 2), 5: (0.0, -h / 2), 6: (b_bot / 2, -h / 2)}
    pts = {k: (y, z - shift) for k, (y, z) in pts.items()}
    walls = [Wall(1, 1, 2, 1), Wall(2, 2, 3, 1), Wall(3, 4, 5, 1), Wall(4, 5, 6, 1), Wall(5, 5, 2, 2)]
    return ThinWalledSectionAnalysis(ThinWalledSection(pts, walls), {1: iso(tf), 2: iso(tw)}, hoop="free", n_sub=2)


def i_beam_properties(h=0.3, b_top=0.15, b_bot=0.15, tf=6e-3, tw=4e-3):
    """Mid-line properties about the centroid (z up), and the Wagner coefficient ``beta``."""
    A_t, A_b, A_w = b_top * tf, b_bot * tf, h * tw
    A = A_t + A_b + A_w
    z_c = (A_t - A_b) * h / 2 / A
    zt, zb = h / 2 - z_c, -h / 2 - z_c
    I_y = A_t * zt**2 + A_b * zb**2 + tw * h**3 / 12 + A_w * z_c**2
    It, Ib = tf * b_top**3 / 12, tf * b_bot**3 / 12
    I_z = It + Ib
    z_s = (It * zt + Ib * zb) / (It + Ib)    # shear centre
    Gamma = It * Ib / (It + Ib) * h**2
    J = (b_top * tf**3 + b_bot * tf**3 + h * tw**3) / 3
    # int z (y^2 + z^2) dA over the mid-lines
    flange = lambda b, z: tf * (b**3 / 12 * z + b * z**3)
    web = tw * ((zt**4 - zb**4) / 4)
    beta = (flange(b_top, zt) + flange(b_bot, zb) + web) / I_y - 2 * z_s
    return dict(A=A, z_c=z_c, z_s=z_s, I_y=I_y, I_z=I_z, J=J, Gamma=Gamma, beta=beta)


# --------------------------------------------------------------------- #
# Direct assembly
# --------------------------------------------------------------------- #
def _ed(etype, n_el, warping):
    o = {"axial": 2, "bending_y": 2, "bending_z": 2, "shear_y": 1, "shear_z": 1, "torsion": 2, "load": 2}
    return {
        "ids": np.arange(n_el),
        "connectivity": np.column_stack([np.arange(n_el), np.arange(1, n_el + 1)]),
        "types": np.array([etype] * n_el),
        "warping": np.full(n_el, int(warping), dtype=np.int8),
        "integration_orders": {k: np.full(n_el, v) for k, v in o.items()},
    }


class Beam:
    """A straight beam along ``x`` assembled directly from one element type."""

    def __init__(self, etype, L, n_el, sa=None, rect=None):
        self.etype, self.L, self.n_el = etype, L, n_el
        Cls = getattr(importlib.import_module(ElementFactory.ELEMENT_CLASS_MAP[etype]), etype)
        ed = _ed(etype, n_el, etype == BS)
        if sa is not None:
            ed["section_analysis_per_element"] = tuple([SectionAnalysis.from_thin_walled(sa)] * n_el)
        b, h = rect if rect is not None else (0.01, 0.2)
        A, Iy, Iz = b * h, b * h**3 / 12, h * b**3 / 12
        J = h * b**3 / 3 * (1 - 0.63 * b / h)
        self.rect = dict(A=A, I_y=Iy, I_z=Iz, J=J)
        mat = {k: np.full(n_el, v) for k, v in dict(E=E, G=G, nu=NU, rho=0.0).items()}
        props = dict(A=A, I_x=Iy + Iz, I_y=Iy, I_z=Iz, J_t=J, kappa=5 / 6, alpha=0.0)
        if etype == BS:  # the warping mesh reads y_sc, z_sc and Gamma (unused with a thin-walled section)
            props.update(y_sc=0.0, z_sc=0.0, Gamma=1e-9)
        sec = {k: np.full(n_el, v) for k, v in props.items()}
        x = np.linspace(0.0, L, n_el + 1)
        coords = np.column_stack([x, 0 * x, 0 * x])
        self._tmp = tempfile.mkdtemp()
        jr = os.path.join(self._tmp, "job")
        for sub in ("element_stiffness_matrices", "element_force_vectors", "logs"):
            os.makedirs(os.path.join(jr, sub), exist_ok=True)
        self.dpn = mesh_dof_per_node(ed)
        self.n = self.dpn * (n_el + 1)
        self.elements, self.maps = [], []
        self.K = np.zeros((self.n, self.n))
        for e in range(n_el):
            el = Cls(element_id=e, element_dictionary=ed, grid_dictionary={"coordinates": coords},
                     material_dictionary=mat, section_dictionary=sec, point_load_array=np.empty((0, 9)),
                     distributed_load_array=np.empty((0, 9)), job_results_dir=jr)
            el.logger_operator = None
            g = np.asarray(el.assemble_global_dof_indices())
            self.K[np.ix_(g, g)] += el.element_stiffness_matrix().K_e
            self.elements.append(el)
            self.maps.append(g)
        shutil.rmtree(self._tmp, ignore_errors=True)

    def dof(self, node, local):
        return self.dpn * node + local

    def K_sigma(self, U):
        Kg = np.zeros((self.n, self.n))
        for el, g in zip(self.elements, self.maps):
            Kg[np.ix_(g, g)] += el.linear_geometric_stiffness_matrix(U[g])
        return Kg

    def fork(self, axial_root=True, mid=False):
        """Fork supports at both ends (``v = w = theta_x = 0``), ``u`` fixed at the root."""
        fixed = [self.dof(0, i) for i in (1, 2, 3)] + [self.dof(self.n_el, i) for i in (1, 2, 3)]
        if axial_root:
            fixed.append(self.dof(0, 0))
        return fixed

    def clamp_root(self):
        return list(range(self.dpn))

    def static(self, F, fixed, K_extra=None):
        free = np.setdiff1d(np.arange(self.n), fixed)
        K = self.K if K_extra is None else self.K + K_extra
        U = np.zeros(self.n)
        U[free] = np.linalg.solve(K[np.ix_(free, free)], F[free])
        return U

    def buckle(self, F, fixed, n_modes=1):
        """Load factors ``lambda`` with ``(K + lambda K_sigma(U_F)) x = 0`` and their modes."""
        U = self.static(F, fixed)
        Kg = self.K_sigma(U)
        free = np.setdiff1d(np.arange(self.n), fixed)
        ev, vec = linalg.eig(self.K[np.ix_(free, free)], -Kg[np.ix_(free, free)])
        ok = np.isfinite(ev) & (np.abs(ev.imag) < 1e-8 * np.abs(ev.real)) & (ev.real > 0)
        order = np.argsort(ev.real[ok])
        lam = ev.real[ok][order][:n_modes]
        modes = np.zeros((self.n, len(lam)))
        modes[free] = vec[:, ok][:, order][:, :n_modes].real
        return lam, modes


def axial(beam, P, y=0.0, z=0.0):
    """End compression ``P`` acting at ``(y, z)`` of both end sections (root reaction at ``u = 0``)."""
    F = np.zeros(beam.n)
    F[beam.dof(beam.n_el, 0)] = -P
    for node, sign in ((beam.n_el, 1.0), (0, -1.0)):
        F[beam.dof(node, 4)] = sign * z * -P      # M_y = z F_x
        F[beam.dof(node, 5)] = sign * -y * -P     # M_z = -y F_x
    return F


def end_moments(beam, M):
    """Uniform right-handed ``M_y`` along the beam."""
    F = np.zeros(beam.n)
    F[beam.dof(0, 4)] = -M
    F[beam.dof(beam.n_el, 4)] = M
    return F


# --------------------------------------------------------------------- #
# Section terms
# --------------------------------------------------------------------- #
def test_geometric_terms_resultants_and_polar_term():
    """N, S_y = -M_z, S_z = M_y and V come back exactly; W/N is I_p/A about the origin."""
    p = channel_properties()
    an = channel()
    H = an.geometric_terms
    np.testing.assert_allclose(H[[0, 1, 2, 4, 5]], np.array([
        [1, 0, 0, 0, 0, 0], [0, 0, -1, 0, 0, 0], [0, 1, 0, 0, 0, 0], [0, 0, 0, 1, 0, 0], [0, 0, 0, 0, 1, 0],
    ], dtype=float), atol=1e-9)
    # an axial force at the centroid: W = N (I_y + I_z + A y_c^2) / A about the web
    F = np.array([1.0, 0.0, -p["y_c"], 0.0, 0.0, 0.0])
    assert (H @ F)[3] == pytest.approx((p["I_y"] + p["I_z"]) / p["A"] + p["y_c"] ** 2, rel=2e-3)


def test_geometric_terms_follow_the_reference_axis():
    """Moving the origin moves N r^2 by the parallel-axis rule and leaves the stresses alone."""
    a0, a1 = channel(), channel(shift=(0.03, -0.02))
    F0 = np.array([1.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    c0 = a0.geometric_terms @ F0
    # the same stresses (an axial force at the old origin) seen from the new origin
    F1 = np.array([1.0, 0.02, 0.03, 0.0, 0.0, 0.0])
    c1 = a1.geometric_terms @ F1
    assert c1[0] == pytest.approx(c0[0])
    assert c1[1] == pytest.approx(c0[1] - 0.03, abs=1e-9)
    assert c1[2] == pytest.approx(c0[2] + 0.02, abs=1e-9)
    assert c1[3] == pytest.approx(c0[3] - 2 * (0.03 * c0[1] - 0.02 * c0[2]) + (0.03**2 + 0.02**2), rel=1e-6)


def test_section_matrix_is_the_rigid_section_energy():
    """g^T S g equals the integral of sigma eps2 + q gamma2 for random terms."""
    rng = np.random.default_rng(1)
    c = rng.normal(size=7)
    g = rng.normal(size=4)
    v1, w1, p1, p0 = g
    N, Sy, Sz, W, Vy, Vz, Qr = c
    energy = N * (v1**2 + w1**2) - 2 * Sz * v1 * p1 + 2 * Sy * w1 * p1 + W * p1**2 \
        + 2 * (-Vz * v1 * p0 + Vy * w1 * p0 + Qr * p0 * p1)
    assert g @ geometric_section_matrix(c) @ g == pytest.approx(energy)


# --------------------------------------------------------------------- #
# Columns
# --------------------------------------------------------------------- #
@pytest.mark.parametrize("etype", [EB, TS, REDDY])
def test_pinned_column_euler_and_engesser(etype):
    L, n_el = 3.0, 48
    beam = Beam(etype, L, n_el, rect=(0.06, 0.1))
    r = beam.rect
    lam, _ = beam.buckle(axial(beam, 1.0), beam.fork())
    Pe = np.pi**2 * E * r["I_z"] / L**2      # weak axis: bending in x-y
    ref = Pe if etype == EB else Pe / (1 + Pe / (5 / 6 * G * r["A"]))
    assert lam[0] == pytest.approx(ref, rel=1e-3)


# --------------------------------------------------------------------- #
# Flexural-torsional buckling (shear-centre and tension-centre offsets)
# --------------------------------------------------------------------- #
def _ft_closed_form(L, Gamma):
    p = channel_properties()
    y0 = p["y_s"] - p["y_c"]
    r02 = (p["I_y"] + p["I_z"]) / p["A"] + y0**2
    Pw = np.pi**2 * E * p["I_y"] / L**2
    Pphi = (G * p["J"] + np.pi**2 * E * Gamma / L**2) / r02
    return float(min(np.roots([1 - y0**2 / r02, -(Pw + Pphi), Pw * Pphi]))), p


@pytest.mark.parametrize("origin", ["web", "centroid", "shear_centre"])
@pytest.mark.parametrize("etype", [TS, BS])
def test_channel_flexural_torsional_buckling(etype, origin):
    """Monosymmetric channel in compression, fork supports, any reference axis (TG 5.4)."""
    L = 4.0
    _, p = _ft_closed_form(L, 0.0)
    dy = {"web": 0.0, "centroid": p["y_c"], "shear_centre": p["y_s"]}[origin]
    beam = Beam(etype, L, 32, sa=channel(shift=(dy, 0.0)))
    # Timoshenko warps freely (Gamma = 0); Benscoter carries E Gamma
    ref, _ = _ft_closed_form(L, p["Gamma"] if etype == BS else 0.0)
    lam, modes = beam.buckle(axial(beam, 1.0, y=p["y_c"] - dy), beam.fork(), n_modes=3)
    # flexure along the symmetry axis (v) is uncoupled; the flexural-torsional mode twists
    twists = [np.abs(modes[3::beam.dpn, k]).max() > 1e-6 * np.abs(modes[:, k]).max() for k in range(3)]
    assert lam[twists.index(True)] == pytest.approx(ref, rel=3e-3)


def test_benscoter_warping_pole_is_the_shear_centre_for_any_origin():
    """E Gamma of a channel does not depend on where the section origin is."""
    p = channel_properties()
    for dy in (0.0, p["y_c"], p["y_s"], -0.05):
        an = channel(shift=(dy, 0.0))
        assert an.benscoter_stiffness[7, 7] == pytest.approx(E * p["Gamma"], rel=2e-3)


# --------------------------------------------------------------------- #
# Lateral-torsional buckling (Wagner term, moment gradient)
# --------------------------------------------------------------------- #
def _ltb_closed_form(L, Gamma, **kw):
    p = i_beam_properties(**kw)
    Pz = np.pi**2 * E * p["I_z"] / L**2
    T = G * p["J"] + np.pi**2 * E * Gamma / L**2
    b = p["beta"]
    return [Pz * b / 2 + s * np.sqrt((Pz * b / 2) ** 2 + Pz * T) for s in (1.0, -1.0)], p


@pytest.mark.parametrize("etype", [EB, TS, BS])
def test_i_beam_lateral_torsional_buckling_uniform_moment(etype):
    L = 6.0
    p = i_beam_properties()
    (ref, _), _ = _ltb_closed_form(L, p["Gamma"] if etype == BS else 0.0)
    beam = Beam(etype, L, 32, sa=i_beam())
    lam, _ = beam.buckle(end_moments(beam, 1.0), beam.fork())
    assert lam[0] == pytest.approx(ref, rel=2e-3)


def test_monosymmetric_i_beam_wagner_term_both_senses():
    """Large flange in compression or in tension: the two critical moments of Trahair eq. 6.3."""
    L, kw = 6.0, dict(b_top=0.2, b_bot=0.1)
    p = i_beam_properties(**kw)
    (m_plus, m_minus), _ = _ltb_closed_form(L, p["Gamma"], **kw)
    beam = Beam(BS, L, 32, sa=i_beam(shift=0.05, **kw))   # origin off the centroid and shear centre
    found = sorted([beam.buckle(end_moments(beam, s), beam.fork())[0][0] for s in (1.0, -1.0)])
    np.testing.assert_allclose(found, sorted([m_plus, -m_minus]), rtol=2e-3)
    assert abs(m_minus) / m_plus > 4  # the Wagner term matters here


@pytest.mark.parametrize("etype", [EB, TS, REDDY])
def test_narrow_rectangle_moment_gradient(etype):
    """Michell cantilever (4.013) and the centrally loaded simple beam (16.94), TG 6.4-6.5."""
    L, n = 2.0, 32
    beam = Beam(etype, L, n, rect=(0.01, 0.2))
    r = beam.rect
    unit = np.sqrt(E * r["I_z"] * G * r["J"]) / L**2
    F = np.zeros(beam.n)
    F[beam.dof(n, 2)] = -1.0
    assert beam.buckle(F, beam.clamp_root())[0][0] / unit == pytest.approx(4.013, rel=2e-3)
    F = np.zeros(beam.n)
    F[beam.dof(n // 2, 2)] = -1.0
    assert beam.buckle(F, beam.fork())[0][0] / unit == pytest.approx(16.94, rel=2e-3)


# --------------------------------------------------------------------- #
# Anisotropic sections and element properties
# --------------------------------------------------------------------- #
def _unbalanced_channel(dy):
    lam = WallLaminate([Ply(CFRP, a, 0.5e-3) for a in (30, 30, 0, 30)])
    return channel(shift=(dy, 0.0), lam=lam)


def test_coupled_composite_buckling_independent_of_reference_axis():
    """Extension-twist and bend-twist coupled walls: the same loads at any origin."""
    L = 3.0
    an0 = _unbalanced_channel(0.0)
    yt, zt = an0.tension_centre
    out = []
    for dy in (0.0, yt, an0.shear_centre[0]):
        beam = Beam(BS, L, 16, sa=_unbalanced_channel(dy))
        out.append(beam.buckle(axial(beam, 1.0, y=yt - dy, z=zt), beam.fork(), n_modes=3)[0])
    np.testing.assert_allclose(out[1], out[0], rtol=1e-6)
    np.testing.assert_allclose(out[2], out[0], rtol=1e-6)


@pytest.mark.parametrize("etype", [EB, TS, BS, REDDY])
def test_k_sigma_symmetric_and_linear_in_state(etype):
    sa = channel() if etype in (TS, BS) else None
    beam = Beam(etype, 1.0, 2, sa=sa)
    rng = np.random.default_rng(3)
    U1, U2 = rng.normal(size=beam.n) * 1e-4, rng.normal(size=beam.n) * 1e-4
    K1, K2, K12 = beam.K_sigma(U1), beam.K_sigma(U2), beam.K_sigma(U1 + 2 * U2)
    np.testing.assert_allclose(K1, K1.T, atol=1e-12 * np.abs(K1).max())
    np.testing.assert_allclose(K12, K1 + 2 * K2, atol=1e-9 * np.abs(K12).max())


def test_gbt_has_no_geometric_stiffness_yet():
    from pre_processing.element_library.linear.beam.first_order_shear_deformation_theory.benscoter.linear_gbt_3D import (
        LinearGBTBeamElement3D,
    )

    with pytest.raises(NotImplementedError, match="Benscoter"):
        LinearGBTBeamElement3D.linear_geometric_stiffness_matrix(None, np.zeros(1))


def test_processing_buckling_path_matches_direct_eigensolve():
    from scipy.sparse import csr_matrix

    from processing.buckling.linear_buckling import (
        apply_buckling_boundary_conditions,
        assemble_global_geometric_stiffness,
        solve_linear_buckling_eigenpairs,
    )

    p = channel_properties()
    beam = Beam(BS, 4.0, 16, sa=channel())
    F = axial(beam, 1.0, y=p["y_c"])
    fixed = beam.fork()
    lam_direct, _ = beam.buckle(F, fixed, n_modes=2)
    U = beam.static(F, fixed)
    Kg = assemble_global_geometric_stiffness(beam.elements, U, beam.n)
    K_mod, Kg_mod, bc = apply_buckling_boundary_conditions(csr_matrix(beam.K), Kg, fixed_dofs=fixed)
    lam, _ = solve_linear_buckling_eigenpairs(K_mod, Kg_mod, 2, constrained_dofs=bc)
    np.testing.assert_allclose(lam, lam_direct, rtol=1e-8)


# --------------------------------------------------------------------- #
# P-delta
# --------------------------------------------------------------------- #
@pytest.mark.parametrize("etype", [EB, TS])
def test_p_delta_beam_column_amplification(etype):
    """Pinned beam-column, central lateral load: delta = Q L^3 / 48 EI * 3 (tan u - u) / u^3."""
    from processing.buckling.p_delta import solve_p_delta

    L, n = 3.0, 32
    beam = Beam(etype, L, n, rect=(0.06, 0.1))
    EI = E * beam.rect["I_z"]
    P = 0.6 * np.pi**2 * EI / L**2
    Q = 100.0
    F = axial(beam, P)
    F[beam.dof(n // 2, 1)] = Q
    res = solve_p_delta(beam.K, F, beam.elements, beam.fork())
    assert res.converged and res.iterations <= 3
    u = L / 2 * np.sqrt(P / EI)
    ref = Q * L**3 / (48 * EI) * 3 * (np.tan(u) - u) / u**3
    shear = Q * L / (4 * 5 / 6 * G * beam.rect["A"]) / (1 - P / (np.pi**2 * EI / L**2)) if etype == TS else 0.0
    assert res.U[beam.dof(n // 2, 1)] == pytest.approx(ref + shear, rel=2e-3)


def test_p_delta_tension_centre_offset_eccentric_column():
    """Compression at the web of a channel is eccentric to its tension centre: delta = e (sec(kL/2) - 1)."""
    from processing.buckling.p_delta import solve_p_delta

    L, n = 6.0, 32   # long enough that flexure along y comes before flexural-torsional buckling
    an = channel()
    beam = Beam(TS, L, n, sa=an)
    EIz = 1.0 / an.flexibility[2, 2]                  # about the tension centre
    P = 0.5 * np.pi**2 * EIz / L**2
    res = solve_p_delta(beam.K, axial(beam, P), beam.elements, beam.fork())
    assert res.converged
    e = -an.tension_centre[0]                         # load at the web, y = 0
    k = np.sqrt(P / EIz)
    v_rel = res.U[beam.dof(n // 2, 1)] - res.U[beam.dof(0, 1)]
    # the centroid line deflects by e (sec(kL/2) - 1) towards the load
    assert v_rel == pytest.approx(-e * (1 / np.cos(k * L / 2) - 1), rel=5e-3)
    lin = beam.static(axial(beam, P), beam.fork())[beam.dof(n // 2, 1)]
    assert abs(v_rel) > 1.5 * abs(lin)


def test_p_delta_past_buckling_is_an_unstable_state():
    from processing.buckling.p_delta import solve_p_delta

    beam = Beam(EB, 3.0, 16, rect=(0.06, 0.1))
    P = 1.3 * np.pi**2 * E * beam.rect["I_z"] / 3.0**2
    F = axial(beam, P)
    F[beam.dof(8, 1)] = 1.0
    res = solve_p_delta(beam.K, F, beam.elements, beam.fork(), max_iterations=10)
    # one fixed point exists but it is the unstable branch: (K + K_sigma) is indefinite
    Kt = beam.K + beam.K_sigma(res.U)
    free = np.setdiff1d(np.arange(beam.n), beam.fork())
    assert np.linalg.eigvalsh(0.5 * (Kt + Kt.T)[np.ix_(free, free)]).min() < 0


# --------------------------------------------------------------------- #
# Example (runs the buckling and static job runners end to end)
# --------------------------------------------------------------------- #
def test_channel_buckling_example(tmp_path):
    spec = importlib.util.spec_from_file_location(
        "channel_buckling_run",
        os.path.join(os.path.dirname(__file__), "..", "examples", "thin_walled", "channel_buckling", "run.py"),
    )
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    s = mod.main(n_el=16, plot=False, work_dir=str(tmp_path))
    assert abs(s["buckling"]["P_ft_error"]) < 5e-3
    assert s["p_delta"]["max_error"] < 1e-2
    for row in s["p_delta"]["rows"]:
        assert row["p_delta"] > row["linear"]

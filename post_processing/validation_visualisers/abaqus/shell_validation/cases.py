# post_processing/validation_visualisers/abaqus/shell_validation/cases.py
"""The Abaqus shell validation cases for the second-order and corotational beams.

Each case is a shell model (``ShellBeamModel``) and the 1d-fem beam job that should
reproduce it. Both are built from the same section description, so the walls, plies,
supports and loads match by construction.

==================  =========================  ===========================================
case                Abaqus analysis            1d-fem beam
==================  =========================  ===========================================
channel_ft          *BUCKLE, fork ends,        LinearBenscoterBeamElement3D buckling
                    compression at the
                    tension centre
ibeam_ltb_buckle    *BUCKLE, fork ends,        LinearBenscoterBeamElement3D buckling
                    uniform major-axis moment
ibeam_ltb_path      NLGEOM static, the same    CorotationalBenscoterBeamElement3D,
                    moment plus a small        one nonlinear job per load level
                    lateral force on the
                    compression flange
ibeam_large_twist   NLGEOM static, clamped     CorotationalBenscoterBeamElement3D
                    root, tip torque           (warping fixed at the root)
composite_box       NLGEOM static, clamped     CorotationalTimoshenkoBeamElement3D
                    root, tip flapwise force,
                    bend-twist coupled CFRP
==================  =========================  ===========================================

Public API
----------
CASES, build_case, steel, cfrp, channel_section, i_section, box_section.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Dict

import numpy as np

from post_processing.validation_visualisers.abaqus.shell_validation.shell_model import (
    EndCondition,
    Laminate,
    Material,
    NodalLoad,
    Ply,
    ShellBeamModel,
    ShellSection,
    Step,
)

E_STEEL, NU_STEEL = 2.0e11, 0.3
G_STEEL = E_STEEL / (2 * (1 + NU_STEEL))
STEEL = Material("STEEL", E_STEEL, E_STEEL, G_STEEL, NU_STEEL, G_STEEL, G_STEEL, 7850.0)
CFRP = Material("CFRP", 135e9, 10e9, 5e9, 0.3, 5e9, 3.8e9, 1550.0)

# channel of the channel_buckling example: web on y = 0, flanges towards +y
CH_H, CH_B, CH_T = 0.2, 0.08, 4e-3
# I-beam of the large_twist example
I_H, I_B, I_TF, I_TW = 0.3, 0.15, 6e-3, 4e-3
# bend-twist coupled CFRP box: 0.1 m square, 4 mm walls
BOX_W, BOX_H, T_PLY = 0.1, 0.1, 0.5e-3


def channel_section(element_size: float = 0.02) -> ShellSection:
    return ShellSection(
        points={1: (CH_B, -CH_H / 2), 2: (0.0, -CH_H / 2), 3: (0.0, CH_H / 2), 4: (CH_B, CH_H / 2)},
        walls=[(1, 1, 2, 1), (2, 2, 3, 1), (3, 3, 4, 1)],
        laminates={1: Laminate((Ply(STEEL, 0.0, CH_T),))},
        element_size=element_size,
    )


def i_section(element_size: float = 0.025) -> ShellSection:
    return ShellSection(
        points={1: (-I_B / 2, I_H / 2), 2: (0.0, I_H / 2), 3: (I_B / 2, I_H / 2),
                4: (-I_B / 2, -I_H / 2), 5: (0.0, -I_H / 2), 6: (I_B / 2, -I_H / 2)},
        walls=[(1, 1, 2, 1), (2, 2, 3, 1), (3, 4, 5, 1), (4, 5, 6, 1), (5, 5, 2, 2)],
        laminates={1: Laminate((Ply(STEEL, 0.0, I_TF),)), 2: Laminate((Ply(STEEL, 0.0, I_TW),))},
        element_size=element_size,
    )


def box_section(element_size: float = 0.0125) -> ShellSection:
    """Circumferentially asymmetric stiffness (CAS): the caps' off-axis plies point the same
    way in global axes (their wall directions are opposite), so flapwise bending twists the box."""
    def plies(angles):
        return Laminate(tuple(Ply(CFRP, a, T_PLY) for a in angles))

    bottom = plies((45, -45, 20, 20, 20, 20, -45, 45))
    top = plies((45, -45, -20, -20, -20, -20, -45, 45))
    web = plies((45, -45, 45, -45, -45, 45, -45, 45))
    return ShellSection(
        points={1: (-BOX_W / 2, -BOX_H / 2), 2: (BOX_W / 2, -BOX_H / 2), 3: (BOX_W / 2, BOX_H / 2),
                4: (-BOX_W / 2, BOX_H / 2)},
        walls=[(1, 1, 2, 1), (2, 2, 3, 3), (3, 3, 4, 2), (4, 4, 1, 3)],
        laminates={1: bottom, 2: top, 3: web},
        element_size=element_size,
    )


@dataclass
class Case:
    name: str
    model: ShellBeamModel
    beam: Dict                                      # 1d-fem job description (see onedfem.py)
    target: float                                   # relative tolerance of the comparison
    quantity: str                                   # what is compared


# -------------------------------------------------------------------- cases
def _channel_ft(fine: float = 1.0) -> Case:
    L = 4.0
    sec = channel_section(0.02 / fine)
    # compression at the centroid (= tension centre of a single-material section)
    from post_processing.validation_visualisers.abaqus.shell_validation.onedfem import tension_centre

    y_t, z_t = tension_centre(sec)
    model = ShellBeamModel(
        name="channel_ft", section=sec, length=L, n_x=int(round(80 * fine)),
        root=EndCondition("fork", (0.0, 0.0), fix_axial=True), tip=EndCondition("fork", (0.0, 0.0)),
        loads=[NodalLoad("root", (1.0, 0, 0, 0, z_t, -y_t)), NodalLoad("tip", (-1.0, 0, 0, 0, -z_t, y_t))],
        step=Step("buckle", n_modes=12),
        stations=tuple(np.linspace(0, 1, 9)),
        description="channel column, fork ends, unit compression at the tension centre",
    )
    beam = dict(kind="buckling", element="LinearBenscoterBeamElement3D", length=L, n_el=32, warping=1,
                supports="fork", loads=[(0.0, [1.0, 0, 0, 0, z_t, -y_t]), (L, [-1.0, 0, 0, 0, -z_t, y_t])])
    return Case("channel_ft", model, beam, 0.03, "lowest flexural-torsional buckling load")


def _ltb_loads(L: float, M: float, perturb: float):
    """End moments ``M`` about +y (bottom flange in compression) and a lateral perturbation
    on the compression flange at midspan."""
    loads = [NodalLoad("root", (0, 0, 0, 0, -M, 0)), NodalLoad("tip", (0, 0, 0, 0, M, 0))]
    beam = [(0.0, [0, 0, 0, 0, -M, 0]), (L, [0, 0, 0, 0, M, 0])]
    if perturb:
        Fy = perturb * 4.0 * M / L
        loads.append(NodalLoad(L / 2, (0, Fy, 0, 0, 0, 0), at=-I_H / 2))
        beam.append((L / 2, [0, Fy, 0, (I_H / 2) * Fy, 0, 0]))    # M_x = y F_z - z F_y at the bottom flange
    return loads, beam


IBEAM_LTB_L = 6.0


def _ibeam_ltb_buckle(fine: float = 1.0) -> Case:
    L = IBEAM_LTB_L
    loads, beam_loads = _ltb_loads(L, 1.0, 0.0)
    model = ShellBeamModel(
        name="ibeam_ltb_buckle", section=i_section(0.025 / fine), length=L, n_x=int(round(120 * fine)),
        root=EndCondition("fork", fix_axial=True), tip=EndCondition("fork"), loads=loads,
        step=Step("buckle", n_modes=10), stations=tuple(np.linspace(0, 1, 9)),
        description="I-beam, fork ends, unit uniform major-axis moment",
    )
    beam = dict(kind="buckling", element="LinearBenscoterBeamElement3D", length=L, n_el=32, warping=1,
                supports="fork", loads=beam_loads)
    return Case("ibeam_ltb_buckle", model, beam, 0.03, "lateral-torsional buckling moment")


def ltb_closed_form() -> float:
    """Timoshenko and Gere's uniform-moment LTB moment of the I-beam (mid-line properties)."""
    I_z = 2 * I_TF * I_B**3 / 12
    J = (2 * I_B * I_TF**3 + I_H * I_TW**3) / 3
    Gamma = I_TF * I_B**3 * I_H**2 / 24
    L = IBEAM_LTB_L
    return float(np.pi / L * np.sqrt(E_STEEL * I_z * G_STEEL * J * (1 + np.pi**2 * E_STEEL * Gamma / (G_STEEL * J * L**2))))


# Below the buckling moment: past it, separate load-controlled runs can land on different
# branches. Near it the perturbation's response is amplified by about 1 / (1 - M / M_cr).
LTB_LEVELS = (0.2, 0.4, 0.6, 0.7, 0.8, 0.85, 0.9, 0.95)
LTB_PERTURBATION = 1e-2


def _ibeam_ltb_path(fine: float = 1.0) -> Case:
    L = IBEAM_LTB_L
    M = 0.95 * ltb_closed_form()                      # load factor 1 = 0.95 M_cr0
    loads, beam_loads = _ltb_loads(L, M, LTB_PERTURBATION)
    levels = tuple(lv / 0.95 for lv in LTB_LEVELS)
    model = ShellBeamModel(
        name="ibeam_ltb_path", section=i_section(0.025 / fine), length=L, n_x=int(round(120 * fine)),
        root=EndCondition("fork", fix_axial=True), tip=EndCondition("fork"), loads=loads,
        step=Step("static", load_points=levels, max_increment=0.02), stations=tuple(np.linspace(0, 1, 9)),
        description=f"I-beam, fork ends, major-axis moment to 0.95 M_cr0 with a lateral force "
                    f"{LTB_PERTURBATION:g} * 4M/L on the compression flange at midspan",
    )
    beam = dict(kind="nonlinear", element="CorotationalBenscoterBeamElement3D", length=L, n_el=24, warping=1,
                supports="fork", loads=beam_loads, levels=levels, increments=10)
    return Case("ibeam_ltb_path", model, beam, 0.03, "midspan twist and lateral deflection along the path")


TWIST_TORQUE = 6000.0
TWIST_LEVELS = (1000 / 6000, 2000 / 6000, 4000 / 6000, 1.0)


def _ibeam_large_twist(fine: float = 1.0) -> Case:
    L = 2.0
    model = ShellBeamModel(
        name="ibeam_large_twist", section=i_section(0.025 / fine), length=L, n_x=int(round(80 * fine)),
        root=EndCondition("clamp"), tip=EndCondition("diaphragm"),
        loads=[NodalLoad("tip", (0, 0, 0, TWIST_TORQUE, 0, 0))],
        step=Step("static", load_points=TWIST_LEVELS, max_increment=0.05), stations=tuple(np.linspace(0, 1, 9)),
        description="I-beam cantilever, root clamped (warping restrained), tip torque to 6 kN m",
    )
    beam = dict(kind="nonlinear", element="CorotationalBenscoterBeamElement3D", length=L, n_el=16, warping=1,
                supports="clamp", loads=[(L, [0, 0, 0, TWIST_TORQUE, 0, 0])], levels=TWIST_LEVELS, increments=8)
    return Case("ibeam_large_twist", model, beam, 0.03, "tip twist")


BOX_L = 6.0
BOX_LEVELS = (0.25, 0.5, 0.75, 1.0)


def box_tip_force() -> float:
    """Tip force for a linear tip deflection of 20 % of the span."""
    from post_processing.validation_visualisers.abaqus.shell_validation.onedfem import section_analysis

    an = section_analysis(box_section())
    EI = 1.0 / an.flexibility[1, 1]                   # flapwise bending (curvature about y)
    return float(0.2 * BOX_L * 3 * EI / BOX_L**3)


def _composite_box(fine: float = 1.0) -> Case:
    L = BOX_L
    P = box_tip_force()
    model = ShellBeamModel(
        name="composite_box", section=box_section(0.0125 / fine), length=L, n_x=int(round(240 * fine)),
        root=EndCondition("clamp"), tip=EndCondition("rigid"),
        loads=[NodalLoad("tip", (0, 0, P, 0, 0, 0))],
        step=Step("static", load_points=BOX_LEVELS, max_increment=0.05), stations=tuple(np.linspace(0, 1, 9)),
        description=f"CAS CFRP box cantilever, tip flapwise force {P:.0f} N (20 % of span linear deflection)",
    )
    beam = dict(kind="nonlinear", element="CorotationalTimoshenkoBeamElement3D", length=L, n_el=24, warping=0,
                supports="clamp", loads=[(L, [0, 0, P, 0, 0, 0])], levels=BOX_LEVELS, increments=8)
    return Case("composite_box", model, beam, 0.03, "tip deflections and twist")


CASES: Dict[str, Callable[..., Case]] = {
    "channel_ft": _channel_ft,
    "ibeam_ltb_buckle": _ibeam_ltb_buckle,
    "ibeam_ltb_path": _ibeam_ltb_path,
    "ibeam_large_twist": _ibeam_large_twist,
    "composite_box": _composite_box,
}


def build_case(name: str, fine: float = 1.0) -> Case:
    """Case ``name``; ``fine`` scales the shell mesh density (1 = the validation mesh)."""
    return CASES[name](fine)


def steel() -> Material:
    return STEEL


def cfrp() -> Material:
    return CFRP


__all__ = ["CASES", "Case", "build_case", "steel", "cfrp", "channel_section", "i_section", "box_section",
           "ltb_closed_form", "box_tip_force"]

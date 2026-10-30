"""A 20.2 m outboard segment of the IEA 15 MW blade (r = 64.5 m to 84.7 m) as one IGK loft.

The 70 % span section idealisation of ``stations.py`` (station 20) is lofted with the
published chord taper (4.07 m to 3.27 m) and twist (1.47 deg to -0.36 deg) of the
r = 64.5 m and 84.7 m stations, the segment composite-laminate-model example 25 uses.
``prismatic=True`` keeps the r = 84.7 m section all along, untwisted.
"""

from __future__ import annotations

import stations as S

R0, CHORD0, TWIST0 = 64.48724137934501, 4.065767388212587, 1.4686902287805776   # r = 64.5 m (station 15)
TIP = S.station(20)
LENGTH = TIP["r"] - R0
MODEL = S.station_model(20)


def geometry(prismatic: bool = False):
    """The segment as an implicit-geometry-kernel ``BladeGeometry``, z from 0 (root) to ``LENGTH``."""
    from implicit_geometry_kernel import Airfoil, BladeGeometry, Station

    af = Airfoil(S.airfoil_points(20), "IEA 15 MW r = 84.7 m")
    if prismatic:
        sts = [Station(z, TIP["chord"], af, pitch_axis=MODEL.pitch_axis) for z in (0.0, LENGTH)]
    else:
        sts = [Station(0.0, CHORD0, af, twist=TWIST0, pitch_axis=MODEL.pitch_axis),
               Station(LENGTH, TIP["chord"], af, twist=TIP["twist_deg"], pitch_axis=MODEL.pitch_axis)]
    t = MODEL.thickness
    return BladeGeometry(sts, MODEL.layout, t["skin"], t["cap"], t["web"], web_positions=MODEL.web_positions)

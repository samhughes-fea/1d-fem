"""IEA 15 MW blade stations as implicit-geometry-kernel sections, and their 1d-fem models.

``data/iea15mw_stations.json`` holds three stations of the IEA Wind 15 MW reference blade
(root, maximum chord and 70 % span) as PreComp sees them: profile, sector layups, web
positions and PreComp's section properties (from composite-laminate-model example 24).

implicit-geometry-kernel (IGK) builds each station with one skin thickness per section.
The idealisation, used by every model compared here, is:

- **Skin:** the sandwich panel layup (gelcoat, triax, foam, triax) all round. The root
  has no foam: gelcoat and 100 mm of triax.
- **Spar caps:** the carbon layer of the cap sector, inside the skin between the webs.
- **Webs:** biax, foam, biax at the station's two web positions (none at the root).
- No trailing-edge reinforcement or trailing-edge web; a trailing-edge gap under 1.5 %
  of the chord is closed over the last 10 % of the chord.

All models run on this geometry, so the comparisons test the mechanics. PreComp, which
sees the full layup, is a secondary check.
"""

from __future__ import annotations

import json
import os
from dataclasses import dataclass
from typing import Dict, List, Tuple

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = json.load(open(os.path.join(HERE, "data", "iea15mw_stations.json"), encoding="utf-8"))
NU23 = 0.3          # transversely isotropic plies, as composite-laminate-model takes them
TE_CLOSE = 0.015    # trailing-edge gaps below this fraction of the chord are closed
NAMES = {0: "root", 6: "max chord", 20: "70 % span"}


@dataclass
class StationModel:
    index: int
    name: str
    chord: float
    pitch_axis: float
    layout: str
    web_positions: Tuple[float, ...]
    laminae: Dict[str, List[Tuple[str, float]]]    # kind -> (material, thickness), inside out

    @property
    def thickness(self) -> Dict[str, float]:
        return {k: sum(t for _, t in v) for k, v in self.laminae.items()}


def station(index: int) -> dict:
    return next(s for s in DATA["stations"] if s["index"] == index)


def _sector(st, side, has):
    """The first sector layup on ``side`` with a lamina of material ``has``."""
    for lam in st[side]["laminae"]:
        if any(m == has and t > 1e-3 for m, t, _ in lam):
            return lam
    raise KeyError(f"station {st['index']}: no {has} sector on the {side} surface")


def station_model(index: int) -> StationModel:
    st = station(index)
    inside_out = lambda lam: [(m, float(t)) for m, t, _ in reversed(lam) if t > 1e-4]  # noqa: E731
    if index == 0:
        skin = st["upper"]["laminae"][0]
        return StationModel(index, NAMES[index], st["chord"], st["pitch_axis"], "0A", (),
                            {"skin": inside_out(skin)})
    skin = _sector(st, "upper", "medium_density_foam")
    cap = [(m, t, a) for m, t, a in _sector(st, "upper", "CarbonUD") if m == "CarbonUD"]
    web = st["webs"]["laminae"][0]
    webs = tuple(float(x) for x in st["webs"]["loc"] if x < 0.99)
    return StationModel(index, NAMES[index], st["chord"], st["pitch_axis"], "2D-F", webs,
                        {"skin": inside_out(skin), "cap": inside_out(cap), "web": inside_out(web)})


def airfoil_points(index: int) -> np.ndarray:
    """Closed profile points (chord units), upper surface from the trailing edge, then lower."""
    p = station(index)["profile"]
    x = np.array(p["x"])
    yu, yl = np.array(p["y_upper"]), np.array(p["y_lower"])
    gap = yu[-1] - yl[-1]
    if 0.0 < abs(gap) < TE_CLOSE:
        ramp = np.clip((x - 0.9) / (x[-1] - 0.9), 0.0, 1.0)
        yu, yl = yu - 0.5 * gap * ramp, yl + 0.5 * gap * ramp
    return np.vstack([np.column_stack([x[::-1], yu[::-1]]), np.column_stack([x[1:], yl[1:]])])


def geometry(index: int, length: float = 1.0):
    """A prismatic IGK ``BladeGeometry`` of the station, ``length`` long (z from 0)."""
    from implicit_geometry_kernel import Airfoil, BladeGeometry, Station

    m = station_model(index)
    af = Airfoil(airfoil_points(index), f"IEA 15 MW station {index}")
    sts = [Station(z, m.chord, af, pitch_axis=m.pitch_axis) for z in (0.0, length)]
    t = m.thickness
    kw = {} if not m.web_positions else {"web_positions": m.web_positions}
    return BladeGeometry(sts, m.layout, t["skin"], t.get("cap", 0.0), t.get("web", 0.0), **kw)


def ply_materials():
    """1d-fem ``PlyMaterial`` per IEA material (transversely isotropic, ``nu23 = 0.3``)."""
    from pre_processing.section_integration.thin_walled import PlyMaterial

    out = {}
    for name, p in DATA["materials"].items():
        out[name] = PlyMaterial(E1=p["E1"], E2=p["E2"], G12=p["G12"], nu12=p["nu12"], G13=p["G12"],
                                G23=p["E2"] / (2.0 * (1.0 + NU23)), rho=p["rho"])
    return out


def laminates(index: int) -> Dict[str, object]:
    """1d-fem ``WallLaminate`` per wall kind."""
    from pre_processing.section_integration.thin_walled import Ply, WallLaminate

    mats = ply_materials()
    return {k: WallLaminate([Ply(mats[m], 0.0, t) for m, t in v]) for k, v in station_model(index).laminae.items()}

# onedfem/sections.py
"""The section provider: where 1d-fem gets its cross-section analyses.

1d-fem's beam elements need a section stiffness, and its pre- and post-processing need
wall laminates, section geometry, 2D section meshes, ply stress recovery and failure
criteria. None of that lives here. It comes from a **section provider**: a Python
package with the layout below. anisotropic-section-model is the default provider and
plugs in as it is; another package with the same modules and names can replace it::

    from onedfem.sections import set_section_provider
    set_section_provider("my_sections")        # a package name, or the module itself

The provider is imported the first time 1d-fem needs a section, so 1d-fem imports no
section package itself (``tests/test_repo_boundaries.py``). Install the default with
``pip install "onedfem[sections]"``.

Axes and order
--------------
A provider works in 1d-fem's beam axes: beam axis ``x``, section axes ``(y, z)``,
strains ``[eps, kappa_y, kappa_z, gamma_xy, gamma_xz, phi]`` and forces
``[N, M_y, M_z, V_y, V_z, T]``, right-handed (contract 2.0: the provider's *engine*
frame, ``SectionAnalysis.frame == "engine"``). ``docs/conventions.md`` gives the
reordering to BeamDyn.

Element signs
-------------
1d-fem's beam ``B`` matrices have ``kappa_y = -theta_y'``, so its elements read each
section with the sign of ``kappa_y`` flipped (``ELEMENT_SIGN``). That is 1d-fem's choice,
not the section's, so it is applied here and nowhere else: ``to_element_section`` re-signs
every signed field of a provider's section once, when an element takes it, and marks the
copy ``frame == "element"``. Station modifiers of the provider (``with_pretwist``,
``with_taper``) work on engine-frame sections: ``section_job`` applies tapers and aligns
GBT bases before any element sees a section, and ``element_pretwist`` re-signs back and
forth around ``with_pretwist`` (the flip is its own inverse, so this is exact).

Contract version
----------------
``SECTION_API`` is the version of the section contract this 1d-fem needs, ``(major,
minor)``. A provider declares the version it meets as ``SECTION_API_VERSION`` (package
level; anisotropic-section-model from 0.5.0). The first time a provider is used its major
number must match and its minor must be at least ``SECTION_API[1]``, or
``SectionProviderError`` is raised; a provider that declares no version is taken as 1.0
with a warning. The contract itself, and the split between the two packages, are in
anisotropic-section-model's ``docs/package_boundary.md``.

What a provider must have
-------------------------
``PROVIDER_MODULES`` lists the modules and the names 1d-fem takes from each. A
provider's objects are used as anisotropic-section-model documents them (for example a
thin-walled analysis has ``stiffness``, ``shear_centre``, ``wall_strains_at`` ...).
``check_section_provider`` reports anything missing.
"""

from __future__ import annotations

import dataclasses
import importlib
import warnings
import weakref
from types import ModuleType
from typing import Dict, List, Optional, Tuple, Union

import numpy as np

#: the package used when no other provider is set
DEFAULT_SECTION_PROVIDER = "anisotropic_section_model"

#: the section-provider contract this 1d-fem needs: (major, minor)
SECTION_API = (2, 1)

#: 1d-fem's element sign on the right-handed strains (``kappa_y = -theta_y'``)
ELEMENT_SIGN = np.diag([1.0, -1.0, 1.0, 1.0, 1.0, 1.0])
#: ``SectionAnalysis.frame`` of a provider's section, and of 1d-fem's re-signed copy
ENGINE_FRAME, ELEMENT_FRAME = "engine", "element"

#: module (relative to the provider) -> names 1d-fem uses from it
PROVIDER_MODULES: Dict[str, Tuple[str, ...]] = {
    "section_analysis": ("SectionAnalysis",),
    "axes": ("MIRROR", "to_engine"),
    "ply_failure": ("CRITERIA", "PlyStrength", "exposures", "ply_stress_from_wall"),
    "section_strength": ("SectionStrength",),
    "gbt_alignment": ("align_gbt_pair", "gbt_signature", "mass_transform"),
    "thin_walled.laminate": ("DEFAULT_SHEAR_CORRECTION", "Ply", "PlyMaterial", "WallLaminate"),
    "thin_walled.thin_walled_geometry": ("ThinWalledSection", "Wall"),
    "thin_walled.section_analysis": ("ThinWalledSectionAnalysis", "T_RATE"),
    "thin_walled.cross_section_modes": ("CrossSectionModes",),
    "thin_walled.stress_recovery": ("KeptField", "SaintVenantRecovery", "WallPointStress", "wall_point_stress"),
    "thin_walled.taper": ("SectionTaper", "section_taper_from_geometry"),
    "thin_walled.from_geometry": ("thin_walled_from_cut",),
    "thin_walled.layup_sensitivity": ("LayupParameter", "gbt_section_derivative", "section_analysis_derivative",
                                    "section_derivative"),
    "solid.material": ("SolidMaterial",),
    "solid.mesh": ("parse_section_mesh_2d",),
    "solid.elements": ("quadrature", "shape"),
    "solid.section_analysis": ("SolidSectionAnalysis",),
    "solid.from_geometry": ("section_mesh_from_cut", "write_section_mesh_2d"),
    "solid.layup_sensitivity": ("solid_section_derivative",),
    "solid.taper": ("solid_section_taper",),
}

_provider: Optional[Union[str, ModuleType]] = None
_checked: "weakref.WeakSet[ModuleType]" = weakref.WeakSet()


class SectionProviderError(ImportError):
    """No section provider is installed, or it lacks a module 1d-fem needs."""


def set_section_provider(provider: Union[str, ModuleType, None]) -> None:
    """Use ``provider`` (a package name or module) for sections; ``None`` restores the default."""
    global _provider
    _provider = provider


def check_api_version(module: ModuleType) -> Tuple[int, int]:
    """The contract version ``module`` declares, checked against ``SECTION_API``.

    Raises ``SectionProviderError`` when the major numbers differ or the provider's minor is
    older than the one needed; warns (once) and returns ``(1, 0)`` when it declares none.
    """
    v = getattr(module, "SECTION_API_VERSION", None)
    if v is None:
        if module not in _checked:
            warnings.warn(f"section provider '{module.__name__}' declares no SECTION_API_VERSION; "
                          f"taking it as contract 1.0", stacklevel=3)
        _checked.add(module)
        return (1, 0)
    major, minor = int(v[0]), int(v[1])
    if major != SECTION_API[0] or minor < SECTION_API[1]:
        raise SectionProviderError(
            f"section provider '{module.__name__}' meets section contract {major}.{minor}; this 1d-fem needs "
            f"{SECTION_API[0]}.x with x >= {SECTION_API[1]} (see onedfem.sections, SECTION_API)")
    _checked.add(module)
    return (major, minor)


def section_provider() -> ModuleType:
    """The section provider package, imported on first use and checked against ``SECTION_API``."""
    global _provider
    p = DEFAULT_SECTION_PROVIDER if _provider is None else _provider
    if isinstance(p, ModuleType):
        if p not in _checked:
            check_api_version(p)
        return p
    try:
        module = importlib.import_module(p)
    except ImportError as exc:
        hint = ' (pip install "onedfem[sections]")' if p == DEFAULT_SECTION_PROVIDER else ""
        raise SectionProviderError(f"1d-fem needs a section provider for cross-sections; '{p}' could not be "
                                   f"imported{hint}: {exc}") from exc
    check_api_version(module)
    _provider = module
    return module


def section_module(name: str) -> ModuleType:
    """The provider's module ``name`` (for example ``"thin_walled.laminate"``)."""
    root = section_provider()
    try:
        return importlib.import_module(f"{root.__name__}.{name}")
    except ImportError as exc:
        raise SectionProviderError(f"section provider '{root.__name__}' has no module '{name}': {exc}") from exc


def section_import(name: str, *names: str):
    """``from <provider>.<name> import <names>``: the objects, in order (one object for one name)."""
    module = section_module(name)
    out = []
    for n in names:
        try:
            out.append(getattr(module, n))
        except AttributeError as exc:
            raise SectionProviderError(f"section provider module '{module.__name__}' has no '{n}'") from exc
    return out[0] if len(out) == 1 else tuple(out)


def check_section_provider(provider: Union[str, ModuleType, None] = None) -> List[str]:
    """Modules and names of ``PROVIDER_MODULES`` that ``provider`` (default: the current one) lacks."""
    previous = _provider
    if provider is not None:
        set_section_provider(provider)
    missing: List[str] = []
    try:
        for name, names in PROVIDER_MODULES.items():
            try:
                module = section_module(name)
            except SectionProviderError:
                missing.append(name)
                continue
            missing += [f"{name}.{n}" for n in names if not hasattr(module, n)]
    finally:
        if provider is not None:
            set_section_provider(previous)
    return missing


# ---------------------------------------------------------------------- element signs
def _resign(sa, frame: str):
    """``sa`` with the sign of ``kappa_y`` flipped in every signed field, marked ``frame``."""
    s6 = np.diag(ELEMENT_SIGN).copy()
    s4 = s6[[0, 1, 2, 5]]
    s8 = np.r_[s6, 1.0, 1.0]

    def sym(A, s):
        return None if A is None else s[:, None] * np.asarray(A, dtype=np.float64) * s[None, :]

    kw = {
        "stiffness": sym(sa.stiffness, s6),
        "classical_stiffness": sym(sa.classical_stiffness, s4),
        "frame": frame,
    }
    for name in ("benscoter_stiffness", "third_order_stiffness"):
        kw[name] = sym(getattr(sa, name, None), s8)
    for name in ("third_order_rate", "geometric_terms"):          # columns: the six strains or forces
        A = getattr(sa, name, None)
        kw[name] = None if A is None else np.asarray(A, dtype=np.float64) * s6[None, :]
    Kg = getattr(sa, "gbt_stiffness", None)
    if Kg is not None:
        sg = np.ones(Kg.shape[0])
        sg[1] = -1.0
        kw["gbt_stiffness"] = sym(Kg, sg)
        G = getattr(sa, "gbt_geometric", None)
        kw["gbt_geometric"] = None if G is None else tuple(np.asarray(T, dtype=np.float64) * sg for T in G)
        Bz = getattr(sa, "gbt_brazier", None)
        if Bz is not None:
            sk = np.array([-1.0, 1.0])       # signs of (kappa_y, kappa_z)
            Tb, Qb = Bz
            kw["gbt_brazier"] = (np.asarray(Tb) * sk[None, :, None] * sg[None, None, :],
                                 np.asarray(Qb) * np.outer(sk, sk)[None, None])
    return dataclasses.replace(sa, **kw)


_element_copies: "Dict[int, Tuple[weakref.ref, object]]" = {}


def to_element_section(sa):
    """The provider's section ``sa`` as 1d-fem's elements read it (``ELEMENT_SIGN`` applied).

    A section already in the element frame is returned as it is; ``None`` stays ``None``.
    One copy is made per provider section, so elements sharing a section share its copy.
    """
    if sa is None or getattr(sa, "frame", ENGINE_FRAME) == ELEMENT_FRAME:
        return sa
    key = id(sa)
    hit = _element_copies.get(key)
    if hit is not None and hit[0]() is sa:
        return hit[1]
    out = _resign(sa, ELEMENT_FRAME)
    try:
        ref = weakref.ref(sa, lambda _r, k=key: _element_copies.pop(k, None))
    except TypeError:                     # not weak-referenceable: no sharing
        return out
    _element_copies[key] = (ref, out)
    return out


def to_engine_section(sa):
    """Inverse of ``to_element_section``: an element's section back in the provider's frame."""
    if sa is None or getattr(sa, "frame", ENGINE_FRAME) != ELEMENT_FRAME:
        return sa
    return _resign(sa, ENGINE_FRAME)


def element_pretwist(sa, rate: float):
    """``with_pretwist(rate)`` of an element-frame section, in the element frame."""
    if not rate:
        return sa
    if getattr(sa, "frame", ENGINE_FRAME) != ELEMENT_FRAME:
        return sa.with_pretwist(rate)
    return _resign(to_engine_section(sa).with_pretwist(rate), ELEMENT_FRAME)


__all__ = [
    "DEFAULT_SECTION_PROVIDER",
    "ELEMENT_FRAME",
    "ELEMENT_SIGN",
    "ENGINE_FRAME",
    "PROVIDER_MODULES",
    "SECTION_API",
    "SectionProviderError",
    "check_api_version",
    "check_section_provider",
    "element_pretwist",
    "section_import",
    "section_module",
    "section_provider",
    "set_section_provider",
    "to_element_section",
    "to_engine_section",
]

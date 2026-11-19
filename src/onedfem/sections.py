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
``[N, M_y, M_z, V_y, V_z, T]``. ``docs/conventions.md`` gives the reordering to BeamDyn.

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

import importlib
import warnings
import weakref
from types import ModuleType
from typing import Dict, List, Optional, Tuple, Union

#: the package used when no other provider is set
DEFAULT_SECTION_PROVIDER = "anisotropic_section_model"

#: the section-provider contract this 1d-fem needs: (major, minor)
SECTION_API = (1, 0)

#: module (relative to the provider) -> names 1d-fem uses from it
PROVIDER_MODULES: Dict[str, Tuple[str, ...]] = {
    "section_analysis": ("SectionAnalysis", "ELEMENT_SIGN"),
    "axes": ("MIRROR", "to_engine"),
    "ply_failure": ("CRITERIA", "PlyStrength", "exposures", "ply_stress_from_wall"),
    "section_strength": ("SectionStrength",),
    "gbt_alignment": ("align_gbt_sections", "mass_transform"),
    "thin_walled.laminate": ("DEFAULT_SHEAR_CORRECTION", "Ply", "PlyMaterial", "WallLaminate"),
    "thin_walled.thin_walled_geometry": ("ThinWalledSection", "Wall"),
    "thin_walled.section_analysis": ("ThinWalledSectionAnalysis", "T_RATE"),
    "thin_walled.cross_section_modes": ("CrossSectionModes",),
    "thin_walled.stress_recovery": ("KeptField", "SaintVenantRecovery", "WallPointStress", "wall_point_stress"),
    "thin_walled.taper": ("SectionTaper", "section_taper_from_geometry"),
    "thin_walled.from_geometry": ("thin_walled_from_cut",),
    "thin_walled.layup_sensitivity": ("LayupParameter", "element_section_derivative", "gbt_section_derivative",
                                    "section_derivative"),
    "solid.material": ("SolidMaterial",),
    "solid.mesh": ("parse_section_mesh_2d",),
    "solid.elements": ("quadrature", "shape"),
    "solid.section_analysis": ("SolidSectionAnalysis",),
    "solid.from_geometry": ("section_mesh_from_cut", "write_section_mesh_2d"),
    "solid.layup_sensitivity": ("solid_section_derivative",),
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


__all__ = [
    "DEFAULT_SECTION_PROVIDER",
    "PROVIDER_MODULES",
    "SECTION_API",
    "SectionProviderError",
    "check_api_version",
    "check_section_provider",
    "section_import",
    "section_module",
    "section_provider",
    "set_section_provider",
]

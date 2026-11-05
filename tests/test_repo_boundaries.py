"""1d-fem's library code imports none of its sibling repositories.

implicit-geometry-kernel, composite-laminate-model and anisotropic-section-model may
appear in examples, validation and tests (as optional references), never in
``src/onedfem``. Blade geometry is accepted duck typed.
"""

import ast
from pathlib import Path

SIBLINGS = ("implicit_geometry_kernel", "composite_laminate_model", "anisotropic_section_model",
            "laminate_fem", "laminate_modelling_package")
SRC = Path(__file__).resolve().parents[1] / "src" / "onedfem"


def test_library_imports_no_sibling_repository():
    offending = []
    for path in SRC.rglob("*.py"):
        for node in ast.walk(ast.parse(path.read_text(encoding="utf-8-sig"))):
            if isinstance(node, ast.Import):
                names = [a.name for a in node.names]
            elif isinstance(node, ast.ImportFrom) and node.module and node.level == 0:
                names = [node.module]
            else:
                continue
            for name in names:
                if name.split(".")[0] in SIBLINGS:
                    offending.append(f"{path.relative_to(SRC)}:{node.lineno} {name}")
    assert not offending, "library code imports a sibling repository:\n" + "\n".join(offending)


def test_blade_adapters_live_in_interfaces():
    """Blade adapters sit in ``onedfem.interfaces.blade``; the old module names still import them."""
    import importlib

    moved = {"onedfem.pre_processing.mesh_library.blade_beam_job": "onedfem.interfaces.blade.blade_beam_job",
             "onedfem.pre_processing.section_integration.thin_walled.from_geometry":
                 "onedfem.interfaces.blade.thin_walled_sections",
             "onedfem.pre_processing.section_integration.solid.from_geometry":
                 "onedfem.interfaces.blade.solid_sections"}
    for old, new in moved.items():
        assert importlib.import_module(old) is importlib.import_module(new)
    core = [p for p in SRC.rglob("*.py") if "interfaces" not in p.relative_to(SRC).parts]
    users = [str(p.relative_to(SRC)) for p in core if "onedfem.interfaces" in p.read_text(encoding="utf-8-sig")
             and p.stem not in ("blade_beam_job", "from_geometry")]
    assert not users, f"core modules import the interfaces: {users}"

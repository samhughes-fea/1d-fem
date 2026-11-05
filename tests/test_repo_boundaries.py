"""1d-fem's library code imports none of its sibling repositories.

implicit-geometry-kernel, laminate-constitutive-model, shell-fem and
anisotropic-section-model may
appear in examples, validation and tests (as optional references), never in
``src/onedfem``. Blade geometry is accepted duck typed.
"""

import ast
from pathlib import Path

SIBLINGS = ("implicit_geometry_kernel", "laminate_constitutive_model", "shell_fem", "anisotropic_section_model",
            "composite_laminate_model", "laminate_fem", "laminate_modelling_package")
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

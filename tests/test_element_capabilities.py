"""Element capabilities are class attributes read through ``ElementFactory``, not name lists.

Each flag is checked against the registry, and the library outside the element library
holds no element type-name literals or name-fragment tests (``"GBT" in type``) that would
decide a capability by name.
"""

import ast
import re
from pathlib import Path

import pytest

from onedfem.pre_processing.element_library.element_factory import ElementFactory as F

SRC = Path(__file__).resolve().parents[1] / "src" / "onedfem"
LIBRARY = SRC / "pre_processing" / "element_library"
TYPE_NAME = re.compile(r"^(Linear|Nonlinear|Corotational|UpdatedLagrangian|GEBT)\w*(Element3D|Beam3D)$")
FRAGMENTS = ("GBT", "Benscoter", "Reddy", "Levinson", "Timoshenko", "EulerBernoulli", "Corotational", "Nonlinear")


def test_flags_on_the_registered_types():
    assert F.types_with("cross_section_modes") == ["CorotationalGBTBeamElement3D", "LinearGBTBeamElement3D"]
    assert set(F.types_with("restrained_warping")) >= {"LinearBenscoterBeamElement3D", "LinearGBTBeamElement3D"}
    assert F.capability("LinearReddyBeamElement3D", "extra_dofs") == 2
    assert F.capability("LinearTimoshenkoBeamElement3D", "extra_dofs") == 0
    assert F.capability("LinearEulerBernoulliBeamElement3D", "bending_interpolation") == "hermite"
    assert F.types_with("reduced_isotropic_section", own=True) == [
        "CorotationalEulerBernoulliBeamElement3D", "CorotationalTimoshenkoBeamElement3D",
        "LinearEulerBernoulliBeamElement3D", "LinearTimoshenkoBeamElement3D",
    ]
    assert F.capability("LinearTimoshenkoBeamElement3D", "nonlinear_counterpart", own=True) \
        == "NonlinearTimoshenkoBeamElement3D"
    assert F.capability("LinearGBTBeamElement3D", "layup_sensitivity", own=True) == "gbt"
    assert F.capability("UnknownElement3D", "extra_dofs", 7) == 7


@pytest.mark.parametrize("name", sorted(F.ELEMENT_CLASS_MAP))
def test_counterparts_are_registered(name):
    twin = F.capability(name, "nonlinear_counterpart", own=True)
    assert twin is None or twin in F.ELEMENT_CLASS_MAP


def test_no_type_names_outside_the_element_library():
    found = []
    for path in SRC.rglob("*.py"):
        if LIBRARY in path.parents:
            continue
        for node in ast.walk(ast.parse(path.read_text(encoding="utf-8-sig"))):
            if isinstance(node, ast.Compare) and any(isinstance(o, (ast.In, ast.NotIn)) for o in node.ops):
                if isinstance(node.left, ast.Constant) and node.left.value in FRAGMENTS:
                    found.append(f"{path.relative_to(SRC)}:{node.lineno} {ast.unparse(node)[:80]}")
            if isinstance(node, (ast.Tuple, ast.List, ast.Set)):
                names = [e.value for e in node.elts if isinstance(e, ast.Constant) and isinstance(e.value, str)
                         and TYPE_NAME.match(e.value)]
                if len(names) >= 2:
                    found.append(f"{path.relative_to(SRC)}:{node.lineno} list of type names {names[:3]}")
    assert not found, "capability decided by element name:\n" + "\n".join(found)

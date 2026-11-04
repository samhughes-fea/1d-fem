# pre_processing/parsing/layup_sensitivity_parser.py
"""Parse ``layup_sensitivity.txt``: the ply parameters and responses of a sensitivity run.

Format::

    [Layup_Parameter]
    [laminate_id] [ply] [kind]
    1 5 angle
    1 5 thickness
    [Response]
    [type] [node_id] [dof] [mode]
    displacement 24 UZ -
    displacement 24 RX -
    compliance - - -
    frequency - - 1
    buckling - - 1

- ``ply`` counts from 1 at the bottom of the laminate (as ``[ply]`` in ``laminate.txt``);
  ``kind`` is ``angle`` (derivative per degree) or ``thickness`` (per metre).
- Responses: ``displacement`` (a nodal DOF ``UX`` ... ``RZ`` or ``Qk``) and ``compliance``
  (``F^T U``) in static jobs; ``frequency`` (mode number from 1, in Hz) in eigen jobs;
  ``buckling`` (the load factor of buckling mode number from 1) in linear buckling jobs.

Public API
----------
LAYUP_SENSITIVITY_FILE, Response, parse_layup_sensitivity_file, format_layup_sensitivity_file.
"""

from __future__ import annotations

import os
from dataclasses import dataclass
from typing import List, Optional, Sequence, Tuple

from onedfem.pre_processing.section_integration.thin_walled.layup_sensitivity import LayupParameter

LAYUP_SENSITIVITY_FILE = "layup_sensitivity.txt"
RESPONSE_TYPES = ("displacement", "compliance", "frequency", "buckling")


@dataclass(frozen=True)
class Response:
    """A response of a sensitivity run (see the module docstring)."""

    type: str
    node_id: Optional[int] = None
    dof: Optional[str] = None
    mode: Optional[int] = None

    def __post_init__(self) -> None:
        if self.type not in RESPONSE_TYPES:
            raise ValueError(f"response type must be one of {RESPONSE_TYPES}, got {self.type!r}")
        if self.type == "displacement" and (self.node_id is None or self.dof is None):
            raise ValueError("a displacement response needs a node_id and a dof")
        if self.type in ("frequency", "buckling") and (self.mode is None or self.mode < 1):
            raise ValueError(f"a {self.type} response needs a mode number from 1")

    @property
    def label(self) -> str:
        if self.type == "displacement":
            return f"{self.dof}@{self.node_id}"
        if self.type == "frequency":
            return f"f{self.mode}"
        if self.type == "buckling":
            return f"lambda{self.mode}"
        return "compliance"


def _opt(v: str, cast):
    return None if v == "-" else cast(v)


def parse_layup_sensitivity_file(path: str) -> Tuple[List[LayupParameter], List[Response]]:
    """Parameters and responses from ``path``; two empty lists when the file is absent."""
    if not os.path.exists(path):
        return [], []
    with open(path, encoding="utf-8") as f:
        lines = [ln.split("#")[0].strip() for ln in f]
    lines = [ln for ln in lines if ln]
    params: List[LayupParameter] = []
    responses: List[Response] = []
    block = None
    for ln in lines:
        low = ln.lower()
        if low == "[layup_parameter]":
            block = "p"
            continue
        if low == "[response]":
            block = "r"
            continue
        if ln.startswith("["):
            continue                                     # column header
        tok = ln.split()
        if block == "p":
            if len(tok) != 3:
                raise ValueError(f"{path}: a parameter row is [laminate_id] [ply] [kind], got {ln!r}")
            ply = int(tok[1])
            if ply < 1:
                raise ValueError(f"{path}: plies count from 1, got {ply}")
            params.append(LayupParameter(int(tok[0]), ply - 1, tok[2].lower()))
        elif block == "r":
            if len(tok) != 4:
                raise ValueError(f"{path}: a response row is [type] [node_id] [dof] [mode], got {ln!r}")
            responses.append(Response(tok[0].lower(), _opt(tok[1], int), _opt(tok[2], str.upper), _opt(tok[3], int)))
        else:
            raise ValueError(f"{path}: rows must follow [Layup_Parameter] or [Response]")
    return params, responses


def format_layup_sensitivity_file(params: Sequence[LayupParameter], responses: Sequence[Response]) -> str:
    """The text of ``layup_sensitivity.txt`` for ``params`` and ``responses``."""
    out = ["[Layup_Parameter]", "[laminate_id] [ply] [kind]"]
    out += [f"{p.laminate_id} {p.ply + 1} {p.kind}" for p in params]
    out += ["[Response]", "[type] [node_id] [dof] [mode]"]
    for r in responses:
        f = lambda v: "-" if v is None else str(v)  # noqa: E731
        out.append(f"{r.type} {f(r.node_id)} {f(r.dof)} {f(r.mode)}")
    return "\n".join(out) + "\n"

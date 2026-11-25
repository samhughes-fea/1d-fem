# onedfem/pre_processing/section_job.py
"""Build section analyses from a job's ``laminate.txt`` and ``section_geometry.txt`` (thin-walled)
and ``section_mesh_2d.txt`` (sections meshed in 2D)."""

from __future__ import annotations

import logging
from typing import Dict

import numpy as np

from onedfem.pre_processing.parsing.laminate_parser import parse_laminate_file
from onedfem.pre_processing.parsing.section_geometry_parser import (
    parse_section_geometry_file,
    parse_section_mesh,
    parse_section_mesh_2d_block,
)
from onedfem.sections import section_import, section_module

log = logging.getLogger(__name__)


def section_analyses_from_files(
    laminate_path: str, geometry_path: str, hoop: str = "free", n_sub: int = 1
) -> Dict[int, ThinWalledSectionAnalysis]:
    """Return ``{section_id: ThinWalledSectionAnalysis}`` for every section in the files.

    ``n_sub`` is the default wall subdivision; a ``[Section_Mesh]`` block overrides it per section.
    Sections listed in a ``[Section_Mesh_2D]`` block are left out (``solid_analyses_from_geometry``).
    """
    ThinWalledSectionAnalysis = section_import("thin_walled.section_analysis", "ThinWalledSectionAnalysis")
    laminates = parse_laminate_file(laminate_path)
    sections = parse_section_geometry_file(geometry_path, laminates.keys())
    mesh = parse_section_mesh(geometry_path)
    two_d = parse_section_mesh_2d_block(geometry_path)
    return {
        sid: ThinWalledSectionAnalysis(sec, laminates, hoop=hoop, n_sub=mesh.get(sid, n_sub))
        for sid, sec in sections.items()
        if sid not in two_d
    }


def solid_analyses_from_geometry(laminate_path: str, geometry_path: str) -> Dict[int, "SolidSectionAnalysis"]:
    """Return ``{section_id: SolidSectionAnalysis}`` for the sections of the ``[Section_Mesh_2D]``
    block of ``section_geometry.txt`` (empty without it).

    Each is meshed in 2D from its walls and ``laminate.txt`` (``SectionMesh.from_thin_walled``),
    so its mesh carries ply tags and the node velocities of every ply thickness: such sections
    take ply-angle and ply-thickness sensitivities.
    """
    two_d = parse_section_mesh_2d_block(geometry_path)
    if not two_d:
        return {}
    SectionMesh = section_import("solid.mesh", "SectionMesh")
    SolidSectionAnalysis = section_import("solid.section_analysis", "SolidSectionAnalysis")
    laminates = parse_laminate_file(laminate_path)
    sections = parse_section_geometry_file(geometry_path, laminates.keys())
    unknown = sorted(set(two_d) - set(sections))
    if unknown:
        raise ValueError(f"[Section_Mesh_2D] names sections {unknown} that section_geometry.txt does not define")
    out = {}
    for sid, (along, rows_per_ply) in two_d.items():
        mesh = SectionMesh.from_thin_walled(sections[sid], laminates, elements_along=along, rows_per_ply=rows_per_ply)
        mesh.section_id = sid
        out[sid] = SolidSectionAnalysis(mesh)
    return out


def solid_tapers(rows_2d, tapers, solid) -> Dict[int, "SectionTaper"]:
    """``{element_id: SectionTaper}`` of the ``[Section_Taper_2D]`` rows of ``section_taper.txt``.

    ``rows_2d`` is ``{element_id: (section_start, section_end, distance)}``, ``tapers`` the
    ``[Section_Taper]`` rows (an element takes one or the other) and ``solid`` the job's 2D
    sections. Each coupling is the provider's ``solid_section_taper`` of the two meshes, made
    once per row value.
    """
    if not rows_2d:
        return {}
    both = sorted(set(rows_2d) & set(tapers))
    if both:
        raise ValueError(f"section_taper.txt gives elements {both} rows in both [Section_Taper] and [Section_Taper_2D]")
    named = {s for a, b, _ in rows_2d.values() for s in (a, b)}
    not_2d = sorted(named - set(solid))
    if not_2d:
        raise ValueError(f"[Section_Taper_2D] names sections {not_2d} that are not meshed in 2D ([Section_Mesh_2D] "
                         "in section_geometry.txt or section_mesh_2d.txt)")
    solid_section_taper = section_import("solid.taper", "solid_section_taper")
    made = {}
    out = {}
    for e, (a, b, d) in rows_2d.items():
        if (a, b, d) not in made:
            try:
                made[a, b, d] = solid_section_taper(solid[a].mesh, solid[b].mesh, d)
            except ValueError as exc:
                raise ValueError(f"[Section_Taper_2D] element {e}: sections {a} and {b}: {exc}") from exc
        out[e] = made[a, b, d]
    return out


def solid_analyses_from_files(mesh_path: str, laminate_path=None) -> Dict[int, "SolidSectionAnalysis"]:
    """Return ``{section_id: SolidSectionAnalysis}`` for every section in ``section_mesh_2d.txt``.

    Material ids not in its ``[Solid_Material]`` block come from the ``[Ply_Material]`` block of
    ``laminate.txt`` (``SolidMaterial.from_ply``).
    """
    from onedfem.pre_processing.parsing.laminate_parser import parse_ply_material_file
    SolidMaterial = section_import("solid.material", "SolidMaterial")
    parse_section_mesh_2d = section_import("solid.mesh", "parse_section_mesh_2d")
    SolidSectionAnalysis = section_import("solid.section_analysis", "SolidSectionAnalysis")

    mats = {}
    if laminate_path is not None:
        mats = {mid: SolidMaterial.from_ply(m) for mid, m in parse_ply_material_file(laminate_path).items()}
    return {sid: SolidSectionAnalysis(mesh) for sid, mesh in parse_section_mesh_2d(mesh_path, mats).items()}


LAMINATE_FILE = "laminate.txt"
SECTION_MESH_2D_FILE = "section_mesh_2d.txt"
SOLID_KEY = "solid_analysis_per_element"
#: each element's ``SectionTaper`` (``None`` where it has none): the layup sensitivities of a
#: tapered 2D section differentiate its coupling through the two meshes it keeps
TAPER_KEY = "section_taper_per_element"
SECTION_GEOMETRY_FILE = "section_geometry.txt"
ELEMENT_SECTION_FILE = "element_section.txt"
WALL_LOAD_FILE = "wall_load.txt"
SECTION_TAPER_FILE = "section_taper.txt"
SECTION_ANALYSIS_KEY = "section_analysis_per_element"
THIN_WALLED_KEY = "thin_walled_analysis_per_element"



def align_gbt_sections(per_element, connectivity, gbt_elements=None):
    """Turn the GBT bases of each element's section to its neighbour's.

    ``per_element`` holds each element's ``SectionAnalysis`` (or ``None``), in the provider's
    frame, ``connectivity`` its two nodes. Elements are visited in order; each is aligned
    (``gbt_alignment.align_gbt_pair``) to an earlier visited element sharing a node whose
    section differs, or shares the (turned) section of a neighbour that has the same one.
    Returns the aligned sections and, per element whose bases were turned,
    ``{"Rw", "Rm", "match", "reference"}``. The walk is 1d-fem's (contract 2.0); the
    provider only turns one section to another.
    """
    align_gbt_pair, gbt_signature = section_import("gbt_alignment", "align_gbt_pair", "gbt_signature")
    MIN_MATCH = getattr(section_module("gbt_alignment"), "MIN_MATCH", 0.9)

    conn = np.asarray(connectivity, dtype=np.int64).reshape(-1, 2)
    use = [sa is not None and getattr(sa, "gbt_stiffness", None) is not None for sa in per_element]
    if gbt_elements is not None:
        use = [u and bool(g) for u, g in zip(use, gbt_elements)]
    out = list(per_element)
    info = {}
    at_node = {}
    sig_cache = {}

    def blocks(sa):
        return sa.gbt_sv_modes + sa.gbt_shear_modes, sa.gbt_modes

    def sig(i):
        if i not in sig_cache:
            sig_cache[i] = gbt_signature(out[i])
        return sig_cache[i]

    for e in range(len(per_element)):
        if use[e]:
            nb = [f for nd in conn[e] for f in at_node.get(int(nd), []) if use[f]]
            same = next((f for f in nb if per_element[f] is per_element[e]), None)
            ref = next((f for f in nb if per_element[f] is not per_element[e]), None)
            if same is not None:                     # the same section: the same (turned) basis
                out[e] = out[same]
                if same in info:
                    info[e] = dict(info[same])
            elif ref is not None and blocks(out[ref]) == blocks(out[e]):
                turned, pair = align_gbt_pair(out[ref], out[e], (sig(ref), sig(e)))
                if pair["match"] < MIN_MATCH:
                    log.warning("GBT bases of element %d do not continue those of element %d", e, ref)
                if turned is not out[e]:
                    out[e] = turned
                    sig_cache.pop(e, None)
                info[e] = {**pair, "reference": ref}
        for nd in conn[e]:
            at_node.setdefault(int(nd), []).append(e)
    return tuple(out), info

def attach_section_analyses(job_dir: str, element_dictionary: dict, hoop: str = "free") -> dict:
    """Return ``element_dictionary`` with per-element ``SectionAnalysis`` entries when the job has thin-walled sections.

    Looks for ``laminate.txt``, ``section_geometry.txt`` and ``element_section.txt``
    in ``job_dir``. With none present the dictionary is returned unchanged. With
    some but not all present, raises ``FileNotFoundError``. Otherwise adds key
    ``"section_analysis_per_element"``: a tuple aligned with
    ``element_dictionary["ids"]`` holding a ``SectionAnalysis`` for each assigned
    element and ``None`` for the rest. Each section is analysed once. Key
    ``"thin_walled_analysis_per_element"`` holds the ``ThinWalledSectionAnalysis``
    behind each, for stress recovery.

    Rows of ``LinearGBTBeamElement3D`` (or rows with ``[modes]`` / ``[shear_modes]``
    columns) make their sections build GBT matrices with as many modes as any element
    asks for. The dictionary then gains ``"gbt_modes"``, ``"gbt_shear_modes"`` and
    ``"extra_dofs"`` (``n_sv + shear_modes + 2 modes`` for those rows, ``n_sv`` the
    section's Saint-Venant warping modes) for the mesh DOF layout.

    Thin-walled sections use ``hoop`` (default ``"free"``, which the 2D and shell
    models match), except sections that a GBT row uses: the cross-section modes carry
    the hoop bending, so those are analysed with ``hoop="restrained"``.

    With ``section_mesh_2d.txt`` the sections it defines are analysed in 2D
    (``SolidSectionAnalysis``, ``SectionAnalysis.from_solid``); they can be assigned in the
    same ``element_section.txt`` as thin-walled ones (section ids must differ). Key
    ``"solid_analysis_per_element"`` then holds the ``SolidSectionAnalysis`` behind each
    (``None`` in ``"thin_walled_analysis_per_element"``). GBT rows on 2D sections use the
    2D section's own modes (``solid.cross_section_modes``) without extra shear-lag modes;
    wall loads need thin-walled sections.

    Sections listed in a ``[Section_Mesh_2D]`` block of ``section_geometry.txt`` are meshed in 2D
    from their walls (``SectionMesh.from_thin_walled``) and analysed like those of
    ``section_mesh_2d.txt``; their meshes also take ply-thickness sensitivities.

    With ``section_matrix.txt`` the sections it defines are given as 6 x 6 stiffness and mass
    tables (``parsing/section_matrix_parser.py``, ``SectionAnalysis.from_matrices``) and are
    assigned in the same ``element_section.txt`` (section ids must differ from the thin-walled
    and 2D ones). They take no GBT modes, taper, wall loads or stress recovery.

    The section provider is imported only when the job has section files, so a job of
    ``section.txt`` sections alone starts without it.

    With ``section_taper.txt`` the listed elements' sections gain the taper (Resal) coupling
    (``SectionAnalysis.with_taper``), GBT rows included.

    GBT rows on more than one section have their sections' warping and mode bases turned
    to continue each other along the beam (``align_gbt_sections``, pairs by the provider's
    ``gbt_alignment.align_gbt_pair``), since
    neighbouring elements share those amplitudes.

    With ``wall_load.txt`` it also gains ``"wall_load_per_element"``: the generalised
    load per length of each loaded element (``None`` for the rest), in the order of the
    section's GBT mass for GBT rows, else ``[N, F_y, F_z, M_x, M_y, M_z, psi load]``.
    """
    import os

    import numpy as np

    from onedfem.pre_processing.parsing.element_section_parser import parse_element_section_modes
    from onedfem.pre_processing.parsing.section_matrix_parser import SECTION_MATRIX_FILE, section_matrix_analyses

    paths = [os.path.join(job_dir, f) for f in (LAMINATE_FILE, SECTION_GEOMETRY_FILE, ELEMENT_SECTION_FILE)]
    present = [os.path.exists(p) for p in paths]
    mesh_path = os.path.join(job_dir, SECTION_MESH_2D_FILE)
    has_mesh = os.path.exists(mesh_path)
    matrix_path = os.path.join(job_dir, SECTION_MATRIX_FILE)
    has_matrix = os.path.exists(matrix_path)
    if not any(present) and not has_mesh and not has_matrix:
        return element_dictionary
    # imported only now: a job without section files needs no section provider
    ThinWalledSectionAnalysis = section_import("thin_walled.section_analysis", "ThinWalledSectionAnalysis")
    SectionAnalysis = section_import("section_analysis", "SectionAnalysis")
    if has_mesh or has_matrix:
        # 2D sections: element_section.txt is needed; laminate.txt only for its ply materials
        # (or with section_geometry.txt for thin-walled sections alongside)
        need = [paths[2]] + ([paths[0]] if present[1] else [])
    else:
        need = paths
    missing = [os.path.basename(p) for p in need if not os.path.exists(p)]
    if missing:
        raise FileNotFoundError(
            "Thin-walled sections need laminate.txt, section_geometry.txt and element_section.txt; "
            "2D sections need section_mesh_2d.txt and element_section.txt, table sections section_matrix.txt and "
            f"element_section.txt; missing {missing}"
        )
    analyses = section_analyses_from_files(paths[0], paths[1], hoop=hoop) if present[1] else {}
    solid = solid_analyses_from_files(mesh_path, paths[0] if present[0] else None) if has_mesh else {}
    meshed = solid_analyses_from_geometry(paths[0], paths[1]) if present[1] else {}
    clash = sorted(set(meshed) & set(solid))
    if clash:
        raise ValueError(f"sections {clash} are meshed in 2D by both [Section_Mesh_2D] in section_geometry.txt "
                         "and section_mesh_2d.txt")
    solid.update(meshed)
    matrices = section_matrix_analyses(matrix_path) if has_matrix else {}
    both = sorted(set(analyses) & set(solid))
    if both:
        raise ValueError(f"sections {both} are defined in both section_geometry.txt and section_mesh_2d.txt")
    twice = sorted((set(analyses) | set(solid)) & set(matrices))
    if twice:
        raise ValueError(f"sections {twice} are defined in section_matrix.txt and in section_geometry.txt or "
                         "section_mesh_2d.txt")
    assignment_modes = parse_element_section_modes(paths[2])
    assignment = {e: v[0] for e, v in assignment_modes.items()}
    modes_of = {e: v[1] for e, v in assignment_modes.items()}
    shear_of = {e: v[2] for e, v in assignment_modes.items()}
    ids = [int(i) for i in np.asarray(element_dictionary["ids"]).ravel()]
    unknown_el = sorted(set(assignment) - set(ids))
    if unknown_el:
        raise ValueError(f"element_section.txt names unknown elements {unknown_el}")
    unknown_sec = sorted(set(assignment.values()) - set(analyses) - set(solid) - set(matrices))
    if unknown_sec:
        raise ValueError(f"element_section.txt names unknown sections {unknown_sec}")
    types = element_dictionary.get("types")
    type_of = {} if types is None else {e: str(t) for e, t in zip(ids, np.asarray(types).ravel())}
    from onedfem.pre_processing.element_library.element_factory import ElementFactory

    def has_modes(e):
        t = type_of.get(e, "")
        return bool(t) and bool(ElementFactory.capability(t, "cross_section_modes", False))

    gbt_rows = {e for e in assignment if has_modes(e) or modes_of[e] or shear_of[e]}
    gbt_matrix = sorted(e for e in gbt_rows if assignment[e] in matrices)
    if gbt_matrix:
        raise ValueError(f"table sections (section_matrix.txt) have no cross-section modes; elements {gbt_matrix} "
                         "ask for GBT")
    gbt_solid_shear = sorted(e for e in gbt_rows if assignment[e] in solid and shear_of[e])
    if gbt_solid_shear:
        raise ValueError(f"2D sections offer no extra shear-lag modes; elements {gbt_solid_shear} ask for [shear_modes]")
    not_gbt = sorted(e for e in gbt_rows if type_of and not has_modes(e))
    if not_gbt:
        raise ValueError(
            f"element_section.txt gives [modes] / [shear_modes] to elements whose type has no cross-section "
            f"modes ({', '.join(ElementFactory.types_with('cross_section_modes'))}): {not_gbt}"
        )
    for sid in {assignment[e] for e in gbt_rows} - set(solid):
        if analyses[sid].hoop != "restrained":
            a = analyses[sid]
            analyses[sid] = ThinWalledSectionAnalysis(a.section, a.laminates, hoop="restrained", n_sub=a.n_sub)
    max_modes = {sid: 0 for sid in (*analyses, *solid)}
    max_shear = {sid: 0 for sid in (*analyses, *solid)}
    gbt_sections = {assignment[e] for e in gbt_rows}
    for e in gbt_rows:
        sid = assignment[e]
        max_modes[sid] = max(max_modes[sid], modes_of[e])
        max_shear[sid] = max(max_shear[sid], shear_of[e])
    per_section = {
        sid: SectionAnalysis.from_thin_walled(
            a, sid, gbt_modes=max_modes[sid], gbt_shear_modes=max_shear[sid], gbt=sid in gbt_sections
        )
        for sid, a in analyses.items()
    }
    per_section.update({
        sid: SectionAnalysis.from_solid(a, sid, gbt_modes=max_modes[sid], gbt=sid in gbt_sections)
        for sid, a in solid.items()
    })
    per_section.update(matrices)
    per_element = tuple(per_section[assignment[e]] if e in assignment else None for e in ids)
    taper_path = os.path.join(job_dir, SECTION_TAPER_FILE)
    tapers = {}
    if os.path.exists(taper_path):
        from onedfem.pre_processing.parsing.section_taper_parser import (
            parse_section_taper_2d,
            parse_section_taper_file,
        )

        tapers = parse_section_taper_file(taper_path)
        rows_2d = parse_section_taper_2d(taper_path)
        unassigned = sorted(e for e in (*tapers, *rows_2d) if e not in assignment)
        if unassigned:
            raise ValueError(f"section_taper.txt names elements without a section in element_section.txt: {unassigned}")
        tabled = sorted(e for e in (*tapers, *rows_2d) if assignment[e] in matrices)
        if tabled:
            raise ValueError(f"section_taper.txt names elements with table sections (section_matrix.txt), whose "
                             f"tables already hold their stiffness: {tabled}")
        tapers.update(solid_tapers(rows_2d, tapers, solid))
        per_element = tuple(sa.with_taper(tapers[e]) if e in tapers else sa for e, sa in zip(ids, per_element))
    align = {}
    if gbt_rows and "connectivity" in element_dictionary and len({assignment[e] for e in gbt_rows}) > 1:
        per_element, align = align_gbt_sections(per_element, element_dictionary["connectivity"],
                                                [e in gbt_rows for e in ids])
    out = {
        **element_dictionary,
        SECTION_ANALYSIS_KEY: per_element,
        THIN_WALLED_KEY: tuple(analyses.get(assignment[e]) if e in assignment else None for e in ids),
    }
    if solid:
        out[SOLID_KEY] = tuple(solid.get(assignment[e]) if e in assignment else None for e in ids)
    if tapers:
        out[TAPER_KEY] = tuple(tapers.get(e) for e in ids)
    if gbt_rows:
        out["gbt_modes"] = np.array([modes_of.get(e, 0) for e in ids], dtype=int)
        out["gbt_shear_modes"] = np.array([shear_of.get(e, 0) for e in ids], dtype=int)
        out["extra_dofs"] = np.array(
            [
                per_section[assignment[e]].gbt_sv_modes + shear_of[e] + 2 * modes_of[e] if e in gbt_rows else 0
                for e in ids
            ],
            dtype=int,
        )

    wl_path = os.path.join(job_dir, WALL_LOAD_FILE)
    if os.path.exists(wl_path):
        from onedfem.pre_processing.parsing.wall_load_parser import parse_wall_load_file
        CrossSectionModes = section_import("thin_walled.cross_section_modes", "CrossSectionModes")

        loads = parse_wall_load_file(wl_path)
        not_thin = sorted(e for e in loads if e not in assignment or assignment[e] not in analyses)
        if not_thin:
            raise ValueError(f"wall_load.txt loads elements without a thin-walled section: {not_thin}")
        modes_cache = {}
        per_load = []
        for e in ids:
            if e not in loads:
                per_load.append(None)
                continue
            sid = assignment[e]
            if e in gbt_rows:
                if sid not in modes_cache:
                    modes_cache[sid] = CrossSectionModes(analyses[sid])
                f = modes_cache[sid].wall_load(loads[e], max_modes[sid], max_shear[sid])
                k = ids.index(e)
                if k in align:                       # the element's bases were turned (gbt_alignment)
                    mass_transform = section_import("gbt_alignment", "mass_transform")

                    f = mass_transform(per_element[k], align[k]["Rw"], align[k]["Rm"]).T @ f
                per_load.append(f)
            else:
                per_load.append(analyses[sid].wall_load_resultant(loads[e]))
        out["wall_load_per_element"] = tuple(per_load)
    return out

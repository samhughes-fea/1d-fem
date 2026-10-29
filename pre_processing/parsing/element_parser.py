# pre_processing\parsing\element_parser.py

from typing import Dict, List, Tuple
import numpy as np
import numpy.typing as npt

class ElementParser:
    """
    Parses an [Element] section and returns a dictionary whose structure is
    identical to the other parsers in the suite:

        {
            "element_dictionary": {
                "ids":               np.ndarray[int64],
                "connectivity":      np.ndarray[int64]  (N, 2),
                "types":             np.ndarray[str_],
                "integration_orders": {
                    "axial":        np.ndarray[int64],
                    "bending_y":    np.ndarray[int64],
                    "bending_z":    np.ndarray[int64],
                    "shear_y":      np.ndarray[int64],
                    "shear_z":      np.ndarray[int64],
                    "torsion":      np.ndarray[int64],
                    "load":         np.ndarray[int64],
                },
                "curvature":         np.ndarray[float64]  (legacy; always present, default 0),
                "warping":           np.ndarray[int8]     (0/1; always present; from optional [warping] column or 0),
                "twist":             np.ndarray[float64]  (section twist about the element axis, degrees; default 0),
            }
        }

    **Column layouts:** 11 standard columns, then any of the optional ``[curvature]``, ``[warping]``
    and ``[twist]`` columns in that order. ``[twist]`` turns the section about the element axis
    (``pre_processing/element_library/element_frame.py``).

    See ``docs/conventions/JOB_INPUT_BEAM_WARPING.md`` for beam warping columns, section tiers, and legacy type migration.
    Removed public type strings are listed in ``docs/conventions/DEPRECATED_ELEMENT_TYPES.md``.
    """

    def __init__(self, filepath: str, job_results_dir: str) -> None:
        self.filepath: str = filepath
        self.job_results_dir: str = job_results_dir
        self.output_filename: str = "elements_parsed.csv"

        self.expected_subheader: List[str] = [
            "[element_id]", "[node1]", "[node2]", "[element_type]",
            "[axial_order]", "[bending_y_order]", "[bending_z_order]",
            "[shear_y_order]", "[shear_z_order]", "[torsion_order]", "[load_order]",
        ]
        # Optional 12th column [curvature]: legacy scalar κ₀ per row; parsed into "curvature" array but
        # not used by element implementations — use straight beam types + precurvature.txt instead.
        self.optional_curvature_header: List[str] = self.expected_subheader + ["[curvature]"]
        # Optional [warping] only (12 columns): 0/1 per element — assemble Vlasov warping stiffness.
        self.optional_warping_header: List[str] = self.expected_subheader + ["[warping]"]
        # Optional [curvature] + [warping] (13 columns), fixed order.
        self.optional_curvature_warping_header: List[str] = (
            self.expected_subheader + ["[curvature]", "[warping]"]
        )

    # --------------------------------------------------------------------- #
    # Internal utilities
    # --------------------------------------------------------------------- #
    @staticmethod
    def _assert_exact_subheader(line: str, expected: List[str]) -> None:
        """
        Ensures the sub-header line matches *expected*, case-insensitively and
        ignoring extra whitespace.
        """
        tokens = [token.lower() for token in line.split()]
        if tokens != [hdr.lower() for hdr in expected]:
            raise ValueError(
                f"Sub-header must match (case-insensitive): {' '.join(expected)}"
            )

    @staticmethod
    def _parse_subheader(line: str) -> List[str]:
        """Return list of subheader tokens (lowercase). Used to detect optional columns."""
        return [t.lower() for t in line.split()]

    @staticmethod
    def _preprocess_lines(filepath: str) -> List[str]:
        """
        Reads a file and returns a list of stripped lines, skipping empty lines
        and those that start with '#' (comments).
        """
        with open(filepath, "r", encoding="utf-8") as fh:
            return [
                ln.strip()
                for ln in fh
                if ln.strip() and not ln.lstrip().startswith("#")
            ]

    OPTIONAL_COLUMNS = ("[curvature]", "[warping]", "[twist]")

    @staticmethod
    def _classify_subheader(subheader_line: str) -> Tuple[str, List[str]]:
        """
        Return ``(mode, optional_columns)`` for the data rows that follow.

        After the 11 standard columns come any of the optional ``[curvature]``,
        ``[warping]`` and ``[twist]`` columns, in that order. ``mode`` is ``base`` when there
        are none, else the optional column names joined by ``+``.
        """
        tokens = [t.lower() for t in subheader_line.split()]
        exp = [h.lower() for h in [
            "[element_id]", "[node1]", "[node2]", "[element_type]",
            "[axial_order]", "[bending_y_order]", "[bending_z_order]",
            "[shear_y_order]", "[shear_z_order]", "[torsion_order]", "[load_order]",
        ]]
        extra = tokens[len(exp):]
        order = list(ElementParser.OPTIONAL_COLUMNS)
        if tokens[:len(exp)] != exp or any(t not in order for t in extra) or \
                [order.index(t) for t in extra] != sorted(order.index(t) for t in extra) or len(set(extra)) != len(extra):
            raise ValueError(
                "Sub-header must be the 11 standard columns, optionally followed by any of "
                "[curvature], [warping] and [twist], in that order."
            )
        return ("base" if not extra else "+".join(extra)), extra

    @staticmethod
    def _parse_warping_token(tok: str) -> int:
        t = str(tok).strip().lower()
        if t in ("1", "true", "yes", "on"):
            return 1
        if t in ("0", "false", "no", "off"):
            return 0
        raise ValueError(f"Invalid [warping] value {tok!r}; use 0/1 or true/false.")

    # --------------------------------------------------------------------- #
    # Public API
    # --------------------------------------------------------------------- #
    def parse(self) -> Dict[str, Dict[str, npt.NDArray]]:
        element_ids:   List[int]        = []
        connectivity:  List[List[int]]  = []
        element_types: List[str]        = []
        integration_orders: Dict[str, List[int]] = {
            "axial":      [],
            "bending_y":  [],
            "bending_z":  [],
            "shear_y":    [],
            "shear_z":    [],
            "torsion":    [],
            "load":       [],
        }
        seen_ids: set[int] = set()

        # ------------------------------------------------------------------ #
        # Read & pre-clean the file
        # ------------------------------------------------------------------ #
        lines = self._preprocess_lines(self.filepath)

        # Locate [Element] (case-insensitive)
        try:
            start_idx = next(
                i for i, ln in enumerate(lines) if ln.lower() == "[element]"
            )
        except StopIteration:
            raise ValueError("Missing [Element] section header.")

        mode, extra = self._classify_subheader(lines[start_idx + 1])

        curvature_list: List[float] = []
        warping_list: List[int] = []
        twist_list: List[float] = []

        # ------------------------------------------------------------------ #
        # Parse each data line
        # ------------------------------------------------------------------ #
        expected_cols = 11 + len(extra)

        for ln in lines[start_idx + 2:]:
            parts = ln.split()
            if len(parts) != expected_cols:
                raise ValueError(f"Malformed element row: {ln!r}")

            try:
                eid               = int(parts[0])
                n1, n2            = map(int, parts[1:3])
                etype             = parts[3]
                orders            = list(map(int, parts[4:11]))
                values = dict(zip(extra, parts[11:]))
                curvature_list.append(float(values.get("[curvature]", 0.0)))
                warping_list.append(self._parse_warping_token(values.get("[warping]", "0")))
                twist_list.append(float(values.get("[twist]", 0.0)))
            except (ValueError, IndexError) as exc:
                raise TypeError(f"Bad data types in line {ln!r} → {exc}") from exc

            if eid in seen_ids:
                raise ValueError(f"Duplicate element_id: {eid}")
            seen_ids.add(eid)

            element_ids.append(eid)
            connectivity.append([n1, n2])
            element_types.append(etype)

            (
                integration_orders["axial"],
                integration_orders["bending_y"],
                integration_orders["bending_z"],
                integration_orders["shear_y"],
                integration_orders["shear_z"],
                integration_orders["torsion"],
                integration_orders["load"],
            ) = [
                lst + [val]  # append to the appropriate list
                for lst, val in zip(integration_orders.values(), orders)
            ]

        # ------------------------------------------------------------------ #
        # Convert everything to NumPy arrays (homogeneous outbound types)
        # ------------------------------------------------------------------ #
        ids_arr   = np.asarray(element_ids,   dtype=np.int64)
        conn_arr  = np.asarray(connectivity,  dtype=np.int64)
        types_arr = np.asarray(element_types, dtype="<U64")      # Long enough for Linear*/Nonlinear* type strings

        integ_np: Dict[str, npt.NDArray[np.int64]] = {
            k: np.asarray(v, dtype=np.int64) for k, v in integration_orders.items()
        }
        curvature_arr = np.asarray(curvature_list, dtype=np.float64)
        warping_arr = np.asarray(warping_list, dtype=np.int8)
        twist_arr = np.asarray(twist_list, dtype=np.float64)

        # ------------------------------------------------------------------ #
        # Return the uniform dictionary structure
        # ------------------------------------------------------------------ #
        return {
            "element_dictionary": {
                "ids":               ids_arr,
                "connectivity":      conn_arr,
                "types":             types_arr,
                "integration_orders": integ_np,
                "curvature":         curvature_arr,
                "warping":           warping_arr,
                "twist":             twist_arr,
            }
        }

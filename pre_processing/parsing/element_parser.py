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
            }
        }

    **Column layouts:** 11 standard columns; optional 12th ``[curvature]`` *or* ``[warping]`` (disambiguated by
    header); optional 13th when both ``[curvature]`` and ``[warping]`` are present (that order).

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

    @staticmethod
    def _classify_subheader(subheader_line: str) -> Tuple[str, List[str]]:
        """
        Return (mode, expected_tokens_lower) for the data rows that follow.

        Modes: ``base`` (11), ``curv`` (12 + curvature), ``warp`` (12 + warping),
        ``both`` (13 + curvature + warping).
        """
        tokens = [t.lower() for t in subheader_line.split()]
        exp = [h.lower() for h in [
            "[element_id]", "[node1]", "[node2]", "[element_type]",
            "[axial_order]", "[bending_y_order]", "[bending_z_order]",
            "[shear_y_order]", "[shear_z_order]", "[torsion_order]", "[load_order]",
        ]]
        if tokens == exp:
            return "base", tokens
        if tokens == exp + ["[curvature]"]:
            return "curv", tokens
        if tokens == exp + ["[warping]"]:
            return "warp", tokens
        if tokens == exp + ["[curvature]", "[warping]"]:
            return "both", tokens
        raise ValueError(
            "Sub-header must be 11 columns, or add optional [curvature], "
            "or [warping], or both [curvature] then [warping] (13 columns)."
        )

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

        mode, _ = self._classify_subheader(lines[start_idx + 1])
        if mode == "base":
            self._assert_exact_subheader(lines[start_idx + 1], self.expected_subheader)
        elif mode == "curv":
            self._assert_exact_subheader(lines[start_idx + 1], self.optional_curvature_header)
        elif mode == "warp":
            self._assert_exact_subheader(lines[start_idx + 1], self.optional_warping_header)
        else:
            self._assert_exact_subheader(lines[start_idx + 1], self.optional_curvature_warping_header)

        curvature_list: List[float] = []
        warping_list: List[int] = []

        # ------------------------------------------------------------------ #
        # Parse each data line
        # ------------------------------------------------------------------ #
        if mode == "base":
            expected_cols = 11
        elif mode in ("curv", "warp"):
            expected_cols = 12
        else:
            expected_cols = 13

        for ln in lines[start_idx + 2:]:
            parts = ln.split()
            if len(parts) != expected_cols:
                raise ValueError(f"Malformed element row: {ln!r}")

            try:
                eid               = int(parts[0])
                n1, n2            = map(int, parts[1:3])
                etype             = parts[3]
                orders            = list(map(int, parts[4:11]))
                curv = 0.0
                warp = 0
                if mode == "curv":
                    curvature_list.append(float(parts[11]))
                    warping_list.append(0)
                elif mode == "warp":
                    curvature_list.append(0.0)
                    warping_list.append(self._parse_warping_token(parts[11]))
                elif mode == "both":
                    curvature_list.append(float(parts[11]))
                    warping_list.append(self._parse_warping_token(parts[12]))
                else:
                    curvature_list.append(0.0)
                    warping_list.append(0)
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
            }
        }

# pre_processing\parsing\element_parser.py

from typing import Dict, List
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
            }
        }
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
        # Optional 12th column for curved beams (Phase 2a)
        self.optional_curvature_header: List[str] = self.expected_subheader + ["[curvature]"]

    # --------------------------------------------------------------------- #
    # Internal utilities
    # --------------------------------------------------------------------- #
    @staticmethod
    def _assert_exact_subheader(line: str, expected: List[str]) -> None:
        """
        Ensures the sub-header line matches `expected`, case-insensitively and
        ignoring extra whitespace.
        """
        tokens = [token.lower() for token in line.split()]
        if tokens != [hdr.lower() for hdr in expected]:
            raise ValueError(
                f"Sub-header must match (case-insensitive): {' '.join(expected)}"
            )

    @staticmethod
    def _parse_subheader(line: str) -> List[str]:
        """Return list of subheader tokens (lowercase). Used to detect optional [curvature] column."""
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

        # Detect sub-header: 11 columns (standard) or 12 with [curvature] (Phase 2a)
        subheader_tokens = self._parse_subheader(lines[start_idx + 1])
        if len(subheader_tokens) == len(self.optional_curvature_header) and subheader_tokens[-1] == "[curvature]":
            self._assert_exact_subheader(lines[start_idx + 1], self.optional_curvature_header)
            use_curvature = True
        else:
            self._assert_exact_subheader(lines[start_idx + 1], self.expected_subheader)
            use_curvature = False

        curvature_list: List[float] = []

        # ------------------------------------------------------------------ #
        # Parse each data line
        # ------------------------------------------------------------------ #
        expected_cols = len(self.optional_curvature_header) if use_curvature else len(self.expected_subheader)
        for ln in lines[start_idx + 2:]:
            parts = ln.split()
            if len(parts) != expected_cols:
                raise ValueError(f"Malformed element row: {ln!r}")

            try:
                eid               = int(parts[0])
                n1, n2            = map(int, parts[1:3])
                etype             = parts[3]
                orders            = list(map(int, parts[4:11]))
                if use_curvature:
                    curvature_list.append(float(parts[11]))
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

        if not use_curvature:
            curvature_list = [0.0] * len(element_ids)

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
            }
        }
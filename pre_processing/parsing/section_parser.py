# pre_processing\parsing\section_parser.py

import os
from typing import Dict, List
import numpy as np
import numpy.typing as npt
import pandas as pd


class SectionParser:
    """
    Parses a [Section] block and returns

        {
            "section_dictionary": {
                "element_id": np.ndarray[int64],
                "A":          np.ndarray[float64],
                "I_x":        np.ndarray[float64],
                "I_y":        np.ndarray[float64],
                "I_z":        np.ndarray[float64],
                "J_t":        np.ndarray[float64],
                "kappa":      np.ndarray[float64],  # optional; shear correction factor (Timoshenko)
                "alpha":      np.ndarray[float64],  # optional; higher-order shear coeff (Levinson)
            }
        }

    Sub-header may be 6 columns (element_id, A, I_x, I_y, I_z, J_t) or 8 columns
    with [kappa] and [alpha] for general-section preprocessing (geometry-derived).
    """

    def __init__(self, filepath: str, job_results_dir: str) -> None:
        self.filepath: str = filepath
        self.job_results_dir: str = job_results_dir
        self.output_filename: str = "section_properties_parsed.csv"
        self.expected_subheader_6: List[str] = [
            "[element_id]", "[A]", "[I_x]", "[I_y]", "[I_z]", "[J_t]"
        ]
        self.expected_subheader_8: List[str] = [
            "[element_id]", "[A]", "[I_x]", "[I_y]", "[I_z]", "[J_t]",
            "[kappa]", "[alpha]"
        ]

    # ------------------------------------------------------------------ #
    @staticmethod
    def _assert_exact_subheader(line: str, expected: List[str]) -> None:
        tokens = [tok.lower() for tok in line.split()]
        if tokens != [hdr.lower() for hdr in expected]:
            raise ValueError(
                f"Sub-header must match (case-insensitive): {' '.join(expected)}"
            )

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

    # ------------------------------------------------------------------ #
    def parse(self) -> Dict[str, Dict[str, npt.NDArray]]:
        lists: Dict[str, List[float]] = {
            "element_id": [],
            "A":  [],
            "I_x": [],
            "I_y": [],
            "I_z": [],
            "J_t": [],
        }
        has_kappa_alpha = False
        seen_ids: set[int] = set()

        # ---- Read & clean ---------------------------------------------- #
        lines = self._preprocess_lines(self.filepath)

        # Locate [Section]
        try:
            start_idx = next(i for i, ln in enumerate(lines) if ln.lower() == "[section]")
        except StopIteration:
            raise ValueError("Missing [Section] section header.")

        subheader_line = lines[start_idx + 1]
        subheader_tokens = [t.lower() for t in subheader_line.split()]
        if subheader_tokens == [h.lower() for h in self.expected_subheader_8]:
            has_kappa_alpha = True
            lists["kappa"] = []
            lists["alpha"] = []
            self._assert_exact_subheader(subheader_line, self.expected_subheader_8)
        elif subheader_tokens == [h.lower() for h in self.expected_subheader_6]:
            self._assert_exact_subheader(subheader_line, self.expected_subheader_6)
        else:
            raise ValueError(
                f"Sub-header must be 6 columns {self.expected_subheader_6} or "
                f"8 columns {self.expected_subheader_8} (case-insensitive)."
            )

        n_cols = 8 if has_kappa_alpha else 6

        # ---- Parse rows ------------------------------------------------- #
        for ln in lines[start_idx + 2:]:
            parts = ln.split()
            if len(parts) != n_cols:
                raise ValueError(f"Malformed section row (expected {n_cols} columns): {ln!r}")

            try:
                eid = int(parts[0])
                A, Ix, Iy, Iz, Jt = map(float, parts[1:6])
            except ValueError as exc:
                raise TypeError(f"Bad data types in line {ln!r} → {exc}") from exc

            if eid in seen_ids:
                raise ValueError(f"Duplicate element_id: {eid}")
            seen_ids.add(eid)

            lists["element_id"].append(eid)
            lists["A"].append(A)
            lists["I_x"].append(Ix)
            lists["I_y"].append(Iy)
            lists["I_z"].append(Iz)
            lists["J_t"].append(Jt)

            if has_kappa_alpha:
                try:
                    kappa, alpha = float(parts[6]), float(parts[7])
                except ValueError as exc:
                    raise TypeError(f"Bad kappa/alpha in line {ln!r} → {exc}") from exc
                lists["kappa"].append(kappa)
                lists["alpha"].append(alpha)

        # ---- Convert to NumPy arrays ----------------------------------- #
        parsed: Dict[str, npt.NDArray] = {
            k: np.asarray(v, dtype=(np.int64 if k == "element_id" else np.float64))
            for k, v in lists.items()
        }

        # ---- Uniform return structure ---------------------------------- #
        return {"section_dictionary": parsed}
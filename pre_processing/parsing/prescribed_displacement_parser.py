# pre_processing\parsing\prescribed_displacement_parser.py

import numpy as np
import logging
import re
import os

# Set up detailed logging
logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")

# DOF name to index mapping
DOF_MAP = {
    "UX": 0, "UY": 1, "UZ": 2,
    "RX": 3, "RY": 4, "RZ": 5
}

def parse_prescribed_displacement(file_path):
    """
    Parses prescribed displacement conditions from a structured text file.
    
    File format:
    [Prescribed Displacement]
    [id] [node_id] [dof] [value] [type] [comment]
    
    Parameters
    ----------
    file_path : str
        Path to the prescribed_displacement.txt file
        
    Returns
    -------
    dict
        Dictionary with keys:
        - 'id': np.ndarray of condition IDs
        - 'node_id': np.ndarray of node IDs
        - 'dof': np.ndarray of DOF names (UX, UY, UZ, RX, RY, RZ)
        - 'dof_index': np.ndarray of DOF indices (0-5)
        - 'value': np.ndarray of prescribed displacement values
        - 'type': np.ndarray of condition types
        - 'global_dof': np.ndarray of global DOF indices (node_id * 6 + dof_index)
    """
    # Check if the file exists
    if not os.path.exists(file_path):
        logging.warning(f"[Prescribed Displacement] File not found: {file_path}. Returning empty arrays.")
        return {
            'id': np.array([], dtype=int),
            'node_id': np.array([], dtype=int),
            'dof': np.array([], dtype=str),
            'dof_index': np.array([], dtype=int),
            'value': np.array([], dtype=float),
            'type': np.array([], dtype=str),
            'global_dof': np.array([], dtype=int)
        }

    logging.info(f"[Prescribed Displacement] Reading file: {file_path}")

    ids = []
    node_ids = []
    dofs = []
    dof_indices = []
    values = []
    types = []
    
    header_pattern = re.compile(r"^\[Prescribed Displacement\]$", re.IGNORECASE)
    current_section = False
    first_data_line_detected = False

    with open(file_path, 'r') as f:
        for line_number, raw_line in enumerate(f, 1):
            line = raw_line.split("#")[0].strip()  # Remove inline comments

            if not line:
                continue  # Skip empty lines

            # Detect the [Prescribed Displacement] section
            if header_pattern.match(line):
                logging.info(f"[Prescribed Displacement] Found section at line {line_number}")
                current_section = True
                continue

            # Skip any data before section
            if not current_section:
                continue

            # Skip header line (contains non-numeric values like [id], [node_id], etc.)
            if not first_data_line_detected:
                # Check if line contains column headers (non-numeric except for potential numeric IDs)
                if any(re.search(r"[^\d\.\-+eE\s]", p) for p in line.split() if p.startswith('[')):
                    logging.debug(f"[Prescribed Displacement] Skipping header line {line_number}")
                    continue
                first_data_line_detected = True

            # Parse data line
            parts = line.split()
            if len(parts) < 5:
                logging.warning(f"[Prescribed Displacement] Line {line_number}: Expected at least 5 values, found {len(parts)}. Skipping.")
                continue

            try:
                condition_id = int(parts[0])
                node_id = int(parts[1])
                dof_name = parts[2].upper()
                value = float(parts[3])
                condition_type = parts[4].lower()

                # Validate DOF name
                if dof_name not in DOF_MAP:
                    logging.warning(f"[Prescribed Displacement] Line {line_number}: Invalid DOF '{parts[2]}'. Expected one of {list(DOF_MAP.keys())}. Skipping.")
                    continue

                dof_index = DOF_MAP[dof_name]
                global_dof = node_id * 6 + dof_index

                ids.append(condition_id)
                node_ids.append(node_id)
                dofs.append(dof_name)
                dof_indices.append(dof_index)
                values.append(value)
                types.append(condition_type)

                logging.debug(f"[Prescribed Displacement] Line {line_number}: Node {node_id}, DOF {dof_name} (index {dof_index}), value {value}, global DOF {global_dof}")

            except (ValueError, IndexError) as e:
                logging.warning(f"[Prescribed Displacement] Line {line_number}: Parsing error: {e}. Skipping.")
                continue

    # Handle case where no valid entries were found
    if not ids:
        logging.warning(f"[Prescribed Displacement] No valid entries found in '{file_path}'. Returning empty arrays.")
        return {
            'id': np.array([], dtype=int),
            'node_id': np.array([], dtype=int),
            'dof': np.array([], dtype=str),
            'dof_index': np.array([], dtype=int),
            'value': np.array([], dtype=float),
            'type': np.array([], dtype=str),
            'global_dof': np.array([], dtype=int)
        }

    # Convert to NumPy arrays
    result = {
        'id': np.array(ids, dtype=int),
        'node_id': np.array(node_ids, dtype=int),
        'dof': np.array(dofs, dtype=str),
        'dof_index': np.array(dof_indices, dtype=int),
        'value': np.array(values, dtype=float),
        'type': np.array(types, dtype=str),
        'global_dof': np.array([node_id * 6 + dof_idx for node_id, dof_idx in zip(node_ids, dof_indices)], dtype=int)
    }

    logging.info(f"[Prescribed Displacement] Successfully parsed {len(ids)} prescribed displacement conditions from '{file_path}'.")
    logging.debug(f"[Prescribed Displacement] Summary: {len(np.unique(result['node_id']))} unique nodes, {len(np.unique(result['dof']))} unique DOFs")

    return result

# Standalone execution for testing
if __name__ == "__main__":
    test_file = r"sandbox\job_new\prescribed_displacement.txt"
    if not os.path.exists(test_file):
        logging.error(f"Test file '{test_file}' not found. Make sure it exists before running.")
    else:
        try:
            output = parse_prescribed_displacement(test_file)
            print("\n-------------Parsed [Prescribed Displacement] Data-------------\n")
            for key, value in output.items():
                print(f"{key}: {value}")
        except Exception as e:
            logging.error(f"Error parsing prescribed displacement file: {e}")


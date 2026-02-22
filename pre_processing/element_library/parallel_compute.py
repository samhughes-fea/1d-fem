# pre_processing/element_library/parallel_compute.py

"""
Parallel computation utilities for element operations.

Provides parallel computation of element stiffness matrices and force vectors
with automatic threshold detection and graceful fallback to sequential processing.
"""

import logging
import multiprocessing
import os
import pickle
from typing import List, Any, Callable, Optional
import numpy as np

# Configure module-level logger
logger = logging.getLogger(__name__)

# Default threshold for parallel processing (minimum number of elements)
DEFAULT_PARALLEL_THRESHOLD = 50


def _worker_stiffness(args):
    """Worker function for parallel stiffness computation."""
    element, index = args
    try:
        if element is None:
            return index, None
        result = element.element_stiffness_matrix()
        return index, result
    except Exception as e:
        logger.error(f"Error computing stiffness for element at index {index}: {e}")
        return index, None


def _worker_force(args):
    """Worker function for parallel force computation."""
    element, index = args
    try:
        if element is None:
            return index, None
        result = element.element_force_vector()
        return index, result
    except Exception as e:
        logger.error(f"Error computing force for element at index {index}: {e}")
        return index, None


def _sequential_stiffness_compute(elements: List[Any]) -> np.ndarray:
    """Sequential fallback for stiffness computation."""
    logger.debug("Using sequential stiffness computation")
    results = []
    for i, elem in enumerate(elements):
        try:
            if elem is None:
                results.append(None)
            else:
                results.append(elem.element_stiffness_matrix())
        except Exception as e:
            logger.error(f"Error computing stiffness for element {i}: {e}")
            results.append(None)
    return np.array(results, dtype=object)


def _sequential_force_compute(elements: List[Any]) -> np.ndarray:
    """Sequential fallback for force computation."""
    logger.debug("Using sequential force computation")
    results = []
    for i, elem in enumerate(elements):
        try:
            if elem is None:
                results.append(None)
            else:
                results.append(elem.element_force_vector())
        except Exception as e:
            logger.error(f"Error computing force for element {i}: {e}")
            results.append(None)
    return np.array(results, dtype=object)


def compute_element_stiffness_parallel(
    elements: List[Any],
    num_processes: Optional[int] = None,
    threshold: int = DEFAULT_PARALLEL_THRESHOLD
) -> np.ndarray:
    """
    Compute element stiffness matrices in parallel.

    Parameters
    ----------
    elements : List[Any]
        List of element objects with element_stiffness_matrix() method.
    num_processes : int, optional
        Number of processes to use. If None or "auto", uses os.cpu_count().
        If 1 or elements < threshold, uses sequential computation.
    threshold : int, optional
        Minimum number of elements to use parallel processing (default 50).

    Returns
    -------
    np.ndarray
        Array of element stiffness matrices (dtype=object), maintaining order.
    """
    num_elements = len(elements)
    
    # Determine if we should use parallel processing
    if num_elements < threshold:
        logger.debug(f"Element count ({num_elements}) < threshold ({threshold}), using sequential")
        return _sequential_stiffness_compute(elements)
    
    # Determine number of processes
    if num_processes is None or num_processes == "auto":
        num_processes = os.cpu_count() or 1
    
    if num_processes == 1:
        logger.debug("num_processes=1, using sequential")
        return _sequential_stiffness_compute(elements)
    
    # Try parallel computation with error handling
    try:
        logger.info(f"Computing {num_elements} element stiffness matrices in parallel ({num_processes} processes)")
        
        # Prepare arguments with indices to maintain order
        args = [(elem, i) for i, elem in enumerate(elements)]
        
        with multiprocessing.Pool(processes=num_processes) as pool:
            results = pool.map(_worker_stiffness, args)
        
        # Sort by index and extract results
        results.sort(key=lambda x: x[0])
        result_array = np.array([r[1] for r in results], dtype=object)
        
        logger.info("Parallel stiffness computation completed successfully")
        return result_array
        
    except (pickle.PicklingError, AttributeError, TypeError) as e:
        logger.warning(f"Parallel computation failed ({type(e).__name__}: {e}), falling back to sequential")
        return _sequential_stiffness_compute(elements)
    except Exception as e:
        logger.error(f"Unexpected error in parallel computation: {e}, falling back to sequential")
        return _sequential_stiffness_compute(elements)


def compute_element_force_parallel(
    elements: List[Any],
    num_processes: Optional[int] = None,
    threshold: int = DEFAULT_PARALLEL_THRESHOLD
) -> np.ndarray:
    """
    Compute element force vectors in parallel.

    Parameters
    ----------
    elements : List[Any]
        List of element objects with element_force_vector() method.
    num_processes : int, optional
        Number of processes to use. If None or "auto", uses os.cpu_count().
        If 1 or elements < threshold, uses sequential computation.
    threshold : int, optional
        Minimum number of elements to use parallel processing (default 50).

    Returns
    -------
    np.ndarray
        Array of element force vectors (dtype=object), maintaining order.
    """
    num_elements = len(elements)
    
    # Determine if we should use parallel processing
    if num_elements < threshold:
        logger.debug(f"Element count ({num_elements}) < threshold ({threshold}), using sequential")
        return _sequential_force_compute(elements)
    
    # Determine number of processes
    if num_processes is None or num_processes == "auto":
        num_processes = os.cpu_count() or 1
    
    if num_processes == 1:
        logger.debug("num_processes=1, using sequential")
        return _sequential_force_compute(elements)
    
    # Try parallel computation with error handling
    try:
        logger.info(f"Computing {num_elements} element force vectors in parallel ({num_processes} processes)")
        
        # Prepare arguments with indices to maintain order
        args = [(elem, i) for i, elem in enumerate(elements)]
        
        with multiprocessing.Pool(processes=num_processes) as pool:
            results = pool.map(_worker_force, args)
        
        # Sort by index and extract results
        results.sort(key=lambda x: x[0])
        result_array = np.array([r[1] for r in results], dtype=object)
        
        logger.info("Parallel force computation completed successfully")
        return result_array
        
    except (pickle.PicklingError, AttributeError, TypeError) as e:
        logger.warning(f"Parallel computation failed ({type(e).__name__}: {e}), falling back to sequential")
        return _sequential_force_compute(elements)
    except Exception as e:
        logger.error(f"Unexpected error in parallel computation: {e}, falling back to sequential")
        return _sequential_force_compute(elements)


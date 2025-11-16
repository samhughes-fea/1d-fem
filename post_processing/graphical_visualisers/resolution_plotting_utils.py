# post_processing/graphical_visualisers/resolution_plotting_utils.py

"""
Resolution-based plotting utilities for FEM visualization.

This module provides functions to plot quantities at different resolution levels:
- Gaussian: Discrete markers at Gauss point locations
- Nodal: Markers at node locations
- Interpolated: Continuous fields using shape function interpolation

Coordinate Systems:
- Natural coordinates: ξ ∈ [-1, 1] (Gauss points, node natural positions)
- Physical coordinates: x ∈ [0, L] (actual spatial positions)
- Transformation: x = x_node1 + (ξ + 1) * L / 2
"""

from typing import List, Optional, Tuple, Dict, Any, Union
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.axes
import importlib


def get_element_node_coords(
    element_id: int,
    element_dictionary: Dict[str, Any],
    grid_dictionary: Dict[str, Any]
) -> np.ndarray:
    """
    Get physical coordinates of element nodes.
    
    Parameters
    ----------
    element_id : int
        Element identifier
    element_dictionary : Dict[str, Any]
        Element connectivity dictionary with "ids" and "connectivity" arrays
    grid_dictionary : Dict[str, Any]
        Grid/node dictionary with coordinates
    
    Returns
    -------
    np.ndarray
        Node coordinates array of shape (n_nodes, 3) with columns [x, y, z]
    """
    # Find element index
    try:
        idx = int(np.where(element_dictionary["ids"] == element_id)[0][0])
    except (IndexError, KeyError) as exc:
        raise ValueError(f"Element {element_id} not found in element_dictionary") from exc
    
    # Get node IDs for this element
    node_ids = element_dictionary["connectivity"][idx]  # Shape: (n_nodes,)
    
    # Get coordinates from grid dictionary
    if "coordinates" in grid_dictionary:
        coords = grid_dictionary["coordinates"]
    elif "grid_dictionary" in grid_dictionary and "coordinates" in grid_dictionary["grid_dictionary"]:
        coords = grid_dictionary["grid_dictionary"]["coordinates"]
    else:
        raise KeyError("Cannot find coordinates in grid_dictionary")
    
    # Extract coordinates for element nodes
    node_coords = coords[node_ids].astype(np.float64, copy=False)
    
    return node_coords


def get_element_length(
    element_id: int,
    element_dictionary: Dict[str, Any],
    grid_dictionary: Dict[str, Any]
) -> float:
    """
    Compute element length from node coordinates.
    
    Parameters
    ----------
    element_id : int
        Element identifier
    element_dictionary : Dict[str, Any]
        Element connectivity dictionary
    grid_dictionary : Dict[str, Any]
        Grid/node dictionary with coordinates
    
    Returns
    -------
    float
        Element length (Euclidean distance between nodes)
    """
    node_coords = get_element_node_coords(element_id, element_dictionary, grid_dictionary)
    
    if len(node_coords) < 2:
        raise ValueError(f"Element {element_id} must have at least 2 nodes")
    
    # Compute length as distance between first and last node
    # For 1D beams, typically use x-coordinate difference or full 3D distance
    # Using 3D distance for generality
    length = np.linalg.norm(node_coords[-1] - node_coords[0])
    
    return float(length)


def get_element_type(
    element_id: int,
    element_dictionary: Dict[str, Any]
) -> str:
    """
    Get element type string for a given element.
    
    Parameters
    ----------
    element_id : int
        Element identifier
    element_dictionary : Dict[str, Any]
        Element dictionary with "ids" and "types" arrays
    
    Returns
    -------
    str
        Element type string (e.g., "EulerBernoulliBeamElement3D", "TimoshenkoBeamElement3D")
    """
    try:
        idx = int(np.where(element_dictionary["ids"] == element_id)[0][0])
    except (IndexError, KeyError) as exc:
        raise ValueError(f"Element {element_id} not found in element_dictionary") from exc
    
    if "types" not in element_dictionary:
        raise KeyError("element_dictionary must contain 'types' array")
    
    element_type = element_dictionary["types"][idx]
    return str(element_type)


def create_shape_function_operator(
    element_type: str,
    element_length: float
) -> Any:
    """
    Create the appropriate ShapeFunctionOperator for an element type.
    
    Maps element type strings to their corresponding ShapeFunctionOperator classes
    and instantiates them with the element length.
    
    Parameters
    ----------
    element_type : str
        Element type string (e.g., "EulerBernoulliBeamElement3D", "TimoshenkoBeamElement3D")
    element_length : float
        Physical length of the element
    
    Returns
    -------
    ShapeFunctionOperator
        Appropriate shape function operator instance
    
    Raises
    ------
    ValueError
        If element type is not recognized or cannot be imported
    """
    # Mapping from element type strings to their shape function module paths
    # This mirrors the structure in element_factory.py
    SHAPE_FUNCTION_MODULE_MAP = {
        "EulerBernoulliBeamElement3D": "pre_processing.element_library.euler_bernoulli.utilities.shape_functions",
        "TimoshenkoBeamElement3D": "pre_processing.element_library.timoshenko.utilities.shape_functions",
        "LevinsonBeamElement3D": "pre_processing.element_library.levinson.utilities.shape_functions",
        "EulerBernoulliBeamElement6DOF": "pre_processing.element_library.co_tide_beam_ML.utilities.shape_functions",
        # Add more mappings as needed
    }
    
    # Try to find the module path
    if element_type not in SHAPE_FUNCTION_MODULE_MAP:
        # Try to infer from element type string
        element_type_lower = element_type.lower()
        if "euler" in element_type_lower or "bernoulli" in element_type_lower:
            module_path = "pre_processing.element_library.euler_bernoulli.utilities.shape_functions"
        elif "timoshenko" in element_type_lower:
            module_path = "pre_processing.element_library.timoshenko.utilities.shape_functions"
        elif "levinson" in element_type_lower:
            module_path = "pre_processing.element_library.levinson.utilities.shape_functions"
        else:
            raise ValueError(
                f"Unknown element type: {element_type}. "
                f"Supported types: {list(SHAPE_FUNCTION_MODULE_MAP.keys())}"
            )
    else:
        module_path = SHAPE_FUNCTION_MODULE_MAP[element_type]
    
    # Import the module and get ShapeFunctionOperator class
    try:
        module = importlib.import_module(module_path)
        if not hasattr(module, "ShapeFunctionOperator"):
            raise AttributeError(f"Module {module_path} does not contain ShapeFunctionOperator class")
        ShapeFunctionOperator = module.ShapeFunctionOperator
    except ImportError as exc:
        raise ValueError(f"Cannot import shape function module {module_path}: {exc}") from exc
    
    # Create and return the operator
    return ShapeFunctionOperator(element_length=element_length)


def get_shape_function_operator_for_element(
    element_id: int,
    element_dictionary: Dict[str, Any],
    grid_dictionary: Dict[str, Any],
    element_object: Optional[Any] = None
) -> Any:
    """
    Get the correct ShapeFunctionOperator for an element.
    
    This function automatically determines the element type and creates
    the appropriate shape function operator. It can use either:
    1. element_object.element_type (if provided)
    2. element_dictionary["types"] (if element_object not provided)
    
    Parameters
    ----------
    element_id : int
        Element identifier
    element_dictionary : Dict[str, Any]
        Element dictionary with connectivity and types
    grid_dictionary : Dict[str, Any]
        Grid/node dictionary for computing element length
    element_object : Optional[ElementObject]
        Optional element object with element_type attribute.
        If provided, uses element_object.element_type instead of looking up in dictionary.
    
    Returns
    -------
    ShapeFunctionOperator
        Appropriate shape function operator for this element
    """
    # Get element type
    if element_object is not None and hasattr(element_object, 'element_type'):
        element_type = element_object.element_type
    else:
        element_type = get_element_type(element_id, element_dictionary)
    
    # Get element length
    element_length = get_element_length(element_id, element_dictionary, grid_dictionary)
    
    # Create and return the operator
    return create_shape_function_operator(element_type, element_length)


def natural_to_physical_coords(
    xi: Union[float, np.ndarray],
    node_coords: np.ndarray
) -> Union[float, np.ndarray]:
    """
    Transform natural coordinates to physical coordinates.
    
    For 1D beam elements: x = x_node1 + (ξ + 1) * L / 2
    where L is the element length.
    
    Parameters
    ----------
    xi : float or np.ndarray
        Natural coordinate(s) ∈ [-1, 1]
    node_coords : np.ndarray
        Node coordinates array of shape (n_nodes, 3)
        First row is node 1, last row is node 2
    
    Returns
    -------
    float or np.ndarray
        Physical coordinate(s) along element axis
        For 1D: returns x-coordinate
        Shape matches input xi
    """
    xi = np.asarray(xi)
    is_scalar = xi.ndim == 0
    xi = np.atleast_1d(xi)
    
    if len(node_coords) < 2:
        raise ValueError("Need at least 2 nodes for coordinate transformation")
    
    # Get first and last node coordinates
    x_node1 = node_coords[0, 0]  # x-coordinate of first node
    x_node2 = node_coords[-1, 0]  # x-coordinate of last node
    L = x_node2 - x_node1
    
    # Transform: x = x_node1 + (ξ + 1) * L / 2
    x_physical = x_node1 + (xi + 1.0) * L / 2.0
    
    if is_scalar:
        return float(x_physical[0])
    return x_physical


def get_gauss_point_physical_coords(
    gauss_data: List[Any],
    node_coords: np.ndarray
) -> np.ndarray:
    """
    Convert Gauss point natural coordinates to physical coordinates.
    
    Parameters
    ----------
    gauss_data : List[GaussPointData]
        List of Gauss point data objects with 'xi' attribute
    node_coords : np.ndarray
        Element node coordinates of shape (n_nodes, 3)
    
    Returns
    -------
    np.ndarray
        Physical coordinates of Gauss points, shape (n_gauss,)
    """
    # Extract natural coordinates
    xi_gauss = [gp.xi for gp in gauss_data]
    xi_array = np.array(xi_gauss)
    
    # Transform to physical coordinates
    x_physical = natural_to_physical_coords(xi_array, node_coords)
    
    return x_physical


def interpolate_field_nodal_to_continuous(
    nodal_values: np.ndarray,
    element_id: int,
    element_dictionary: Dict[str, Any],
    grid_dictionary: Dict[str, Any],
    shape_function_operator: Optional[Any] = None,
    element_object: Optional[Any] = None,
    n_points: int = 100
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Interpolate nodal field values to continuous field using shape functions.
    
    Parameters
    ----------
    nodal_values : np.ndarray
        Field values at nodes, shape (n_nodes, n_components) or (n_nodes,)
    element_id : int
        Element identifier
    element_dictionary : Dict[str, Any]
        Element connectivity dictionary
    grid_dictionary : Dict[str, Any]
        Grid/node dictionary
    shape_function_operator : Optional[ShapeFunctionOperator]
        Shape function operator for the element type.
        If None, will be automatically retrieved based on element type.
    element_object : Optional[ElementObject]
        Optional element object. If provided and shape_function_operator is None,
        will be used to determine element type.
    n_points : int, optional
        Number of interpolation points (default: 100)
    
    Returns
    -------
    Tuple[np.ndarray, np.ndarray]
        (x_coords, interpolated_values)
        - x_coords: Physical coordinates, shape (n_points,)
        - interpolated_values: Interpolated field, shape (n_points, n_components) or (n_points,)
    """
    # Get shape function operator if not provided
    if shape_function_operator is None:
        shape_function_operator = get_shape_function_operator_for_element(
            element_id, element_dictionary, grid_dictionary, element_object
        )
    # Get node coordinates
    node_coords = get_element_node_coords(element_id, element_dictionary, grid_dictionary)
    n_nodes = len(node_coords)
    
    # Get natural coordinates of nodes
    if n_nodes == 2:
        xi_nodes = np.array([-1.0, 1.0])
    elif n_nodes == 3:
        xi_nodes = np.array([-1.0, 0.0, 1.0])
    else:
        # General case: equally spaced in [-1, 1]
        xi_nodes = np.linspace(-1.0, 1.0, n_nodes)
    
    # Create fine grid of natural coordinates for interpolation
    xi_interp = np.linspace(-1.0, 1.0, n_points)
    
    # Evaluate shape functions at interpolation points
    N, _, _ = shape_function_operator.natural_coordinate_form(xi_interp)
    
    # Determine which DOF to use for interpolation
    # Shape function matrix N has shape (n_points, 12, 6):
    # - 12 DOFs: 6 per node (u_x, u_y, u_z, θ_x, θ_y, θ_z) for 2 nodes
    # - 6 components: physical components
    # For node i, DOF j: N[:, i*6 + j, j] gives shape function for that DOF/component
    
    nodal_values = np.asarray(nodal_values)
    if nodal_values.ndim == 1:
        # Scalar field: shape (n_nodes,)
        # Assume scalar field corresponds to first component (axial displacement, u_x)
        # Use shape functions for axial DOF: N[:, node*6 + 0, 0]
        N_interp = np.zeros((n_points, n_nodes))
        for i in range(n_nodes):
            dof_idx = i * 6  # Axial DOF (u_x) for node i
            component_idx = 0  # Axial component
            N_interp[:, i] = N[:, dof_idx, component_idx]
        
        # Interpolate: u(xi) = sum(N_i(xi) * u_i)
        interpolated = np.dot(N_interp, nodal_values)
    else:
        # Vector field: shape (n_nodes, n_components)
        n_components = nodal_values.shape[1]
        interpolated = np.zeros((n_points, n_components))
        
        # Map component index to DOF index within a node
        # Component 0: u_x (DOF 0)
        # Component 1: u_y (DOF 1)
        # Component 2: u_z (DOF 2)
        # Component 3: θ_x (DOF 3)
        # Component 4: θ_y (DOF 4)
        # Component 5: θ_z (DOF 5)
        
        for comp in range(n_components):
            N_interp = np.zeros((n_points, n_nodes))
            for i in range(n_nodes):
                dof_idx = i * 6 + comp  # DOF index for node i, component comp
                N_interp[:, i] = N[:, dof_idx, comp]
            
            interpolated[:, comp] = np.dot(N_interp, nodal_values[:, comp])
    
    # Transform natural coordinates to physical
    x_coords = natural_to_physical_coords(xi_interp, node_coords)
    
    return x_coords, interpolated


def plot_gauss_points(
    ax: matplotlib.axes.Axes,
    gauss_coords: np.ndarray,
    values: np.ndarray,
    marker: str = 'o',
    color: str = 'red',
    size: float = 50.0,
    alpha: float = 0.8,
    label: Optional[str] = None,
    **kwargs
) -> None:
    """
    Plot discrete Gauss point markers.
    
    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Matplotlib axes to plot on
    gauss_coords : np.ndarray
        Physical coordinates of Gauss points, shape (n_gauss,)
    values : np.ndarray
        Field values at Gauss points, shape (n_gauss,) or (n_gauss, n_components)
    marker : str, optional
        Marker style (default: 'o')
    color : str, optional
        Marker color (default: 'red')
    size : float, optional
        Marker size (default: 50.0)
    alpha : float, optional
        Marker transparency (default: 0.8)
    label : str, optional
        Label for legend
    **kwargs
        Additional arguments passed to scatter()
    """
    values = np.asarray(values)
    
    if values.ndim == 1:
        # Scalar field
        ax.scatter(
            gauss_coords,
            values,
            marker=marker,
            color=color,
            s=size,
            alpha=alpha,
            label=label,
            **kwargs
        )
    else:
        # Vector field: plot each component
        n_components = values.shape[1]
        colors = plt.cm.get_cmap('tab10')(np.linspace(0, 1, n_components))
        
        for comp in range(n_components):
            comp_label = f"{label} (comp {comp})" if label else f"Gauss (comp {comp})"
            ax.scatter(
                gauss_coords,
                values[:, comp],
                marker=marker,
                color=colors[comp],
                s=size,
                alpha=alpha,
                label=comp_label,
                **kwargs
            )


def plot_nodal_points(
    ax: matplotlib.axes.Axes,
    node_coords: np.ndarray,
    values: np.ndarray,
    marker: str = 's',
    color: str = 'blue',
    size: float = 80.0,
    alpha: float = 0.9,
    label: Optional[str] = None,
    edgecolors: Optional[str] = None,
    linewidths: Optional[float] = None,
    **kwargs
) -> None:
    """
    Plot nodal markers.
    
    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Matplotlib axes to plot on
    node_coords : np.ndarray
        Physical coordinates of nodes, shape (n_nodes, 3) or (n_nodes,)
        For 1D, uses first column (x-coordinate)
    values : np.ndarray
        Field values at nodes, shape (n_nodes,) or (n_nodes, n_components)
    marker : str, optional
        Marker style (default: 's' for square)
    color : str, optional
        Marker color (default: 'blue')
    size : float, optional
        Marker size (default: 80.0)
    alpha : float, optional
        Marker transparency (default: 0.9)
    label : str, optional
        Label for legend
    **kwargs
        Additional arguments passed to scatter()
    """
    # Extract x-coordinates for 1D plotting
    if node_coords.ndim > 1:
        x_coords = node_coords[:, 0]
    else:
        x_coords = node_coords
    
    values = np.asarray(values)
    
    if values.ndim == 1:
        # Scalar field
        scatter_kwargs = {
            'marker': marker,
            'color': color,
            's': size,
            'alpha': alpha,
            'label': label,
        }
        if edgecolors is not None:
            scatter_kwargs['edgecolors'] = edgecolors
        if linewidths is not None:
            scatter_kwargs['linewidths'] = linewidths
        scatter_kwargs.update(kwargs)
        ax.scatter(x_coords, values, **scatter_kwargs)
    else:
        # Vector field: plot each component
        n_components = values.shape[1]
        colors = plt.cm.get_cmap('tab10')(np.linspace(0, 1, n_components))
        
        for comp in range(n_components):
            comp_label = f"{label} (comp {comp})" if label else f"Node (comp {comp})"
            scatter_kwargs = {
                'marker': marker,
                'color': colors[comp],
                's': size,
                'alpha': alpha,
                'label': comp_label,
            }
            if edgecolors is not None:
                scatter_kwargs['edgecolors'] = edgecolors
            if linewidths is not None:
                scatter_kwargs['linewidths'] = linewidths
            scatter_kwargs.update(kwargs)
            ax.scatter(x_coords, values[:, comp], **scatter_kwargs)


def plot_interpolated_field(
    ax: matplotlib.axes.Axes,
    x_coords: np.ndarray,
    interpolated_values: np.ndarray,
    linestyle: str = '-',
    linewidth: float = 2.0,
    alpha: float = 0.7,
    label: Optional[str] = None,
    **kwargs
) -> None:
    """
    Plot continuous interpolated field.
    
    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Matplotlib axes to plot on
    x_coords : np.ndarray
        Physical coordinates for interpolation, shape (n_points,)
    interpolated_values : np.ndarray
        Interpolated field values, shape (n_points,) or (n_points, n_components)
    linestyle : str, optional
        Line style (default: '-')
    linewidth : float, optional
        Line width (default: 2.0)
    alpha : float, optional
        Line transparency (default: 0.7)
    label : str, optional
        Label for legend
    **kwargs
        Additional arguments passed to plot()
    """
    interpolated_values = np.asarray(interpolated_values)
    
    if interpolated_values.ndim == 1:
        # Scalar field
        ax.plot(
            x_coords,
            interpolated_values,
            linestyle=linestyle,
            linewidth=linewidth,
            alpha=alpha,
            label=label,
            **kwargs
        )
    else:
        # Vector field: plot each component
        n_components = interpolated_values.shape[1]
        colors = plt.cm.get_cmap('tab10')(np.linspace(0, 1, n_components))
        
        for comp in range(n_components):
            comp_label = f"{label} (comp {comp})" if label else f"Interpolated (comp {comp})"
            ax.plot(
                x_coords,
                interpolated_values[:, comp],
                linestyle=linestyle,
                linewidth=linewidth,
                alpha=alpha,
                color=colors[comp],
                label=comp_label,
                **kwargs
            )


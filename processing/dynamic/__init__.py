"""Deprecated compatibility shim for transient processing.

Prefer [`processing.transient`](../transient/__init__.py).
"""

__all__ = ["assemble_global_system", "apply_boundary_conditions", "newmark_integrate"]


def assemble_global_system(*args, **kwargs):
    from processing.transient.assembly import assemble_global_system as _impl

    return _impl(*args, **kwargs)


def apply_boundary_conditions(*args, **kwargs):
    from processing.transient.boundary_conditions import apply_boundary_conditions as _impl

    return _impl(*args, **kwargs)


def newmark_integrate(*args, **kwargs):
    from processing.transient.time_integration import newmark_integrate as _impl

    return _impl(*args, **kwargs)

__all__ = ["assemble_global_system", "apply_boundary_conditions", "newmark_integrate"]

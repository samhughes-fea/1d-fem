# pre_processing/element_library/nonlinear/large_rotations/corotational/utilities/rotation.py
"""Rotation-vector operators for corotational frames (Rodrigues / quaternion log)."""

from __future__ import annotations

import numpy as np


def skew(a: np.ndarray) -> np.ndarray:
    return np.array([[0.0, -a[2], a[1]], [a[2], 0.0, -a[0]], [-a[1], a[0], 0.0]])


def rotation_exp(theta: np.ndarray) -> np.ndarray:
    """Rotation matrix of the rotation vector ``theta`` (Rodrigues)."""
    t = float(np.linalg.norm(theta))
    S = skew(theta)
    if t < 1e-8:
        return np.eye(3) + S + 0.5 * S @ S
    return np.eye(3) + np.sin(t) / t * S + (1.0 - np.cos(t)) / (t * t) * S @ S


def rotation_log(R: np.ndarray) -> np.ndarray:
    """Rotation vector of ``R`` (through the quaternion, robust up to pi)."""
    tr = np.trace(R)
    w2 = np.array([1.0 + tr, 1.0 + 2 * R[0, 0] - tr, 1.0 + 2 * R[1, 1] - tr, 1.0 + 2 * R[2, 2] - tr])
    k = int(np.argmax(w2))
    q = np.zeros(4)
    q[k] = 0.5 * np.sqrt(w2[k])
    if k == 0:
        q[1] = (R[2, 1] - R[1, 2]) / (4 * q[0])
        q[2] = (R[0, 2] - R[2, 0]) / (4 * q[0])
        q[3] = (R[1, 0] - R[0, 1]) / (4 * q[0])
    elif k == 1:
        q[0] = (R[2, 1] - R[1, 2]) / (4 * q[1])
        q[2] = (R[0, 1] + R[1, 0]) / (4 * q[1])
        q[3] = (R[0, 2] + R[2, 0]) / (4 * q[1])
    elif k == 2:
        q[0] = (R[0, 2] - R[2, 0]) / (4 * q[2])
        q[1] = (R[0, 1] + R[1, 0]) / (4 * q[2])
        q[3] = (R[1, 2] + R[2, 1]) / (4 * q[2])
    else:
        q[0] = (R[1, 0] - R[0, 1]) / (4 * q[3])
        q[1] = (R[0, 2] + R[2, 0]) / (4 * q[3])
        q[2] = (R[1, 2] + R[2, 1]) / (4 * q[3])
    if q[0] < 0:
        q = -q
    v = q[1:]
    s = float(np.linalg.norm(v))
    if s < 1e-12:
        return 2.0 * v
    return 2.0 * np.arctan2(s, q[0]) / s * v


def exp_batch(th: np.ndarray) -> np.ndarray:
    """Rotation matrices (m, 3, 3) of rotation vectors (m, 3); complex-step safe."""
    t2 = np.sum(th * th, axis=1)
    small = np.abs(t2) < 1e-12
    t = np.sqrt(np.where(small, 1.0, t2))
    a = np.where(small, 1.0 - t2 / 6.0, np.sin(t) / t)
    b = np.where(small, 0.5 - t2 / 24.0, (1.0 - np.cos(t)) / (t * t))
    S = np.zeros((th.shape[0], 3, 3), dtype=th.dtype)
    S[:, 0, 1], S[:, 0, 2], S[:, 1, 2] = -th[:, 2], th[:, 1], -th[:, 0]
    S[:, 1, 0], S[:, 2, 0], S[:, 2, 1] = th[:, 2], -th[:, 1], th[:, 0]
    return np.eye(3) + a[:, None, None] * S + b[:, None, None] * (S @ S)


def log_batch(R: np.ndarray) -> np.ndarray:
    """Rotation vectors (m, 3) of rotations (m, 3, 3) below pi/2 in angle."""
    v = 0.5 * np.stack([R[:, 2, 1] - R[:, 1, 2], R[:, 0, 2] - R[:, 2, 0], R[:, 1, 0] - R[:, 0, 1]], axis=1)
    s2 = np.sum(v * v, axis=1)
    c = 0.5 * (np.trace(R, axis1=1, axis2=2) - 1.0)
    small = np.abs(s2) < 1e-24
    s = np.sqrt(np.where(small, 1.0, s2))
    f = np.where(small, 1.0 / c + s2 / 6.0, np.arctan(s / c) / s)
    return f[:, None] * v

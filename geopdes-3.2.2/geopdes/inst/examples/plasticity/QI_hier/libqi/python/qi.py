"""
qi.py — Python (ctypes) binding for libqi (THB-spline Quasi-Interpolant).

Build the shared library first:

    cd ..; make python

Then in Python::

    import numpy as np
    from qi import compute_coeffs, eval_tp, HSpace

    # Single-level tensor-product space, degree 3 in each direction
    hs = HSpace(degree=(3, 3), nel=(16, 16))
    x = np.random.rand(4000)
    y = np.random.rand(4000)
    f = (np.tanh(9*y - 9*x) + 1)/9 \
        + np.exp(-((10*x - 6)**2 + (10*y + 7)**2)) / 1.5
    coeffs = compute_coeffs(hs, x, y, f)

    # Evaluate on a grid
    gx, gy = np.meshgrid(np.linspace(0, 1, 100), np.linspace(0, 1, 100),
                         indexing='ij')
    fhat = eval_tp(hs, 0, coeffs, gx.ravel(), gy.ravel()).reshape(gx.shape)

For genuinely hierarchical spaces, pass the same dictionary structure as the
MATLAB MEX wrapper: ``HSpace(degree=(d, d), nel=[(nx0, ny0), (nx1, ny1), ...],
                              active=[idx_lev0_1based, idx_lev1_1based, ...])``
"""

from __future__ import annotations

import ctypes as ct
import os
import sys
from dataclasses import dataclass, field
from typing import Optional, Sequence, Union

import numpy as np


# --------------------------------------------------------------------------- #
# Locate the shared library
# --------------------------------------------------------------------------- #
def _find_libqi() -> ct.CDLL:
    here = os.path.dirname(os.path.abspath(__file__))
    candidates = [
        os.path.join(here, "..", "build", "libqi.so"),
        os.path.join(here, "..", "build", "libqi.dylib"),
        os.path.join(here, "libqi.so"),
        os.path.join(here, "libqi.dylib"),
        "libqi.so",
        "libqi.dylib",
    ]
    if sys.platform.startswith("win"):
        candidates = [
            os.path.join(here, "..", "build", "libqi.dll"),
            os.path.join(here, "libqi.dll"),
            "libqi.dll",
        ]
    for c in candidates:
        try:
            return ct.CDLL(c)
        except OSError:
            continue
    raise OSError(
        "libqi shared library not found.  Build it first:\n"
        "    cd .. && make"
    )


_lib = _find_libqi()


# --------------------------------------------------------------------------- #
# C structures
# --------------------------------------------------------------------------- #
class _CHSpace(ct.Structure):
    _fields_ = [
        ("degree_x",         ct.c_int),
        ("degree_y",         ct.c_int),
        ("nlevels",          ct.c_int),
        ("nel_x",            ct.POINTER(ct.c_int)),
        ("nel_y",            ct.POINTER(ct.c_int)),
        ("breaks_x_data",    ct.POINTER(ct.c_double)),
        ("breaks_x_offset",  ct.POINTER(ct.c_int)),
        ("breaks_y_data",    ct.POINTER(ct.c_double)),
        ("breaks_y_offset",  ct.POINTER(ct.c_int)),
        ("ndof",             ct.c_int),
        ("ndof_per_level",   ct.POINTER(ct.c_int)),
        ("active_offset",    ct.POINTER(ct.c_int)),
        ("active_indices",   ct.POINTER(ct.c_int)),
    ]


class _CData(ct.Structure):
    _fields_ = [
        ("npoints", ct.c_int),
        ("ncomp",   ct.c_int),
        ("x",       ct.POINTER(ct.c_double)),
        ("y",       ct.POINTER(ct.c_double)),
        ("f",       ct.POINTER(ct.c_double)),
    ]


class _COptions(ct.Structure):
    _fields_ = [
        ("max_degree",       ct.c_int),
        ("sigma_threshold",  ct.c_double),
        ("n_threads",        ct.c_int),
        ("verbose",          ct.c_int),
    ]


class _CDiag(ct.Structure):
    _fields_ = [
        ("deg_used",       ct.POINTER(ct.c_int)),
        ("nlocal",         ct.POINTER(ct.c_int)),
        ("cond_estimate",  ct.POINTER(ct.c_double)),
        ("level",          ct.POINTER(ct.c_int)),
    ]


# --------------------------------------------------------------------------- #
# Function prototypes
# --------------------------------------------------------------------------- #
_lib.qi_options_defaults.restype  = None
_lib.qi_options_defaults.argtypes = [ct.POINTER(_COptions)]

_lib.qi_compute_coeffs.restype  = ct.c_int
_lib.qi_compute_coeffs.argtypes = [ct.POINTER(_CHSpace), ct.POINTER(_CData),
                                   ct.POINTER(_COptions), ct.POINTER(ct.c_double),
                                   ct.POINTER(_CDiag)]

_lib.qi_eval_tp.restype  = ct.c_int
_lib.qi_eval_tp.argtypes = [ct.POINTER(_CHSpace), ct.c_int,
                            ct.POINTER(ct.c_double), ct.c_int, ct.c_int,
                            ct.POINTER(ct.c_double), ct.POINTER(ct.c_double),
                            ct.POINTER(ct.c_double)]

_lib.qi_to_finest.restype  = ct.c_int
_lib.qi_to_finest.argtypes = [ct.POINTER(_CHSpace), ct.POINTER(ct.c_double),
                              ct.c_int, ct.POINTER(ct.c_double)]

_lib.qi_version.restype = ct.c_char_p


# --------------------------------------------------------------------------- #
# Status -> exception
# --------------------------------------------------------------------------- #
_STATUS = {
    0: "OK",
    1: "INVALID_ARG",
    2: "OUT_OF_MEMORY",
    3: "NOT_ENOUGH_DATA",
    4: "NO_ADMISSIBLE_DEGREE",
    5: "SINGULAR_MATRIX",
}


class QIError(RuntimeError):
    pass


def _check(st: int) -> None:
    if st != 0:
        raise QIError(f"qi error: {_STATUS.get(st, st)} ({st})")


def version() -> str:
    return _lib.qi_version().decode()


# --------------------------------------------------------------------------- #
# Hierarchical-space description
# --------------------------------------------------------------------------- #
@dataclass
class HSpace:
    """Lightweight description of a (hierarchical) tensor-product B-spline space.

    Parameters
    ----------
    degree : (int, int)
        Polynomial degrees (d_x, d_y).
    nel : (int, int) | sequence of (int, int)
        Number of cells per level in (x, y).  Pass a single tuple for a
        single-level tensor-product space.
    active : sequence of array_like, optional
        Per-level 1-based active index arrays (matching MATLAB convention).
        If None, every basis function on every level is active.
    breaks : sequence of (array_like, array_like), optional
        Per-level breakpoints in (x, y).  If None, uniform breaks 0..1.
    """

    degree: tuple
    nel: Union[tuple, Sequence]
    active: Optional[Sequence[np.ndarray]] = None
    breaks: Optional[Sequence[tuple]] = None
    _kept_alive: list = field(default_factory=list, repr=False)

    def _as_c(self) -> _CHSpace:
        # normalise nel
        nel_list = self.nel
        if isinstance(nel_list, tuple) and len(nel_list) == 2 \
                and isinstance(nel_list[0], int):
            nel_list = [nel_list]
        nlevels = len(nel_list)

        nel_x = np.asarray([n[0] for n in nel_list], dtype=np.int32)
        nel_y = np.asarray([n[1] for n in nel_list], dtype=np.int32)

        ndof_per_level = np.empty(nlevels, dtype=np.int32)
        for lv in range(nlevels):
            ndof_per_level[lv] = (nel_x[lv] + self.degree[0]) \
                                  * (nel_y[lv] + self.degree[1])

        # active indices: 0-based, flat
        if self.active is None:
            active_offset = np.zeros(nlevels + 1, dtype=np.int32)
            for lv in range(nlevels):
                active_offset[lv+1] = active_offset[lv] + ndof_per_level[lv]
            ndof = int(active_offset[-1])
            active_idx = np.arange(ndof, dtype=np.int32)
        else:
            assert len(self.active) == nlevels
            offsets = [0]
            for a in self.active:
                offsets.append(offsets[-1] + len(a))
            active_offset = np.asarray(offsets, dtype=np.int32)
            ndof = int(active_offset[-1])
            active_idx = np.empty(ndof, dtype=np.int32)
            for lv, a in enumerate(self.active):
                active_idx[active_offset[lv]:active_offset[lv+1]] = \
                    np.asarray(a, dtype=np.int32) - 1   # 1-based -> 0-based

        # breaks
        bxd = byd = None
        bxo = byo = None
        if self.breaks is not None:
            bx_off = [0]; by_off = [0]
            bx_chunks = []; by_chunks = []
            for bx, by in self.breaks:
                bx = np.asarray(bx, dtype=np.float64)
                by = np.asarray(by, dtype=np.float64)
                bx_chunks.append(bx); by_chunks.append(by)
                bx_off.append(bx_off[-1] + len(bx))
                by_off.append(by_off[-1] + len(by))
            bxd = np.concatenate(bx_chunks).astype(np.float64)
            byd = np.concatenate(by_chunks).astype(np.float64)
            bxo = np.asarray(bx_off, dtype=np.int32)
            byo = np.asarray(by_off, dtype=np.int32)

        # keep numpy buffers alive while the struct is in use
        self._kept_alive = [nel_x, nel_y, ndof_per_level,
                            active_offset, active_idx,
                            bxd, byd, bxo, byo]

        H = _CHSpace()
        H.degree_x = int(self.degree[0])
        H.degree_y = int(self.degree[1])
        H.nlevels  = nlevels
        H.nel_x    = nel_x.ctypes.data_as(ct.POINTER(ct.c_int))
        H.nel_y    = nel_y.ctypes.data_as(ct.POINTER(ct.c_int))
        H.ndof     = ndof
        H.ndof_per_level = ndof_per_level.ctypes.data_as(ct.POINTER(ct.c_int))
        H.active_offset  = active_offset.ctypes.data_as(ct.POINTER(ct.c_int))
        H.active_indices = active_idx.ctypes.data_as(ct.POINTER(ct.c_int))
        if bxd is not None:
            H.breaks_x_data   = bxd.ctypes.data_as(ct.POINTER(ct.c_double))
            H.breaks_x_offset = bxo.ctypes.data_as(ct.POINTER(ct.c_int))
            H.breaks_y_data   = byd.ctypes.data_as(ct.POINTER(ct.c_double))
            H.breaks_y_offset = byo.ctypes.data_as(ct.POINTER(ct.c_int))
        else:
            H.breaks_x_data = None
            H.breaks_x_offset = None
            H.breaks_y_data = None
            H.breaks_y_offset = None
        return H

    @property
    def ndof(self) -> int:
        H = self._as_c()
        return H.ndof


# --------------------------------------------------------------------------- #
# Public API
# --------------------------------------------------------------------------- #
def compute_coeffs(hspace: HSpace,
                   x: np.ndarray,
                   y: np.ndarray,
                   f: np.ndarray,
                   *,
                   max_degree: int = 0,
                   sigma_threshold: float = 0.0,
                   n_threads: int = 0,
                   verbose: int = 0,
                   return_diag: bool = False) -> np.ndarray:
    """Compute the THB-spline QI coefficients."""
    x = np.ascontiguousarray(x, dtype=np.float64)
    y = np.ascontiguousarray(y, dtype=np.float64)
    f = np.ascontiguousarray(f, dtype=np.float64)
    if x.ndim != 1 or y.ndim != 1 or x.shape != y.shape:
        raise ValueError("x and y must be 1-D and same length")
    if f.ndim == 1:
        if f.shape[0] != x.shape[0]:
            raise ValueError("f.shape[0] must equal len(x)")
        ncomp = 1
        f_row = f.reshape(-1, 1)
    elif f.ndim == 2:
        if f.shape[0] != x.shape[0]:
            raise ValueError("f.shape[0] must equal len(x)")
        ncomp = f.shape[1]
        f_row = f
    else:
        raise ValueError("f must be 1-D or 2-D")
    f_row = np.ascontiguousarray(f_row, dtype=np.float64)

    H = hspace._as_c()
    D = _CData(npoints=x.size, ncomp=ncomp,
               x=x.ctypes.data_as(ct.POINTER(ct.c_double)),
               y=y.ctypes.data_as(ct.POINTER(ct.c_double)),
               f=f_row.ctypes.data_as(ct.POINTER(ct.c_double)))

    O = _COptions()
    _lib.qi_options_defaults(ct.byref(O))
    if max_degree:      O.max_degree = max_degree
    if sigma_threshold: O.sigma_threshold = sigma_threshold
    if n_threads:       O.n_threads = n_threads
    O.verbose = verbose

    coeffs = np.zeros(H.ndof * ncomp, dtype=np.float64)

    if return_diag:
        deg_used = np.zeros(H.ndof, dtype=np.int32)
        nlocal   = np.zeros(H.ndof, dtype=np.int32)
        cond_est = np.zeros(H.ndof, dtype=np.float64)
        level    = np.zeros(H.ndof, dtype=np.int32)
        diag = _CDiag(
            deg_used=deg_used.ctypes.data_as(ct.POINTER(ct.c_int)),
            nlocal=nlocal.ctypes.data_as(ct.POINTER(ct.c_int)),
            cond_estimate=cond_est.ctypes.data_as(ct.POINTER(ct.c_double)),
            level=level.ctypes.data_as(ct.POINTER(ct.c_int)),
        )
        diag_ptr = ct.byref(diag)
    else:
        diag_ptr = None

    st = _lib.qi_compute_coeffs(ct.byref(H), ct.byref(D), ct.byref(O),
                                coeffs.ctypes.data_as(ct.POINTER(ct.c_double)),
                                diag_ptr)
    _check(st)

    out = coeffs.reshape(H.ndof, ncomp)
    if ncomp == 1:
        out = out.ravel()
    if return_diag:
        return out, dict(deg_used=deg_used, nlocal=nlocal,
                         cond_estimate=cond_est, level=level)
    return out


def eval_tp(hspace: HSpace,
            level: int,
            coeffs_lev: np.ndarray,
            x: np.ndarray,
            y: np.ndarray) -> np.ndarray:
    """Evaluate a tensor-product B-spline at scattered points."""
    x = np.ascontiguousarray(x, dtype=np.float64)
    y = np.ascontiguousarray(y, dtype=np.float64)
    coeffs_lev = np.ascontiguousarray(coeffs_lev, dtype=np.float64)
    if coeffs_lev.ndim == 1:
        ncomp = 1
        cflat = coeffs_lev
    else:
        ncomp = coeffs_lev.shape[1]
        cflat = coeffs_lev.reshape(-1)
    H = hspace._as_c()
    out = np.zeros(x.size * ncomp, dtype=np.float64)
    st = _lib.qi_eval_tp(ct.byref(H), level,
                         cflat.ctypes.data_as(ct.POINTER(ct.c_double)),
                         ncomp, x.size,
                         x.ctypes.data_as(ct.POINTER(ct.c_double)),
                         y.ctypes.data_as(ct.POINTER(ct.c_double)),
                         out.ctypes.data_as(ct.POINTER(ct.c_double)))
    _check(st)
    if ncomp == 1:
        return out
    return out.reshape(x.size, ncomp)


__all__ = ["HSpace", "compute_coeffs", "eval_tp", "version", "QIError"]

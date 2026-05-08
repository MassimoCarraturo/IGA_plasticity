# libqi — Hierarchical THB-spline Quasi-Interpolant (C / OpenMP)

A C99 + OpenMP implementation of the hierarchical quasi-interpolant
(QI) for THB-splines described in

> C. Bracco, C. Giannelli, A. Sestini, *"Adaptive scattered data fitting by
> extension of local approximations to hierarchical splines"*, CAGD 52–53
> (2017).
>
> C. Bracco, C. Giannelli, F. Mazzia, A. Sestini, *"Adaptive fitting with
> THB-splines: Error analysis and industrial applications"*, CAGD 62 (2018),
> 239–252.

The library replaces the reference MATLAB implementation `get_QI_coeffs.m`
with a self-contained C library that can be linked from MATLAB (MEX),
Python (ctypes) and C++.

The outer loop over basis functions is parallelised with **OpenMP** —
each basis function is solved by an independent thread with its own
scratch space.

---

## Layout

```
qi_lib/
├── include/qi.h            # Public C API
├── src/
│   ├── qi.c                # Main algorithm (qi_compute_coeffs, qi_eval_tp, ...)
│   ├── qi_linalg.[ch]      # Pinv via Jacobi eig of A'A; LU; gemm
│   └── qi_bspline.[ch]     # findspan, basisfun, spcol  (NURBS Book Alg 2.1/2.2)
├── test/
│   ├── qi_test.c           # End-to-end correctness tests
│   └── qi_bench.c          # Parallel-scaling benchmark
├── matlab/
│   ├── qi_mex.c            # MATLAB MEX wrapper
│   └── get_QI_coeffs_c.m   # Drop-in replacement for get_QI_coeffs.m
├── python/
│   ├── qi.py               # ctypes binding (no external deps but numpy)
│   └── test_qi.py          # Python test driver
├── cpp/
│   ├── qi.hpp              # Header-only modern-C++ wrapper
│   └── qi_test.cpp         # C++ test driver
└── Makefile
```

---

## Build

### Linux / macOS

Requires gcc/clang with OpenMP and a POSIX `make`.

```bash
cd qi_lib
make                     # builds build/libqi.a, build/libqi.so (or .dylib),
                         # and build/qi_test
make test                # runs the C test suite
make bench               # parallel-scaling benchmark
```

### Windows

The C99 sources are portable. Use MSYS2/MinGW or WSL to build. With the
Microsoft compiler, replace `-fopenmp` with `/openmp` and `-fPIC` with
nothing.

---

## API summary (C)

```c
#include "qi.h"

qi_options_t opts; qi_options_defaults(&opts);

qi_status_t st = qi_compute_coeffs(&hspace,    /* qi_hspace_t  */
                                   &data,      /* qi_data_t    */
                                   &opts,
                                   coeffs,     /* [ndof*ncomp] row-major */
                                   &diag);     /* optional     */
```

### `qi_hspace_t`

Describes a (possibly hierarchical) tensor-product B-spline space:

* `degree_x`, `degree_y`
* `nlevels`
* `nel_x[ℓ]`, `nel_y[ℓ]` — number of elements per direction at level ℓ
* `breaks_x_data + breaks_x_offset[ℓ..ℓ+1]` — non-uniform breakpoints
  (set to `NULL` for the uniform default)
* `active_indices` — flat list of active basis-function indices, indexed
  by `active_offset[ℓ..ℓ+1]`. Indices are 0-based and use the column-major
  layout `idx = k_y * ndof_dir_x[ℓ] + k_x` (matches MATLAB
  `ind2sub([ndof_x, ndof_y], idx+1)`).

All arrays are passed by pointer; the library never takes ownership.

### `qi_data_t`

Plain scattered data: `npoints`, `ncomp`, `x[npoints]`, `y[npoints]`,
`f[npoints*ncomp]` (row-major).

### `qi_options_t`

* `max_degree` — overrides the default `d_x`
* `sigma_threshold` — degree-admissibility threshold σ (default 1e2,
  matching `soglia` in the reference MATLAB code).
  A degree p_d is accepted iff 1/σ_min(A_J) ≤ σ.
* `n_threads` — OpenMP thread count (0 = OpenMP default)
* `verbose` — 0 / 1 / 2

### Diagnostics

`qi_diag_t` exposes per-basis-function values: chosen polynomial degree,
number of local data points, the chosen `1/σ_min`, and the basis-function
level. Pass `NULL` for fields you don't need.

### Helpers

* `qi_eval_tp` — evaluate a tensor-product B-spline coefficient vector
  on a list of scattered points.
* `qi_to_finest` — embed THB coefficients into the finest tensor-product
  basis (untruncated). Useful for fast evaluation pipelines.

---

## MATLAB

```bash
cd qi_lib
make mex                  # or, from MATLAB:
                          % mex -R2018a -Iinclude -Isrc \
                          %    CFLAGS='$CFLAGS -fopenmp -O3' \
                          %    LDFLAGS='$LDFLAGS -fopenmp' \
                          %    matlab/qi_mex.c src/qi.c \
                          %    src/qi_linalg.c src/qi_bspline.c
```

Then in MATLAB:

```matlab
addpath('qi_lib/matlab')
coeffs = get_QI_coeffs_c(hspace, hmsh, data);
% identical signature/output as get_QI_coeffs.m, ~10-50x faster on
% multi-core CPUs.
```

The wrapper accepts the standard GeoPDEs `hspace` / `hmsh` structures.
Output `coeffs` is `[ndof × ncomp]` in the same basis-function ordering
as the MATLAB reference (active indices flattened across levels).

---

## Python

Build the shared library, then:

```python
import sys, numpy as np
sys.path.insert(0, 'qi_lib/python')
from qi import HSpace, compute_coeffs, eval_tp

hs = HSpace(degree=(3, 3), nel=(16, 16))
n = 4000
x = np.random.rand(n); y = np.random.rand(n)
f = (np.tanh(9*y - 9*x) + 1)/9 \
    + np.exp(-((10*x - 6)**2 + (10*y + 7)**2))/1.5
coeffs = compute_coeffs(hs, x, y, f)

gx, gy = np.meshgrid(np.linspace(0,1,100), np.linspace(0,1,100), indexing='ij')
fhat = eval_tp(hs, 0, coeffs, gx.ravel(), gy.ravel()).reshape(gx.shape)
```

For a hierarchical mesh, pass active indices as 1-based per-level lists:

```python
hs = HSpace(degree=(3, 3),
            nel=[(16, 16), (32, 32)],
            active=[active_lev0_1based, active_lev1_1based])
```

`compute_coeffs(..., return_diag=True)` returns the per-basis diagnostics
dictionary as well.

The wrapper depends only on `numpy` and the Python standard library
(ctypes).

---

## C++

```cpp
#include "qi.hpp"
#include <vector>

qi::HSpace hs({3, 3}, {{16, 16}});

std::vector<double> x, y, f;          // ... fill ...
auto coeffs = qi::compute_coeffs(hs, x, y, f);

std::vector<double> gx, gy;           // ... fill ...
auto fhat = qi::eval_tp(hs, 0, coeffs, gx, gy);
```

Build:

```bash
g++ -O2 -std=c++17 -fopenmp -Iinclude my_app.cpp build/libqi.a -lm -fopenmp
```

The wrapper is header-only (`cpp/qi.hpp`) and provides:

* `qi::HSpace`           — RAII wrapper that owns marshalled arrays
* `qi::Options`          — algorithm options
* `qi::compute_coeffs`   — main entry, returns `std::vector<double>`
* `qi::eval_tp`          — evaluate at scattered points
* `qi::to_finest`        — convert to fine-grid TP coefficients
* `qi::Error`            — exception raised on `QI_ERR_*`

---

## Algorithm notes

For each active basis function T^ℓ_J the code:

1. Locates its parametric support `[x_mu, x_nu] × [y_mu, y_nu]` from the
   level-ℓ break vectors.
2. Selects the local data set `F^ℓ_J` inside the disk centred at the
   support centre, with radius adaptively enlarged until at least
   `(p_d + 1)(p_d + 2)/2` points are gathered.
3. Tries polynomial degrees `p_d = d, d−1, …, 0` (descending). For each
   candidate degree it builds the local LS design matrix on the scaled
   coordinates `(u, v) ∈ [0,1]²`, computes the (Moore–Penrose) pseudoinverse
   via the Jacobi eigendecomposition of `A^T A`, and rejects the degree
   if `1/σ_min(A) > σ` (the threshold from `qi_options_t.sigma_threshold`).
4. Converts the chosen polynomial coefficients to the (x,y)-monomial
   basis and then, by collocation, to the local B-spline tensor product
   basis.
5. Extracts the row corresponding to (k_x, k_y).

Each basis function is independent — the outer `#pragma omp parallel for`
in `qi_compute_coeffs` distributes them across threads with dynamic
scheduling. Per-thread scratch buffers are reused across all basis
functions assigned to the same thread.

The serial (1-thread) numerical results match the MATLAB reference to
machine precision on polynomial test cases.

### Performance notes

* Parallel speedup is observable from a few hundred basis functions up.
  For very small problems (few dozen basis functions) the OpenMP overhead
  dominates.
* Scaling is dynamic because the per-basis cost varies (different local
  data set sizes, different chosen degrees). `schedule(dynamic, 16)` keeps
  load balanced.
* Memory usage is `O(n_threads × max_local_buffer)` — typically < 1 MB
  per thread for d ≤ 4.

---

## Testing

```bash
make test
# THB-spline QI test suite — version libqi 1.0.0  (Hierarchical Quasi-Interpolant; OpenMP)
# OpenMP enabled — max threads = …
#
# === Test 1: polynomial reproduction ===
#   max polynomial reproduction error = 3.952e-14   (machine precision)
#   PASS
#
# === Test 2: smooth function approximation ===
#   ndof = 361, npts = 4000, time = 0.03 s
#   degree histogram: d3=361
#   max approx error  = 1.088e-02
#   L2  approx error  = 3.918e-03
#   PASS
#
# === Test 3: vector-valued data (ncomp = 3) ===
#   component errors: 3.952e-14 3.952e-14 1.887e-14
#   PASS
```

The test suite exercises:

1. **Polynomial reproduction** — feeds an exact polynomial of degree ≤ d
   and checks the QI reproduces it to machine precision.
2. **Smooth approximation** — Bracco et al. (2018) Example 1 test
   function on a 16×16 mesh with d=3.
3. **Vector-valued data** — verifies that all components of a multi-valued
   function are processed correctly and independently.

The same tests are also runnable through the Python and C++ wrappers
(`python python/test_qi.py`, `./cpp/qi_test_cpp`).

The whole codebase is clean under AddressSanitizer + UndefinedBehaviorSanitizer
(both single-threaded and 4-thread runs).

---

## License

Implementation generated to match the algorithm described in the cited
papers and the reference MATLAB code. Use accordingly.

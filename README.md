# IGA Plasticity — Adaptive THB-Spline J2 Elastoplasticity

[![MATLAB CI](https://github.com/MassimoCarraturo/IGA_plasticity/actions/workflows/matlab-ci.yml/badge.svg?branch=development)](https://github.com/MassimoCarraturo/IGA_plasticity/actions/workflows/matlab-ci.yml)

MATLAB implementation of adaptive **J2 (von Mises) elastoplasticity** on Truncated Hierarchical B-spline (THB-spline) meshes, built as an extension of the [GeoPDEs](http://rafavzqz.github.io/geopdes/) isogeometric analysis library.

## Objectives

This project addresses two key challenges in isogeometric computational plasticity:

1. **Adaptive mesh refinement and coarsening for elastoplasticity** — the elastic-plastic front introduces sharp gradients that demand local mesh refinement, while regions that unload can benefit from coarsening. The framework implements a full SOLVE &rarr; PROJECT &rarr; ESTIMATE &rarr; MARK &rarr; REFINE/COARSEN adaptive loop on hierarchical THB-spline meshes.

2. **Accurate transfer of history variables across mesh changes** — in plasticity, internal variables (plastic strain, stress) must be projected onto the new mesh after refinement or coarsening without introducing spurious oscillations or violating yield constraints. The project implements and compares six projection operators, arranged as a target-space by fitting-rule design, together with two choices of transfer scope.

## Extensions from GeoPDEs

This project extends the [GeoPDEs 3.2.2](http://rafavzqz.github.io/geopdes/) framework in the following areas:

### Adaptive J2 Plasticity Solver
- **`adaptivity_J2_plasticity.m`** — Main adaptive solver implementing incremental load stepping with Newton-Raphson iteration, radial return mapping for J2 plasticity, and a full adaptive refinement/coarsening cycle at each load step.
- Support for **symmetry**, **slider** (frictionless contact), **Neumann pressure**, and **Dirichlet** boundary conditions.
- Admissible mesh refinement and coarsening strategies for THB-splines.

### History Variable Projection Methods
Three strategies for projecting plastic strain and stress onto refined/coarsened meshes:

Six operators, differing in the target space (maximal regularity `C^{p-1}` or a `C^0` space obtained by raising interior knots to multiplicity `p`) and in the fitting rule:

| Method | Target space | Description | Cost |
|--------|--------------|-------------|------|
| **L2** | `C^{p-1}` | Global mass-matrix projection (`M \ rhs`). Smooth but may oscillate near discontinuities. | Sparse solve per component |
| **DLSQ** | `C^{p-1}` | Global collocated least squares at the quadrature points, solved once by sparse QR (Hennig et al. 2018). | Sparse QR per component |
| **QI** | `C^{p-1}` | Local quasi-interpolant. On a hierarchy it uses the construction of Speleers and Manni (2016), with each functional supported on the region where its own level is finest. | O(N) per component |
| **QI_C0** | `C^0` | The same rule on the `C^0` space. Represents the kink at the elastic-plastic front exactly. | O(N) per component |
| **Bezier** | `C^{p-1}` | Element-local fit reconciled by a support-weighted average, through multi-level Bezier extraction (Thomas et al. 2015, D'Angella et al. 2018). | Dense solve per element |
| **Bezier_C0** | `C^0` | The same rule on the `C^0` space, where the average reduces to a single term. | Dense solve per element |

A seventh option, `DLSQ_W`, is the weighted local least squares of Hennig et al. (2018). At the analysis quadrature of `p+1` points per direction the element fit is square, so the measure cancels and it coincides with `Bezier`. It remains available for comparison.

### Scope of the History Transfer
History variables live at quadrature points and must be recomputed when the mesh changes. Two scopes are implemented, selected with `method_data.type_transfer`:

- **`PROJECT`** (default) — the transferred field is written to every active element.
- **`PROJECT_NEW`** — only the newly activated elements are filled, and elements surviving the refinement keep the data they already carried (Hennig et al. 2018, Eq. 40).

The distinction is immaterial on a `C^0` space, where an element carries its own functions, and consequential on `C^{p-1}`, where a basis function spans `p+1` elements. Two further options, **`CPT`** (closest point transfer) and **`WPLSQ`** (weighted patch least squares), transfer the quadrature data directly without passing through control values.

### C/OpenMP Backend (libqi)
A compiled C library with OpenMP parallelisation for the quasi-interpolant projections:
- **`qi.c`** — THB-spline quasi-interpolant: polynomial LS fit + basis change, with uniform-grid spatial indexing for fast point queries.
- **`qi_localls.c`** — Local least-squares B-spline projection: Tikhonov-penalised normal equations with LU pivoting and pseudoinverse fallback.
- **`qi_spatial.c`** — Counting-sort spatial hash for O(1) box queries.
- Compiled as MEX binaries via `build_mex.m`; automatic fallback to pure MATLAB if MEX is unavailable.

### Error Estimators
Two element-level error indicators to drive adaptive refinement:

- **Equilibrium residual** (`div_sigma`) — measures local violation of equilibrium: `eta_K = h_K * ||f + div(sigma_h)||`.
- **Stress gradient** (`stress_gradient`) — stress-gradient seminorm: `eta_K = h_K * sqrt(int_K ||grad sigma||_F^2)`. Concentrates refinement at elastic-plastic fronts.

### Hierarchical Mesh Extensions
- Extended `hierarchical_mesh` and `hierarchical_space` classes with support for coarsening and multi-level state variable management.
- Level-wise evaluation and projection utilities for history variables.

## Project Structure

```
IGA_plasticity/
├── geopdes_hierarchical_plasticity/     # Main project code
│   ├── adaptivity_iga/                  # Adaptive solver + refine/coarsen
│   ├── examples/plasticity/             # Benchmarks and examples
│   ├── hierarchical_classes/            # THB-spline mesh/space classes
│   ├── initialize/                      # Problem initialisation
│   └── quasi_interpolant_hierarchical/  # QI projection + libqi C backend
│       └── libqi/
│           ├── src/                     # C source (qi.c, qi_localls.c, ...)
│           ├── matlab/                  # MEX binaries + MATLAB wrappers
│           └── test/                    # Unit tests
├── geopdes-3.2.2/                       # GeoPDEs base library
├── nurbs-1.4.3/                         # NURBS toolbox
├── ci/                                  # CI test runner
└── .github/workflows/                   # GitHub Actions CI pipeline
```

## Examples

### 3D: Thick-Walled Sphere Under Internal Pressure
**`examples/plasticity/ex_plastic_sphere_hier.m`**

Eighth-sphere benchmark (de Souza Neto, Peric & Owen) with steel material properties (E = 210 GPa, nu = 0.3, sigma_y = 240 MPa). Five load steps to ~99% of the analytical limit pressure. Results are validated against **Hill's analytical solution** for radial and tangential stress distributions.

### 2D: Plane-Strain Thick-Walled Ring
**`examples/plasticity/ex_plane_strain_hill_plasticity_hier.m`**

Quarter-cylinder plane-strain problem with the same material, 10 load steps. Symmetry boundary conditions and Neumann pressure loading. Validated against Hill's analytical solution.

### Projector Regression Tests
Fast checks that need no plasticity, no Newton iteration and no load loop. Shared mesh and space construction lives in **`projector_testbed.m`**.

- **`test_projector_reproduction.m`** — projects a spline that already lies in the space and checks it comes back, which is the defining property of a projector.
- **`test_projector_convergence.m`** — one projection of a smooth function on a uniform family, checking the optimal rate `p+1`.
- **`test_projector_stability.m`** — repeated application on a refine/coarsen cycle, checking the operator does not drift.
- **`test_projector_kink.m`** — one projection of a kinked field, reporting each operator's error relative to the best approximation in its own target space. This is the test that separates the operators.
- **`test_bezier_vs_dlsqw.m`** — measures the gap between `Bezier` and `DLSQ_W` against degree and quadrature order.

### Comparative Studies
- **`study_front_metrics_cylinder.m`** / **`study_front_metrics_sphere.m`** — the operator comparison on identical, solution-independent meshes. Degree, Poisson ratio, operator list, refinement depth and transfer scope are overridable through the `FM_DEGREE`, `FM_NU`, `FM_METHODS`, `FM_LEVELS` and `FM_TRANSFER` environment variables.
- **`benchmark_qi_optimized.m`** — L2 vs QI vs QI_C0 accuracy and timing comparison.
- **`compare_estimators.m`** — Full factorial study: 2 estimators x 3 projection methods.
- **`adaptivity_max_level_study.m`** / **`adaptivity_max_level_study_2D.m`** — Convergence studies with increasing refinement levels.

## Requirements

- **MATLAB** R2018a or later
- **C compiler** configured for MEX (`mex -setup C`):
  - Windows: Microsoft Visual C++ (MSVC)
  - Linux/macOS: GCC or Clang with OpenMP support
- Dependencies (included in the repository):
  - [GeoPDEs 3.2.2](http://rafavzqz.github.io/geopdes/)
  - [NURBS Toolbox 1.4.3](https://octave.sourceforge.io/nurbs/)

## Quick Start

```matlab
% 1. Add paths
addpath(genpath('nurbs-1.4.3/nurbs-1.4.3/inst'));
addpath(genpath('geopdes-3.2.2/geopdes/inst'));
addpath(genpath('geopdes_hierarchical_plasticity'));

% 2. Build MEX binaries (once)
cd geopdes_hierarchical_plasticity/quasi_interpolant_hierarchical/libqi
build_mex();
cd ../../..

% 3. Run the sphere benchmark
cd geopdes_hierarchical_plasticity/examples/plasticity
ex_plastic_sphere_hier
```

## CI/CD

The project uses GitHub Actions with the [MATLAB Actions](https://github.com/matlab-actions/setup-matlab) toolchain. On every push to `main` or `development`:

1. **Build** — compiles MEX binaries from C source with GCC + OpenMP
2. **Unit tests** — verifies MEX output matches MATLAB reference (6 checks)
3. **Smoke test** — runs a full adaptive plasticity solve on a small mesh

Run tests locally:
```matlab
cd ci
run_ci_tests           % build + unit + smoke
run_ci_tests('unit')   % unit tests only
run_ci_tests('all')    % everything including benchmark
```

## References

- R. Vazquez, *"A new design for the implementation of isogeometric analysis in Octave and Matlab: GeoPDEs 3.0"*, Computers & Mathematics with Applications, 2016.
- E.A. de Souza Neto, D. Peric, D.R.J. Owen, *"Computational Methods for Plasticity: Theory and Applications"*, Wiley, 2008.
- R. Hill, *"The Mathematical Theory of Plasticity"*, Oxford University Press, 1950.
- C. Giannelli, B. Juettler, H. Speleers, *"THB-splines: The truncated basis for hierarchical splines"*, CAGD, 2012.

## License

This project is provided for academic and research purposes.

/*
 * qi.h — Public C API for the Hierarchical Quasi-Interpolant (THB-spline QI)
 *
 * C99 / OpenMP implementation of the algorithm described in
 *   Bracco, Giannelli, Sestini, "Adaptive scattered data fitting by extension
 *   of local approximations to hierarchical splines", CAGD 52-53 (2017),
 *   and Bracco et al., "Adaptive fitting with THB-splines: Error analysis and
 *   industrial applications", CAGD 62 (2018), pp. 239-252.
 *
 * The algorithm computes one coefficient λ^ℓ_J ∈ R^m per active hierarchical
 * basis function (the THB-spline T^ℓ_J).  Each computation is independent:
 *
 *   1. Locate the parametric support box [x_mu,x_nu]×[y_mu,y_nu] of the mother
 *      B-spline B^ℓ_J at level ℓ.
 *   2. Select a local data set F^ℓ_J inside the disk centred at the support
 *      centre, with radius adaptively enlarged until enough data is gathered.
 *   3. Try polynomial degrees p_d = d, d-1, ..., 0 until the local LS design
 *      matrix A is well-conditioned (σ_min(A) > 1/σ_thresh).
 *   4. Solve the local polynomial least-squares problem; convert to the local
 *      B-spline basis by collocation; extract the row corresponding to (k_x, k_y).
 *
 * Outer loop is parallelised with OpenMP — each thread treats one basis
 * function with its own scratch space.
 *
 * The interface is plain C99 with primitive arrays, so it can be called
 * verbatim from C, C++, MATLAB (MEX) and Python (ctypes / cffi / cython).
 *
 * Author: Generated for the project — fully compatible with the reference
 *         MATLAB implementation get_QI_coeffs.m
 */

#ifndef QI_H
#define QI_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ---------- Versioning ---------------------------------------------------- */
#define QI_VERSION_MAJOR 1
#define QI_VERSION_MINOR 0
#define QI_VERSION_PATCH 0

/* ---------- Status codes -------------------------------------------------- */
typedef enum {
    QI_OK                     = 0,
    QI_ERR_INVALID_ARG        = 1,
    QI_ERR_OUT_OF_MEMORY      = 2,
    QI_ERR_NOT_ENOUGH_DATA    = 3,   /* fewer local points than (1) basis  */
    QI_ERR_NO_ADMISSIBLE_DEG  = 4,   /* every degree failed the σ-check    */
    QI_ERR_SINGULAR_MATRIX    = 5
} qi_status_t;

/* ---------- Hierarchical space description -------------------------------- */
/*
 * The hierarchical mesh / space is described by primitive arrays so it can
 * be marshalled from MATLAB / Python / C++ without any translation layer.
 *
 * All indices passed in are 0-based (subtract 1 from MATLAB indices).
 *
 * Layout conventions:
 *   - Tensor product index of basis function (k_x, k_y) at level ℓ is the
 *     column-major linear index   idx = k_y * ndof_dir_x[ℓ] + k_x
 *     where ndof_dir_x[ℓ] = nel_x[ℓ] + degree_x and similarly for y.
 *     This matches GeoPDEs/MATLAB ind2sub([ndof_dir_x ndof_dir_y], idx+1).
 *
 *   - Knot/break vectors are stored flat with offsets:
 *         breaks_x_data[breaks_x_offset[ℓ] .. breaks_x_offset[ℓ+1]-1]
 *     gives the (nel_x[ℓ]+1) breakpoints at level ℓ.
 *     breaks_x_offset has length (nlevels+1), with breaks_x_offset[0]=0.
 *
 *   - Active basis indices at level ℓ live in
 *         active_indices[active_offset[ℓ] .. active_offset[ℓ+1]-1]
 *     each entry is the column-major idx (0-based) into the level-ℓ tensor
 *     product space described above.
 */
typedef struct {
    int    degree_x;            /* d_x (e.g. 2, 3, 4)                       */
    int    degree_y;            /* d_y                                      */
    int    nlevels;             /* M, number of hierarchical levels         */

    /* Per-level mesh */
    const int    *nel_x;        /* [nlevels] cells in x at each level       */
    const int    *nel_y;        /* [nlevels] cells in y at each level       */

    /* Optional non-uniform breakpoints. If breaks_x_data == NULL the code
       falls back to uniform breaks 0, 1/nel, 2/nel, ..., 1 at each level. */
    const double *breaks_x_data;
    const int    *breaks_x_offset;   /* [nlevels+1] or NULL                 */
    const double *breaks_y_data;
    const int    *breaks_y_offset;   /* [nlevels+1] or NULL                 */

    /* Active basis */
    int           ndof;              /* total number of THB-splines         */
    const int    *ndof_per_level;    /* [nlevels]                           */
    const int    *active_offset;     /* [nlevels+1]                         */
    const int    *active_indices;    /* [ndof]                              */
} qi_hspace_t;

/* ---------- Scattered data ------------------------------------------------ */
typedef struct {
    int           npoints;       /* n                                       */
    int           ncomp;         /* m, vector dimension of f_i              */
    const double *x;             /* [npoints]   x-coordinates in [0,1]      */
    const double *y;             /* [npoints]   y-coordinates in [0,1]      */
    const double *f;             /* [npoints*ncomp] row-major:
                                    f[i*ncomp + k] is the k-th component
                                    of f_i.                                 */
} qi_data_t;

/* ---------- Algorithm options -------------------------------------------- */
typedef struct {
    int    max_degree;           /* default: degree_x (= d). Pass 0 to use  *
                                    the default.                            */
    double sigma_threshold;      /* σ in the paper (e.g. 1e2 .. 1e8).
                                    A degree p_d is admissible iff
                                    1/σ_min(A_J) < sigma_threshold.
                                    Pass 0.0 for default (1e2 — same as
                                    "soglia" in the reference MATLAB code). */
    int    n_threads;            /* OpenMP threads. 0 = OpenMP default.     */
    int    verbose;              /* 0 = quiet, 1 = progress, 2 = debug      */
} qi_options_t;

/* fill an options struct with defaults */
void qi_options_defaults(qi_options_t *opts);

/* ---------- Diagnostics returned by qi_compute_coeffs --------------------- */
typedef struct {
    /* per-basis-function output diagnostics; arrays have length ndof.
       Pass NULL pointers if you don't want them.                           */
    int    *deg_used;            /* polynomial degree actually used         */
    int    *nlocal;              /* number of local data points used        */
    double *cond_estimate;       /* 1/σ_min(A) chosen at that degree        */
    int    *level;               /* level ℓ of each basis function          */
} qi_diag_t;

/* ---------- Main entry point --------------------------------------------- */
/*
 *  Compute the QI coefficients λ^ℓ_J for every active basis function.
 *
 *  Output:
 *    coeffs      -- buffer of length ndof * ncomp, row-major:
 *                   coeffs[k*ncomp + c] is the c-th component of the
 *                   coefficient associated with the k-th THB-spline.
 *                   The ordering of basis functions matches the order of
 *                   active_indices flattened over levels (same as the
 *                   MATLAB reference).
 *    diag        -- optional diagnostics (any field may be NULL).
 *
 *  Returns QI_OK on success, otherwise an error code.
 */
qi_status_t qi_compute_coeffs(const qi_hspace_t  *hspace,
                              const qi_data_t    *data,
                              const qi_options_t *opts,
                              double             *coeffs,
                              qi_diag_t          *diag);

/* ---------- Evaluation of the QI on scattered points --------------------- */
/*
 *  Helper: evaluate the THB-spline approximation at a list of points.
 *  This mirrors sp_eval_alt.m but uses the basis-function level and index
 *  to select the correct mother B-spline.
 *
 *  Implementation note: this evaluates each basis function on the level it
 *  was *defined* — without applying the truncation operator — which is
 *  enough for the "no truncation needed in single-level" case used by
 *  the test driver.  For multi-level meshes, you should pass the result
 *  through the standard hspace_subdivision_matrix pipeline (call
 *  qi_compute_coeffs and use the existing MATLAB / Python toolchain).
 *
 *  This function is provided mostly so that the C-only test driver can
 *  validate the implementation against synthetic data on a single-level
 *  tensor-product space.
 */
qi_status_t qi_eval_tp(const qi_hspace_t *hspace,
                       int                level,        /* 0-based  */
                       const double      *coeffs_lev,   /*[ndof_lev*ncomp]*/
                       int                ncomp,
                       int                npts,
                       const double      *x,
                       const double      *y,
                       double            *out);         /*[npts*ncomp] */

/* ---------- Utility: subdivision (fine-level conversion) ----------------- *
 *
 *  Convert a list of THB-spline coefficients λ^ℓ_J into the equivalent
 *  list of coefficients on the finest tensor-product B-spline basis
 *  B^{M-1}.  This corresponds to the matrix C{nlevels} returned by
 *  hspace_subdivision_matrix(...) in GeoPDEs.
 *
 *  Output buffer must have length ndof_finest * ncomp where
 *  ndof_finest = (nel_x[M-1]+d_x) * (nel_y[M-1]+d_y).
 *
 *  Note: we apply only the "untruncated" hierarchical basis embedding
 *  (preservation of coefficients), which is sufficient for evaluation
 *  on the finest tensor-product grid.
 */
qi_status_t qi_to_finest(const qi_hspace_t *hspace,
                         const double      *coeffs,    /*[ndof*ncomp]*/
                         int                ncomp,
                         double            *coeffs_finest);

/* ---------- Library version string --------------------------------------- */
const char *qi_version(void);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* QI_H */

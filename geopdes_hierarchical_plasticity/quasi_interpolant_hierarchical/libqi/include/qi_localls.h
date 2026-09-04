/*
 * qi_localls.h — Hierarchical local-LS B-spline projection (parallel C/OpenMP).
 *
 * For each active hierarchical THB-spline T^lev_kl, this routine solves a small
 * penalised normal-equations system over the local set of data points whose
 * coordinates fall in T^lev_kl's parametric support box, and returns the row
 * of the solution corresponding to kl.  This is a faithful port of
 * `getcoeff_localLS_Bspl.m` from the GeoPDEs hierarchical-plasticity codebase.
 *
 * The two heavy MATLAB-side ingredients (the level-wise gradgrad mass matrix
 * and the level-wise collocation matrix at the data points) are computed
 * once per call by the MATLAB wrapper using the standard GeoPDEs operators
 * (op_gradgradu_gradgradv_tp, basisfun_multi) and passed to libqi as CSC
 * sparse arrays.  The C side handles the per-DOF support / index lookup
 * and the parallel local solve loop.
 *
 * Supports parametric dimensions 2 and 3.  Outputs match the MATLAB
 * reference to numerical precision.
 */
#ifndef QI_LOCALLS_H
#define QI_LOCALLS_H

#include <stddef.h>
#include "qi.h"   /* for qi_status_t */

#ifdef __cplusplus
extern "C" {
#endif

/* ---------- Hierarchical space (parametric-dim agnostic) ----------------- */
/*
 * Indexing conventions (match MATLAB GeoPDEs after subtracting 1 from indices):
 *   - active_indices[k] is the column-major linear index of the k-th active
 *     basis function within its level's tensor-product space.
 *     For dim = 3:  idx = i + j*ndof_dir[0] + k*ndof_dir[0]*ndof_dir[1]
 *     (MATLAB ind2sub([ndof_dir(1) ndof_dir(2) ndof_dir(3)], idx+1) returns
 *      [i+1, j+1, k+1].)
 *   - breaks_data + breaks_offset[lev*par_dim + d ..] gives nel_dir[lev,d]+1
 *     non-uniform breakpoints in direction d at level lev.
 */
typedef struct {
    int           par_dim;         /* 2 or 3                                 */
    int           nlevels;
    const int    *degree;          /* [par_dim]                              */
    const int    *nel_dir;         /* [nlevels * par_dim] row-major          */
    const int    *ndof_dir;        /* [nlevels * par_dim] row-major          */
    const double *breaks_data;     /* flat                                    */
    const int    *breaks_offset;   /* [nlevels * par_dim + 1]                */

    /* Per-direction-per-level cell connectivity:
     * conn_first_data[conn_first_offset[lv*par_dim + d] + c] is the 0-based
     * index of the first 1D basis function active on cell c.  Each cell has
     * exactly (degree[d] + 1) active 1D basis functions, conn_first[c] ..
     * conn_first[c] + degree[d].  This handles arbitrary regularity (in
     * particular regularity = 0 used by the QI_ref projection space). */
    const int    *conn_first_data;
    const int    *conn_first_offset; /* [nlevels * par_dim + 1] */

    int           ndof;            /* total active hierarchical DOFs          */
    const int    *ndof_per_level;  /* [nlevels]                              */
    const int    *active_offset;   /* [nlevels + 1]                          */
    const int    *active_indices;  /* [ndof] 0-based level-local TP indices  */
} qi_ls_hspace_t;

/* ---------- Scattered data ---------------------------------------------- */
typedef struct {
    int           npoints;
    int           ncomp;          /* vector dimension of f                  */
    int           par_dim;
    const double *coords;         /* [npoints * par_dim] row-major          */
    const double *f;              /* [npoints * ncomp]   row-major          */
    const double *weight;         /* [npoints] or NULL for unit weights     */
} qi_ls_data_t;

/* ---------- CSC sparse matrix ------------------------------------------- */
/*
 * MATLAB's sparse matrices are CSC.  This struct mirrors that layout:
 *   col_ptr has length (ncols + 1)
 *   row_ind has length nnz = col_ptr[ncols]
 *   values  has length nnz
 * All indices are 0-based.
 */
typedef struct {
    int           nrows;
    int           ncols;
    const int    *col_ptr;
    const int    *row_ind;
    const double *values;
} qi_csc_t;

/* Per-level pre-computed matrices: gradgrad mass and collocation (at data). */
typedef struct {
    qi_csc_t      mass;           /* nrows = ncols = ndof_lev               */
    qi_csc_t      col;            /* nrows = ndof_lev, ncols = npoints       */
} qi_ls_level_precomp_t;

/* ---------- Options ----------------------------------------------------- */
typedef struct {
    double        lambda;         /* Tikhonov coefficient on mass matrix     */
    int           n_threads;      /* OpenMP threads (0 = OpenMP default)     */
    int           verbose;        /* 0/1                                     */
} qi_ls_options_t;

/* Initialise an options struct with defaults (lambda=0, OMP default). */
void qi_ls_options_defaults(qi_ls_options_t *opts);

/* ---------- Main entry point -------------------------------------------- */
/*
 *   coeffs:  output buffer of length ndof * ncomp, row-major
 *            coeffs[k*ncomp + c] is the c-th component of the
 *            coefficient associated with the k-th hierarchical basis
 *            function (matching MATLAB's QI_coeff(k+1, c+1)).
 *
 * Returns QI_OK on success.
 */
qi_status_t qi_localLS_compute(const qi_ls_hspace_t        *hs,
                               const qi_ls_data_t          *dat,
                               const qi_ls_level_precomp_t *level_pre,
                               const qi_ls_options_t       *opts,
                               double                      *coeffs);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* QI_LOCALLS_H */

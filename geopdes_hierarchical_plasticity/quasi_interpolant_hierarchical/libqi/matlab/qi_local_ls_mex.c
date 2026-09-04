/*
 * qi_local_ls_mex.c — MATLAB MEX wrapper around qi_localLS_compute.
 *
 * Usage:
 *   coeffs = qi_local_ls_mex(hspace_info, data_coords, fvals, lambda, ...
 *                            weight_or_empty, mass_cell, col_cell, opts_or_empty);
 *
 * Inputs:
 *   hspace_info      struct with fields:
 *     .par_dim         scalar (2 or 3)
 *     .degree          [1 x par_dim] integer
 *     .nlevels         scalar
 *     .ndof            scalar (total active hierarchical DOFs)
 *     .ndof_per_level  [nlevels x 1] or [1 x nlevels]
 *     .active_indices  [ndof x 1] 1-based level-local TP linear indices
 *     .active_offset   [nlevels+1 x 1] offsets into active_indices (0-based)
 *     .ndof_dir        [nlevels x par_dim] (row lev = ndof_dir at level lev)
 *     .nel_dir         [nlevels x par_dim]
 *     .breaks          {nlevels x 1} cell, breaks{lv} = {par_dim x 1} cell of breakpoints
 *
 *   data_coords      [npoints x par_dim] doubles
 *   fvals            [npoints x ncomp]   doubles
 *   lambda           scalar double
 *   weight_or_empty  [npoints x 1] doubles, or [] for unit weights
 *   mass_cell        {nlevels x 1} cell, each entry sparse double of size ndof_lev x ndof_lev
 *                    (use [] for levels with ndof_per_level == 0)
 *   col_cell         {nlevels x 1} cell, each entry sparse double of size ndof_lev x npoints
 *   opts_or_empty    struct with optional fields .n_threads, .verbose
 *
 * Output:
 *   coeffs           [ndof x ncomp] doubles
 */

#include "mex.h"
#include "matrix.h"
#include "qi.h"
#include "qi_localls.h"
#include <string.h>
#include <stdlib.h>

/* ---------- Field readers ------------------------------------------------ */
static int field_scalar_int(const mxArray *s, const char *name, int dflt)
{
    mxArray *f = mxGetField(s, 0, name);
    if (!f || mxIsEmpty(f)) return dflt;
    return (int)mxGetScalar(f);
}

static double field_scalar_double(const mxArray *s, const char *name, double dflt)
{
    mxArray *f = mxGetField(s, 0, name);
    if (!f || mxIsEmpty(f)) return dflt;
    return mxGetScalar(f);
}

/* Copy a double-or-int32 vector field to a freshly-allocated int array.
 * subtract_one: 1 = treat input as 1-based and convert to 0-based. */
static int *field_int_array(const mxArray *s, const char *name,
                            int *n_out, int subtract_one)
{
    mxArray *f = mxGetField(s, 0, name);
    if (!f || mxIsEmpty(f)) { *n_out = 0; return NULL; }
    size_t n = mxGetNumberOfElements(f);
    int *out = (int*)mxMalloc(n * sizeof(int));
    if (mxIsDouble(f)) {
        double *p = mxGetPr(f);
        for (size_t i = 0; i < n; ++i) out[i] = (int)p[i] - (subtract_one ? 1 : 0);
    } else if (mxIsInt32(f)) {
        int *p = (int*)mxGetData(f);
        for (size_t i = 0; i < n; ++i) out[i] = p[i] - (subtract_one ? 1 : 0);
    } else {
        mxFree(out);
        mexErrMsgIdAndTxt("qi_ls:badtype", "field %s must be double or int32", name);
    }
    *n_out = (int)n;
    return out;
}

/* Read a [nrows x ncols] numeric matrix from a struct field, copying to a new
 * row-major int array (dim2 = ncols).  MATLAB stores column-major, so we
 * transpose on the fly. */
static int *field_matrix_int_rowmajor(const mxArray *s, const char *name,
                                      int nrows, int ncols, int subtract_one)
{
    mxArray *f = mxGetField(s, 0, name);
    if (!f || mxIsEmpty(f))
        mexErrMsgIdAndTxt("qi_ls:missing", "field %s is required", name);
    if ((int)mxGetM(f) != nrows || (int)mxGetN(f) != ncols)
        mexErrMsgIdAndTxt("qi_ls:badsize", "field %s must be %dx%d", name, nrows, ncols);

    int *out = (int*)mxMalloc((size_t)nrows * ncols * sizeof(int));
    if (mxIsDouble(f)) {
        double *p = mxGetPr(f);  /* column-major */
        for (int i = 0; i < nrows; ++i)
            for (int j = 0; j < ncols; ++j)
                out[i * ncols + j] = (int)p[i + j * nrows] - (subtract_one ? 1 : 0);
    } else if (mxIsInt32(f)) {
        int *p = (int*)mxGetData(f);
        for (int i = 0; i < nrows; ++i)
            for (int j = 0; j < ncols; ++j)
                out[i * ncols + j] = p[i + j * nrows] - (subtract_one ? 1 : 0);
    } else {
        mxFree(out);
        mexErrMsgIdAndTxt("qi_ls:badtype", "field %s must be double or int32", name);
    }
    return out;
}

/* Flatten breaks{lv}{d} (cell of cell of vectors) into a single double array
 * with associated offsets [nlevels*par_dim + 1]. */
static void flatten_breaks(const mxArray *brk, int nlevels, int par_dim,
                           double **data_out, int **offset_out)
{
    int *off = (int*)mxMalloc((size_t)(nlevels * par_dim + 1) * sizeof(int));
    off[0] = 0;
    /* First pass: total length */
    for (int lv = 0; lv < nlevels; ++lv) {
        mxArray *clv = mxGetCell(brk, lv);
        if (!clv || !mxIsCell(clv) || (int)mxGetNumberOfElements(clv) < par_dim)
            mexErrMsgIdAndTxt("qi_ls:breaks", "breaks{%d} must be a 1x%d cell", lv+1, par_dim);
        for (int d = 0; d < par_dim; ++d) {
            mxArray *bd = mxGetCell(clv, d);
            if (!bd || !mxIsDouble(bd))
                mexErrMsgIdAndTxt("qi_ls:breaks", "breaks{%d}{%d} must be double", lv+1, d+1);
            off[lv * par_dim + d + 1] = off[lv * par_dim + d] + (int)mxGetNumberOfElements(bd);
        }
    }
    int total = off[nlevels * par_dim];
    double *dat = (double*)mxMalloc((size_t)(total > 0 ? total : 1) * sizeof(double));
    /* Second pass: copy */
    for (int lv = 0; lv < nlevels; ++lv) {
        mxArray *clv = mxGetCell(brk, lv);
        for (int d = 0; d < par_dim; ++d) {
            mxArray *bd = mxGetCell(clv, d);
            int n = (int)mxGetNumberOfElements(bd);
            memcpy(dat + off[lv * par_dim + d], mxGetPr(bd), (size_t)n * sizeof(double));
        }
    }
    *data_out   = dat;
    *offset_out = off;
}

/* Convert a MATLAB sparse double matrix to a qi_csc_t (with int-typed indices).
 * MATLAB's mxGetIr/mxGetJc return mwIndex (size_t).  We allocate freshly-typed
 * int arrays via mxMalloc for compatibility with libqi. */
static void copy_sparse_to_csc(const mxArray *m, qi_csc_t *out,
                               int **alloc_col_ptr, int **alloc_row_ind)
{
    if (!m || !mxIsSparse(m) || !mxIsDouble(m))
        mexErrMsgIdAndTxt("qi_ls:sparse", "expected a sparse double matrix");
    out->nrows = (int)mxGetM(m);
    out->ncols = (int)mxGetN(m);
    mwIndex *jc = mxGetJc(m);
    mwIndex *ir = mxGetIr(m);
    double  *pr = mxGetPr(m);
    int nnz = (int)jc[out->ncols];

    int *cp = (int*)mxMalloc((size_t)(out->ncols + 1) * sizeof(int));
    int *ri = (int*)mxMalloc((size_t)(nnz > 0 ? nnz : 1) * sizeof(int));
    for (int j = 0; j <= out->ncols; ++j) cp[j] = (int)jc[j];
    for (int q = 0; q < nnz; ++q)         ri[q] = (int)ir[q];

    out->col_ptr = cp;
    out->row_ind = ri;
    out->values  = pr;            /* shared with mxArray, no copy */
    *alloc_col_ptr = cp;
    *alloc_row_ind = ri;
}

/* ------------------------------------------------------------------------ */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs < 7)
        mexErrMsgIdAndTxt("qi_ls:nrhs",
            "Usage: coeffs = qi_local_ls_mex(hspace, data, f, lambda, weight, mass_cell, col_cell [, opts])");

    const mxArray *hsp_in     = prhs[0];
    const mxArray *coords_in  = prhs[1];
    const mxArray *f_in       = prhs[2];
    const mxArray *lambda_in  = prhs[3];
    const mxArray *weight_in  = prhs[4];
    const mxArray *mass_in    = prhs[5];
    const mxArray *col_in     = prhs[6];
    const mxArray *opts_in    = (nrhs >= 8) ? prhs[7] : NULL;

    if (!mxIsStruct(hsp_in))     mexErrMsgIdAndTxt("qi_ls:type", "hspace must be a struct");
    if (!mxIsCell(mass_in))      mexErrMsgIdAndTxt("qi_ls:type", "mass_cell must be a cell array");
    if (!mxIsCell(col_in))       mexErrMsgIdAndTxt("qi_ls:type", "col_cell must be a cell array");

    /* --- hspace_info --- */
    int par_dim = field_scalar_int(hsp_in, "par_dim", 0);
    if (par_dim != 2 && par_dim != 3)
        mexErrMsgIdAndTxt("qi_ls:par_dim", "par_dim must be 2 or 3");
    int nlevels = field_scalar_int(hsp_in, "nlevels", 0);
    int ndof    = field_scalar_int(hsp_in, "ndof", 0);
    if (nlevels < 1 || ndof < 1)
        mexErrMsgIdAndTxt("qi_ls:dims", "nlevels and ndof must be positive");

    int n_deg = 0;
    int *deg = field_int_array(hsp_in, "degree", &n_deg, 0);
    if (n_deg < par_dim)
        mexErrMsgIdAndTxt("qi_ls:degree", "degree must have par_dim entries");

    int n_dpl = 0;
    int *ndof_per_level = field_int_array(hsp_in, "ndof_per_level", &n_dpl, 0);
    if (n_dpl < nlevels)
        mexErrMsgIdAndTxt("qi_ls:ndof_per_level", "ndof_per_level too short");

    int n_act_off = 0;
    int *active_offset = field_int_array(hsp_in, "active_offset", &n_act_off, 0);
    if (n_act_off < nlevels + 1)
        mexErrMsgIdAndTxt("qi_ls:active_offset", "active_offset must have nlevels+1 entries");

    int n_act = 0;
    int *active_indices = field_int_array(hsp_in, "active_indices", &n_act, 1); /* 1->0 based */
    if (n_act < ndof)
        mexErrMsgIdAndTxt("qi_ls:active_indices", "active_indices too short");

    int *ndof_dir = field_matrix_int_rowmajor(hsp_in, "ndof_dir", nlevels, par_dim, 0);
    int *nel_dir  = field_matrix_int_rowmajor(hsp_in, "nel_dir",  nlevels, par_dim, 0);

    mxArray *brk = mxGetField(hsp_in, 0, "breaks");
    if (!brk || !mxIsCell(brk))
        mexErrMsgIdAndTxt("qi_ls:breaks", "hspace.breaks (cell) is required");
    double *breaks_data = NULL; int *breaks_offset = NULL;
    flatten_breaks(brk, nlevels, par_dim, &breaks_data, &breaks_offset);

    /* Per-cell connectivity (first 1D basis index, 0-based) */
    mxArray *conn_in = mxGetField(hsp_in, 0, "conn_first");
    if (!conn_in || !mxIsCell(conn_in))
        mexErrMsgIdAndTxt("qi_ls:conn_first", "hspace.conn_first (cell) is required");
    int *conn_first_offset = (int*)mxMalloc((size_t)(nlevels * par_dim + 1) * sizeof(int));
    conn_first_offset[0] = 0;
    /* first pass: total length */
    for (int lv = 0; lv < nlevels; ++lv) {
        mxArray *clv = mxGetCell(conn_in, lv);
        if (!clv || !mxIsCell(clv) || (int)mxGetNumberOfElements(clv) < par_dim)
            mexErrMsgIdAndTxt("qi_ls:conn_first", "conn_first{%d} must be a 1x%d cell", lv+1, par_dim);
        for (int d = 0; d < par_dim; ++d) {
            mxArray *cd = mxGetCell(clv, d);
            if (!cd)
                mexErrMsgIdAndTxt("qi_ls:conn_first", "conn_first{%d}{%d} is empty", lv+1, d+1);
            int nel_lv_d = nel_dir[lv * par_dim + d];
            if ((int)mxGetNumberOfElements(cd) < nel_lv_d)
                mexErrMsgIdAndTxt("qi_ls:conn_first",
                    "conn_first{%d}{%d} length %d < nel_dir %d",
                    lv+1, d+1, (int)mxGetNumberOfElements(cd), nel_lv_d);
            conn_first_offset[lv * par_dim + d + 1] =
                conn_first_offset[lv * par_dim + d] + nel_lv_d;
        }
    }
    int conn_total = conn_first_offset[nlevels * par_dim];
    int *conn_first_data = (int*)mxMalloc((size_t)(conn_total > 0 ? conn_total : 1) * sizeof(int));
    /* second pass: copy and convert from 1-based to 0-based */
    for (int lv = 0; lv < nlevels; ++lv) {
        mxArray *clv = mxGetCell(conn_in, lv);
        for (int d = 0; d < par_dim; ++d) {
            mxArray *cd = mxGetCell(clv, d);
            int nel_lv_d = nel_dir[lv * par_dim + d];
            int off = conn_first_offset[lv * par_dim + d];
            if (mxIsDouble(cd)) {
                double *p = mxGetPr(cd);
                for (int c = 0; c < nel_lv_d; ++c) conn_first_data[off + c] = (int)p[c] - 1;
            } else if (mxIsInt32(cd)) {
                int *p = (int*)mxGetData(cd);
                for (int c = 0; c < nel_lv_d; ++c) conn_first_data[off + c] = p[c] - 1;
            } else {
                mexErrMsgIdAndTxt("qi_ls:conn_first", "conn_first{%d}{%d} must be double or int32", lv+1, d+1);
            }
        }
    }

    /* --- data --- */
    if (!mxIsDouble(coords_in) || mxIsSparse(coords_in))
        mexErrMsgIdAndTxt("qi_ls:coords", "data_coords must be dense double");
    int npoints = (int)mxGetM(coords_in);
    if ((int)mxGetN(coords_in) != par_dim)
        mexErrMsgIdAndTxt("qi_ls:coords", "data_coords must be npoints x par_dim");

    /* MATLAB column-major coords -> row-major */
    double *coords_cm = mxGetPr(coords_in);
    double *coords_rm = (double*)mxMalloc((size_t)npoints * par_dim * sizeof(double));
    for (int i = 0; i < npoints; ++i)
        for (int d = 0; d < par_dim; ++d)
            coords_rm[(size_t)i * par_dim + d] = coords_cm[i + (size_t)d * npoints];

    if (!mxIsDouble(f_in) || mxIsSparse(f_in))
        mexErrMsgIdAndTxt("qi_ls:f", "f must be dense double");
    int ncomp = (int)mxGetN(f_in);
    if ((int)mxGetM(f_in) != npoints)
        mexErrMsgIdAndTxt("qi_ls:f", "f must be npoints x ncomp");
    double *f_cm = mxGetPr(f_in);
    double *f_rm = (double*)mxMalloc((size_t)npoints * ncomp * sizeof(double));
    for (int i = 0; i < npoints; ++i)
        for (int c = 0; c < ncomp; ++c)
            f_rm[(size_t)i * ncomp + c] = f_cm[i + (size_t)c * npoints];

    double lambda = mxGetScalar(lambda_in);

    double *weight = NULL;
    if (!mxIsEmpty(weight_in)) {
        if (!mxIsDouble(weight_in) || (int)mxGetNumberOfElements(weight_in) != npoints)
            mexErrMsgIdAndTxt("qi_ls:weight", "weight must be empty or [npoints x 1] double");
        weight = mxGetPr(weight_in);
    }

    /* --- per-level mass and col matrices --- */
    if ((int)mxGetNumberOfElements(mass_in) < nlevels ||
        (int)mxGetNumberOfElements(col_in)  < nlevels)
        mexErrMsgIdAndTxt("qi_ls:cells", "mass_cell / col_cell must each have nlevels entries");

    qi_ls_level_precomp_t *level_pre = (qi_ls_level_precomp_t*)mxMalloc((size_t)nlevels * sizeof(qi_ls_level_precomp_t));
    int **alloc_mp = (int**)mxMalloc((size_t)nlevels * 4 * sizeof(int*));
    /* alloc_mp layout: per level lv, four pointers
     *   [4*lv + 0] mass.col_ptr
     *   [4*lv + 1] mass.row_ind
     *   [4*lv + 2] col.col_ptr
     *   [4*lv + 3] col.row_ind */
    for (int lv = 0; lv < nlevels; ++lv) {
        memset(&level_pre[lv], 0, sizeof(qi_ls_level_precomp_t));
        if (ndof_per_level[lv] <= 0) {
            alloc_mp[4*lv + 0] = NULL;
            alloc_mp[4*lv + 1] = NULL;
            alloc_mp[4*lv + 2] = NULL;
            alloc_mp[4*lv + 3] = NULL;
            continue;
        }
        mxArray *mlv = mxGetCell(mass_in, lv);
        mxArray *clv = mxGetCell(col_in,  lv);
        copy_sparse_to_csc(mlv, &level_pre[lv].mass,
                           &alloc_mp[4*lv + 0], &alloc_mp[4*lv + 1]);
        copy_sparse_to_csc(clv, &level_pre[lv].col,
                           &alloc_mp[4*lv + 2], &alloc_mp[4*lv + 3]);
    }

    /* --- options --- */
    qi_ls_options_t opts; qi_ls_options_defaults(&opts);
    opts.lambda = lambda;
    if (opts_in && mxIsStruct(opts_in)) {
        opts.n_threads = field_scalar_int(opts_in, "n_threads", opts.n_threads);
        opts.verbose   = field_scalar_int(opts_in, "verbose",   opts.verbose);
    }

    /* --- pack and call --- */
    qi_ls_hspace_t H = (qi_ls_hspace_t){
        .par_dim            = par_dim,
        .nlevels            = nlevels,
        .degree             = deg,
        .nel_dir            = nel_dir,
        .ndof_dir           = ndof_dir,
        .breaks_data        = breaks_data,
        .breaks_offset      = breaks_offset,
        .conn_first_data    = conn_first_data,
        .conn_first_offset  = conn_first_offset,
        .ndof               = ndof,
        .ndof_per_level     = ndof_per_level,
        .active_offset      = active_offset,
        .active_indices     = active_indices
    };
    qi_ls_data_t D = (qi_ls_data_t){
        .npoints = npoints, .ncomp = ncomp, .par_dim = par_dim,
        .coords  = coords_rm, .f = f_rm, .weight = weight
    };

    /* Output: ndof x ncomp, MATLAB stores column-major */
    plhs[0] = mxCreateDoubleMatrix(ndof, ncomp, mxREAL);
    double *coeffs_cm = mxGetPr(plhs[0]);

    /* libqi writes row-major; we copy to column-major */
    double *coeffs_rm = (double*)mxMalloc((size_t)ndof * ncomp * sizeof(double));
    qi_status_t st = qi_localLS_compute(&H, &D, level_pre, &opts, coeffs_rm);
    if (st != QI_OK) {
        mexErrMsgIdAndTxt("qi_ls:compute", "qi_localLS_compute returned %d", (int)st);
    }
    for (int k = 0; k < ndof; ++k)
        for (int c = 0; c < ncomp; ++c)
            coeffs_cm[(size_t)c * ndof + k] = coeffs_rm[(size_t)k * ncomp + c];
    mxFree(coeffs_rm);

    /* cleanup */
    for (int lv = 0; lv < nlevels; ++lv) {
        if (alloc_mp[4*lv + 0]) mxFree(alloc_mp[4*lv + 0]);
        if (alloc_mp[4*lv + 1]) mxFree(alloc_mp[4*lv + 1]);
        if (alloc_mp[4*lv + 2]) mxFree(alloc_mp[4*lv + 2]);
        if (alloc_mp[4*lv + 3]) mxFree(alloc_mp[4*lv + 3]);
    }
    mxFree(alloc_mp);
    mxFree(level_pre);
    mxFree(coords_rm); mxFree(f_rm);
    mxFree(deg); mxFree(ndof_per_level); mxFree(active_offset);
    mxFree(active_indices); mxFree(ndof_dir); mxFree(nel_dir);
    mxFree(breaks_data); mxFree(breaks_offset);
    mxFree(conn_first_data); mxFree(conn_first_offset);
}

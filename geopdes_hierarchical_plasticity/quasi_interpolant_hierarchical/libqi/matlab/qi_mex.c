/*
 * qi_mex.c — MATLAB MEX wrapper around qi_compute_coeffs.
 *
 * Usage from MATLAB:
 *
 *   coeffs = qi_mex(hspace, hmsh, data);
 *   coeffs = qi_mex(hspace, hmsh, data, opts);
 *
 *   [coeffs, diag] = qi_mex(...);
 *
 * Where:
 *
 *   hspace.degree         = [d_x d_y]
 *   hspace.nlevels        = M
 *   hspace.ndof           = total number of THB-splines
 *   hspace.ndof_per_level = [M x 1] integer
 *   hspace.active{ℓ}      = column-vector of 1-based active indices at level ℓ
 *
 *   hmsh.nel_dir{ℓ}       = [nel_x nel_y] at level ℓ
 *   hmsh.breaks{ℓ}{1}     = breakpoints in x at level ℓ  (optional)
 *   hmsh.breaks{ℓ}{2}     = breakpoints in y at level ℓ  (optional)
 *
 *   data.x  = [n x 1] in [0,1]
 *   data.y  = [n x 1] in [0,1]
 *   data.f  = [n x m] (m = ncomp)
 *
 *   opts (optional struct):
 *     opts.max_degree
 *     opts.sigma_threshold
 *     opts.n_threads
 *     opts.verbose
 *
 * Build:
 *
 *   make mex
 *   % or, from MATLAB:
 *   mex -R2018a -Iinclude -Isrc CFLAGS='$CFLAGS -fopenmp -O3' \
 *       LDFLAGS='$LDFLAGS -fopenmp' \
 *       matlab/qi_mex.c src/qi.c src/qi_linalg.c src/qi_bspline.c
 */

#include "mex.h"
#include "matrix.h"
#include "qi.h"
#include <string.h>
#include <stdlib.h>

/* helper: read a struct field that must be a real scalar, return as int */
static int field_as_int(const mxArray *s, const char *name, int dflt)
{
    mxArray *f = mxGetField(s, 0, name);
    if (!f || mxIsEmpty(f)) return dflt;
    return (int)mxGetScalar(f);
}

/* helper: read an integer 1-based vector field, copy & convert to 0-based */
static int *field_as_int_array_1based(const mxArray *s, const char *name,
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
        mexErrMsgIdAndTxt("qi:badtype", "field %s must be double or int32", name);
        return NULL;
    }
    *n_out = (int)n;
    return out;
}

/* extract per-level data from cell-array fields */
static int *gather_active(const mxArray *hspace, int nlevels,
                          int **offsets_out, int *ndof_total)
{
    mxArray *active = mxGetField(hspace, 0, "active");
    if (!active || !mxIsCell(active))
        mexErrMsgIdAndTxt("qi:badactive",
                          "hspace.active must be a cell array");

    int *offsets = (int*)mxMalloc((nlevels+1) * sizeof(int));
    offsets[0] = 0;
    int total = 0;
    for (int lv = 0; lv < nlevels; ++lv) {
        mxArray *c = mxGetCell(active, lv);
        size_t n = c ? mxGetNumberOfElements(c) : 0;
        total += (int)n;
        offsets[lv+1] = total;
    }
    int *idx = (int*)mxMalloc((total > 0 ? total : 1) * sizeof(int));
    int j = 0;
    for (int lv = 0; lv < nlevels; ++lv) {
        mxArray *c = mxGetCell(active, lv);
        if (!c) continue;
        size_t n = mxGetNumberOfElements(c);
        if (mxIsDouble(c)) {
            double *p = mxGetPr(c);
            for (size_t i = 0; i < n; ++i) idx[j++] = (int)p[i] - 1;
        } else if (mxIsInt32(c)) {
            int *p = (int*)mxGetData(c);
            for (size_t i = 0; i < n; ++i) idx[j++] = p[i] - 1;
        } else {
            mexErrMsgIdAndTxt("qi:badtype", "active{l} must be double or int32");
        }
    }
    *offsets_out = offsets;
    *ndof_total  = total;
    return idx;
}

/* gather per-level mesh sizes from hmsh.nel_dir{lv} = [nelx nely] */
static void gather_mesh(const mxArray *hmsh, int nlevels,
                        int **nel_x_out, int **nel_y_out)
{
    mxArray *nel = mxGetField(hmsh, 0, "nel_dir");
    if (!nel)
        mexErrMsgIdAndTxt("qi:badmesh", "hmsh.nel_dir is required");

    int *nx = (int*)mxMalloc(nlevels * sizeof(int));
    int *ny = (int*)mxMalloc(nlevels * sizeof(int));

    if (mxIsCell(nel)) {
        if ((int)mxGetNumberOfElements(nel) < nlevels)
            mexErrMsgIdAndTxt("qi:badmesh",
                              "hmsh.nel_dir cell length < nlevels");
        for (int lv = 0; lv < nlevels; ++lv) {
            mxArray *c = mxGetCell(nel, lv);
            if (!c || mxGetNumberOfElements(c) < 2)
                mexErrMsgIdAndTxt("qi:badmesh",
                                  "hmsh.nel_dir{%d} needs 2 entries", lv+1);
            double *p = mxGetPr(c);
            nx[lv] = (int)p[0]; ny[lv] = (int)p[1];
        }
    } else {
        /* matrix: nlevels x 2 */
        if ((int)mxGetM(nel) < nlevels || mxGetN(nel) < 2)
            mexErrMsgIdAndTxt("qi:badmesh",
                              "hmsh.nel_dir matrix must be nlevels x 2");
        double *p = mxGetPr(nel);
        size_t M = mxGetM(nel);
        for (int lv = 0; lv < nlevels; ++lv) {
            nx[lv] = (int)p[lv];
            ny[lv] = (int)p[lv + M];
        }
    }
    *nel_x_out = nx;
    *nel_y_out = ny;
}

/* gather optional non-uniform breaks from hmsh.breaks{lv}{1|2} */
static void gather_breaks(const mxArray *hmsh, int nlevels,
                          double **dx_out, int **off_x_out,
                          double **dy_out, int **off_y_out)
{
    *dx_out = NULL; *off_x_out = NULL;
    *dy_out = NULL; *off_y_out = NULL;
    mxArray *brk = mxGetField(hmsh, 0, "breaks");
    if (!brk || !mxIsCell(brk)) return;
    if ((int)mxGetNumberOfElements(brk) < nlevels) return;

    int *off_x = (int*)mxMalloc((nlevels+1)*sizeof(int));
    int *off_y = (int*)mxMalloc((nlevels+1)*sizeof(int));
    off_x[0] = 0; off_y[0] = 0;

    /* first pass: total sizes */
    for (int lv = 0; lv < nlevels; ++lv) {
        mxArray *c = mxGetCell(brk, lv);
        if (!c || !mxIsCell(c) || mxGetNumberOfElements(c) < 2) {
            mxFree(off_x); mxFree(off_y);
            return;   /* fall back to uniform */
        }
        mxArray *bx = mxGetCell(c, 0);
        mxArray *by = mxGetCell(c, 1);
        off_x[lv+1] = off_x[lv] + (int)mxGetNumberOfElements(bx);
        off_y[lv+1] = off_y[lv] + (int)mxGetNumberOfElements(by);
    }

    double *dx = (double*)mxMalloc(off_x[nlevels]*sizeof(double));
    double *dy = (double*)mxMalloc(off_y[nlevels]*sizeof(double));
    for (int lv = 0; lv < nlevels; ++lv) {
        mxArray *c  = mxGetCell(brk, lv);
        mxArray *bx = mxGetCell(c, 0);
        mxArray *by = mxGetCell(c, 1);
        memcpy(dx + off_x[lv], mxGetPr(bx),
               (off_x[lv+1]-off_x[lv])*sizeof(double));
        memcpy(dy + off_y[lv], mxGetPr(by),
               (off_y[lv+1]-off_y[lv])*sizeof(double));
    }
    *dx_out = dx; *off_x_out = off_x;
    *dy_out = dy; *off_y_out = off_y;
}

/* ------------------------------------------------------------------------ */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs < 3 || nrhs > 4)
        mexErrMsgIdAndTxt("qi:nrhs",
            "Usage: coeffs = qi_mex(hspace, hmsh, data [, opts])");

    const mxArray *hspace = prhs[0];
    const mxArray *hmsh   = prhs[1];
    const mxArray *data   = prhs[2];
    const mxArray *optms  = (nrhs >= 4) ? prhs[3] : NULL;

    if (!mxIsStruct(hspace) || !mxIsStruct(hmsh) || !mxIsStruct(data))
        mexErrMsgIdAndTxt("qi:type", "hspace, hmsh, data must be structs");

    /* --- hspace --- */
    int nlevels = field_as_int(hspace, "nlevels", -1);
    if (nlevels < 1) mexErrMsgIdAndTxt("qi:nlevels", "hspace.nlevels invalid");

    int n_deg = 0;
    int *deg = field_as_int_array_1based(hspace, "degree", &n_deg, 0);
    if (n_deg < 2) mexErrMsgIdAndTxt("qi:degree", "hspace.degree must have 2 entries");
    int dx = deg[0], dy = deg[1];

    int n_dofperlev = 0;
    int *dofperlev = field_as_int_array_1based(hspace, "ndof_per_level", &n_dofperlev, 0);
    if (n_dofperlev < nlevels)
        mexErrMsgIdAndTxt("qi:ndof", "hspace.ndof_per_level too short");

    int *active_offset = NULL, ndof = 0;
    int *active = gather_active(hspace, nlevels, &active_offset, &ndof);

    /* --- hmsh --- */
    int *nel_x = NULL, *nel_y = NULL;
    gather_mesh(hmsh, nlevels, &nel_x, &nel_y);

    double *dxbreak = NULL, *dybreak = NULL;
    int    *offx    = NULL, *offy    = NULL;
    gather_breaks(hmsh, nlevels, &dxbreak, &offx, &dybreak, &offy);

    /* --- data --- */
    mxArray *xfld = mxGetField(data, 0, "x");
    mxArray *yfld = mxGetField(data, 0, "y");
    mxArray *ffld = mxGetField(data, 0, "f");
    if (!xfld || !yfld || !ffld)
        mexErrMsgIdAndTxt("qi:data", "data must have fields x, y, f");
    int npts  = (int)mxGetNumberOfElements(xfld);
    int ncomp = (int)mxGetN(ffld);
    if ((int)mxGetM(ffld) != npts) {
        if ((int)mxGetN(ffld) == npts && (int)mxGetM(ffld) > 0) {
            /* row-major / column orientation: transpose */
            ncomp = (int)mxGetM(ffld);
        } else
            mexErrMsgIdAndTxt("qi:data", "data.f must be n x ncomp");
    }
    double *xs = mxGetPr(xfld), *ys = mxGetPr(yfld);
    double *fs_in = mxGetPr(ffld);
    /* convert column-major MATLAB f (n x ncomp) -> row-major (ncomp inner) */
    double *fs = (double*)mxMalloc((size_t)npts * ncomp * sizeof(double));
    for (int i = 0; i < npts; ++i)
        for (int c = 0; c < ncomp; ++c)
            fs[i*ncomp + c] = fs_in[c*npts + i];

    /* --- options --- */
    qi_options_t opts; qi_options_defaults(&opts);
    if (optms && mxIsStruct(optms)) {
        opts.max_degree      = field_as_int(optms, "max_degree", opts.max_degree);
        opts.n_threads       = field_as_int(optms, "n_threads",  opts.n_threads);
        opts.verbose         = field_as_int(optms, "verbose",    opts.verbose);
        mxArray *st = mxGetField(optms, 0, "sigma_threshold");
        if (st && !mxIsEmpty(st)) opts.sigma_threshold = mxGetScalar(st);
    }

    /* --- pack hspace --- */
    qi_hspace_t H;
    H.degree_x         = dx;
    H.degree_y         = dy;
    H.nlevels          = nlevels;
    H.nel_x            = nel_x;
    H.nel_y            = nel_y;
    H.breaks_x_data    = dxbreak;
    H.breaks_x_offset  = offx;
    H.breaks_y_data    = dybreak;
    H.breaks_y_offset  = offy;
    H.ndof             = ndof;
    H.ndof_per_level   = dofperlev;
    H.active_offset    = active_offset;
    H.active_indices   = active;

    qi_data_t D = { npts, ncomp, xs, ys, fs };

    /* --- output: coeffs (ndof x ncomp) column-major MATLAB matrix --- */
    plhs[0] = mxCreateDoubleMatrix(ndof, ncomp, mxREAL);
    double *coeffs_col = mxGetPr(plhs[0]);

    /* compute into row-major buffer, then copy out as column-major */
    double *coeffs_row = (double*)mxMalloc((size_t)ndof * ncomp * sizeof(double));

    qi_diag_t diag = {0};
    int *deg_used  = NULL;
    int *nlocal    = NULL;
    double *condest = NULL;
    int *level_arr  = NULL;
    if (nlhs >= 2) {
        deg_used  = (int*)mxMalloc(ndof * sizeof(int));
        nlocal    = (int*)mxMalloc(ndof * sizeof(int));
        condest   = (double*)mxMalloc(ndof * sizeof(double));
        level_arr = (int*)mxMalloc(ndof * sizeof(int));
        diag.deg_used = deg_used;
        diag.nlocal   = nlocal;
        diag.cond_estimate = condest;
        diag.level    = level_arr;
    }

    qi_status_t st = qi_compute_coeffs(&H, &D, &opts, coeffs_row, &diag);
    if (st != QI_OK)
        mexErrMsgIdAndTxt("qi:compute",
                          "qi_compute_coeffs returned %d", (int)st);

    for (int k = 0; k < ndof; ++k)
        for (int c = 0; c < ncomp; ++c)
            coeffs_col[c*ndof + k] = coeffs_row[k*ncomp + c];

    if (nlhs >= 2) {
        const char *fields[] = {"deg_used","nlocal","cond_estimate","level"};
        plhs[1] = mxCreateStructMatrix(1, 1, 4, fields);
        mxArray *a;
        a = mxCreateDoubleMatrix(ndof, 1, mxREAL);
        for (int k = 0; k < ndof; ++k) mxGetPr(a)[k] = (double)deg_used[k];
        mxSetField(plhs[1], 0, "deg_used", a);
        a = mxCreateDoubleMatrix(ndof, 1, mxREAL);
        for (int k = 0; k < ndof; ++k) mxGetPr(a)[k] = (double)nlocal[k];
        mxSetField(plhs[1], 0, "nlocal", a);
        a = mxCreateDoubleMatrix(ndof, 1, mxREAL);
        memcpy(mxGetPr(a), condest, ndof*sizeof(double));
        mxSetField(plhs[1], 0, "cond_estimate", a);
        a = mxCreateDoubleMatrix(ndof, 1, mxREAL);
        for (int k = 0; k < ndof; ++k) mxGetPr(a)[k] = (double)(level_arr[k] + 1);
        mxSetField(plhs[1], 0, "level", a);

        mxFree(deg_used); mxFree(nlocal); mxFree(condest); mxFree(level_arr);
    }

    /* cleanup */
    mxFree(coeffs_row);
    mxFree(fs);
    mxFree(deg);
    mxFree(dofperlev);
    mxFree(active_offset);
    mxFree(active);
    mxFree(nel_x); mxFree(nel_y);
    if (dxbreak) mxFree(dxbreak);
    if (dybreak) mxFree(dybreak);
    if (offx) mxFree(offx);
    if (offy) mxFree(offy);
}

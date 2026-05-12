/*
 * qi.c — Main implementation of the THB-spline quasi-interpolant.
 *
 * Contains the parallel loop over active basis functions and the per-basis
 * computation (support box, local data selection, polynomial LS,
 * polynomial-to-B-spline basis change, coefficient extraction).
 *
 * v1.1 optimisations:
 *   - Spatial grid index for O(N_local) data selection (was O(N_total)).
 *   - Pre-computed clamped knot vectors per level (avoid per-DOF rebuild).
 *   - Power-table pre-computation in basis change (avoid O(p) inner loops).
 *   - Factor Nx*Ny from component loop in qi_eval_tp.
 *   - restrict qualifiers throughout.
 */
#include "../include/qi.h"
#include "qi_linalg.h"
#include "qi_bspline.h"
#include "qi_spatial.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef _OPENMP
#include <omp.h>
#endif

/* ------------------------------------------------------------------------ */
const char *qi_version(void)
{
    return "libqi 1.1.0  (Hierarchical Quasi-Interpolant; OpenMP; spatial-index)";
}

void qi_options_defaults(qi_options_t *opts)
{
    opts->max_degree      = 0;
    opts->sigma_threshold = 0.0;
    opts->n_threads       = 0;
    opts->verbose         = 0;
}

/* ------------------------------------------------------------------------ */
/* helpers                                                                  */
/* ------------------------------------------------------------------------ */

static inline int qi_locate_level(const int *ndof_offset, int nlevels,
                                  int k, int *kloc_out)
{
    int lo = 0, hi = nlevels - 1;
    while (lo < hi) {
        int mid = (lo + hi + 1) / 2;
        if (ndof_offset[mid] <= k) lo = mid;
        else                       hi = mid - 1;
    }
    *kloc_out = k - ndof_offset[lo];
    return lo;
}

static inline double qi_break_x(const qi_hspace_t *h, int lev, int i)
{
    if (h->breaks_x_data)
        return h->breaks_x_data[h->breaks_x_offset[lev] + i];
    return (double)i / (double)h->nel_x[lev];
}
static inline double qi_break_y(const qi_hspace_t *h, int lev, int i)
{
    if (h->breaks_y_data)
        return h->breaks_y_data[h->breaks_y_offset[lev] + i];
    return (double)i / (double)h->nel_y[lev];
}

static void qi_make_clamped_knots(const qi_hspace_t *h, int lev, int axis,
                                  double *U_out, int *nknots_out)
{
    int p   = (axis == 0) ? h->degree_x : h->degree_y;
    int nel = (axis == 0) ? h->nel_x[lev] : h->nel_y[lev];
    int nb  = nel + 1;
    int nknots = nb + 2*p;
    int idx = 0;
    double left = (axis == 0) ? qi_break_x(h, lev, 0) : qi_break_y(h, lev, 0);
    for (int j = 0; j <= p; ++j) U_out[idx++] = left;
    for (int j = 1; j < nb - 1; ++j)
        U_out[idx++] = (axis == 0) ? qi_break_x(h, lev, j) : qi_break_y(h, lev, j);
    double right = (axis == 0) ? qi_break_x(h, lev, nb-1) : qi_break_y(h, lev, nb-1);
    for (int j = 0; j <= p; ++j) U_out[idx++] = right;
    *nknots_out = nknots;
}

/* Pre-computed knot vectors for all levels (built once before parallel loop) */
typedef struct {
    double *Ux;   int nUx;
    double *Uy;   int nUy;
} qi_level_knots_t;

static int enumerate_monomials(int pd, int *dx_arr, int *dy_arr)
{
    int n = 0;
    for (int i = 0; i <= pd; ++i)
        for (int j = 0; j + i <= pd; ++j) {
            dx_arr[n] = i;
            dy_arr[n] = j;
            ++n;
        }
    return n;
}

static long long binomial_int(int n, int k)
{
    if (k < 0 || k > n) return 0;
    if (k == 0 || k == n) return 1;
    if (k > n - k) k = n - k;
    long long r = 1;
    for (int i = 1; i <= k; ++i)
        r = r * (n - i + 1) / i;
    return r;
}

/* ------------------------------------------------------------------------ */
/* Per-basis-function computation                                           */
/* ------------------------------------------------------------------------ */

#define QI_MAX_DEGREE 4
#define QI_MAX_NPB    ((QI_MAX_DEGREE+1)*(QI_MAX_DEGREE+2)/2)
#define QI_MAX_NL1D   (2*QI_MAX_DEGREE + 1)
#define QI_MAX_NL2D   (QI_MAX_NL1D*QI_MAX_NL1D)

typedef struct {
    int    *loc_idx;  int loc_idx_cap;
    double *xloc;     int xloc_cap;
    double *yloc;     int yloc_cap;
    double *floc;     int floc_cap;
    double *A;        int A_cap;
    double *coef_p;   int coef_p_cap;
    double *coef_p_xy;int coef_p_xy_cap;
    double *aux_x;    int aux_x_cap;
    double *aux_y;    int aux_y_cap;
    double *bbasis;   int bbasis_cap;
    double *poly_emat;int poly_emat_cap;
    double *bbasis_T; int bbasis_T_cap;
    double *rhs;      int rhs_cap;
    int    *piv;      int piv_cap;
    double *work;     int work_cap;
} qi_scratch_t;

static void qi_scratch_init(qi_scratch_t *s) { memset(s, 0, sizeof(*s)); }

static void *qi_xrealloc(void *p, int cap_old, int cap_new, size_t elem)
{
    if (cap_new <= cap_old) return p;
    return realloc(p, (size_t)cap_new * elem);
}
#define QI_ENSURE(buf, cap, need, elem)                              \
    do {                                                             \
        if ((cap) < (need)) {                                        \
            int new_cap = (need) * 2;                                \
            void *ptr = qi_xrealloc((buf), (cap), new_cap, (elem));  \
            if (!ptr) return QI_ERR_OUT_OF_MEMORY;                   \
            (buf) = ptr; (cap) = new_cap;                            \
        }                                                            \
    } while (0)

static void qi_scratch_free(qi_scratch_t *s)
{
    free(s->loc_idx);
    free(s->xloc); free(s->yloc); free(s->floc);
    free(s->A);
    free(s->coef_p); free(s->coef_p_xy);
    free(s->aux_x); free(s->aux_y);
    free(s->bbasis); free(s->poly_emat);
    free(s->bbasis_T); free(s->rhs);
    free(s->piv); free(s->work);
}

/* ------------------------------------------------------------------------ */
/* Solve the QI for one basis function k.                                   */
/* ------------------------------------------------------------------------ */
static int qi_solve_one_basis(const qi_hspace_t       *hs,
                              const qi_data_t         *dat,
                              const qi_options_t      *opts,
                              const int               *level_offset,
                              const qi_spatial_t      *spatial,
                              const qi_level_knots_t  *lev_knots,
                              int                      k,
                              int                      max_degree,
                              double                   sigma_threshold,
                              qi_scratch_t            *S,
                              double                  *out_coeff,
                              int                     *deg_used,
                              int                     *nlocal_out,
                              double                  *cond_out,
                              int                     *level_out)
{
    const int dx = hs->degree_x;
    const int dy = hs->degree_y;
    const int ncomp = dat->ncomp;

    int kloc;
    int lev = qi_locate_level(level_offset, hs->nlevels, k, &kloc);
    if (level_out) *level_out = lev;

    int act_idx = hs->active_indices[level_offset[lev] + kloc];
    int ndof_dir_x = hs->nel_x[lev] + dx;
    int ndof_dir_y = hs->nel_y[lev] + dy;
    int k_x = act_idx % ndof_dir_x;
    int k_y = act_idx / ndof_dir_x;

    int nx = hs->nel_x[lev];
    int ny = hs->nel_y[lev];

    /* support box (0-based break indices) */
    int mux = (k_x - dx > 0) ? (k_x - dx) : 0;
    int nux = (k_x + 1 < nx) ? (k_x + 1) : nx;
    int muy = (k_y - dy > 0) ? (k_y - dy) : 0;
    int nuy = (k_y + 1 < ny) ? (k_y + 1) : ny;

    double x_mu = qi_break_x(hs, lev, mux);
    double x_nu = qi_break_x(hs, lev, nux);
    double y_mu = qi_break_y(hs, lev, muy);
    double y_nu = qi_break_y(hs, lev, nuy);
    double Dx   = x_nu - x_mu;
    double Dy   = y_nu - y_mu;

    if (Dx <= 0.0 || Dy <= 0.0) return QI_ERR_INVALID_ARG;

    double cx = 0.5 * (x_mu + x_nu);
    double cy = 0.5 * (y_mu + y_nu);
    double r0 = 0.5 * sqrt(Dx*Dx + Dy*Dy);

    int ref_lev = (lev >= 1) ? (lev - 1) : lev;
    double ref_nx = (double)hs->nel_x[ref_lev];
    double ref_ny = (double)hs->nel_y[ref_lev];
    double ratio_d = sqrt(((double)(dx+1)/ref_nx)*((double)(dx+1)/ref_nx)
                        + ((double)(dy+1)/ref_ny)*((double)(dy+1)/ref_ny))
                     / (2.0 * r0);
    int ratio_max = (int)ceil(ratio_d);
    if (ratio_max < 1) ratio_max = 1;
    int max_iter = 2*ratio_max + 1;

    int pd = max_degree;
    int npb_target = (pd+1)*(pd+2)/2;

    /* ---- local data selection via spatial index ---- */
    QI_ENSURE(S->loc_idx, S->loc_idx_cap, dat->npoints, sizeof(int));

    int nloc = 0;
    int iter = 1;
    double r = r0;
    while (1) {
        nloc = qi_spatial_query_disk_2d(spatial, cx, cy, r,
                                        dat->x, dat->y,
                                        S->loc_idx, S->loc_idx_cap);
        while (nloc < (pd+1)*(pd+2)/2 && pd > 0) --pd;
        npb_target = (pd+1)*(pd+2)/2;
        if (nloc >= npb_target) break;
        if (iter > max_iter) break;
        ++iter;
        r = (double)iter * r0;
    }

    if (nloc < npb_target) {
        nloc = dat->npoints;
        for (int i = 0; i < nloc; ++i) S->loc_idx[i] = i;
        while (nloc < (pd+1)*(pd+2)/2 && pd > 0) --pd;
        npb_target = (pd+1)*(pd+2)/2;
        if (nloc < npb_target) return QI_ERR_NOT_ENOUGH_DATA;
    }
    if (nlocal_out) *nlocal_out = nloc;

    /* extract local data, scaled to [0,1]^2 */
    QI_ENSURE(S->xloc, S->xloc_cap, nloc,        sizeof(double));
    QI_ENSURE(S->yloc, S->yloc_cap, nloc,        sizeof(double));
    QI_ENSURE(S->floc, S->floc_cap, nloc*ncomp,  sizeof(double));

    for (int i = 0; i < nloc; ++i) {
        int idx = S->loc_idx[i];
        S->xloc[i] = (dat->x[idx] - x_mu) / Dx;
        S->yloc[i] = (dat->y[idx] - y_mu) / Dy;
        for (int c = 0; c < ncomp; ++c)
            S->floc[i*ncomp + c] = dat->f[idx*ncomp + c];
    }

    /* ---- use pre-computed knot vectors ---- */
    const double *Ux = lev_knots[lev].Ux;
    const double *Uy = lev_knots[lev].Uy;
    int nUx = lev_knots[lev].nUx;
    int nUy = lev_knots[lev].nUy;

    /* local 2D B-spline tile */
    int nl_x = nux - mux + dx;
    int nl_y = nuy - muy + dy;
    int ncpx = nl_x;
    int ncpy = nl_y;
    int npts_grid = ncpx * ncpy;
    int nl2d  = nl_x * nl_y;

    /* collocation points */
    double col_x[64], col_y[64];
    if (ncpx >= 64 || ncpy >= 64) return QI_ERR_INVALID_ARG;
    if (ncpx == 1) col_x[0] = 0.5*(x_mu + x_nu);
    else for (int i = 0; i < ncpx; ++i)
        col_x[i] = x_mu + (x_nu - x_mu) * (double)i / (double)(ncpx - 1);
    if (ncpy == 1) col_y[0] = 0.5*(y_mu + y_nu);
    else for (int i = 0; i < ncpy; ++i)
        col_y[i] = y_mu + (y_nu - y_mu) * (double)i / (double)(ncpy - 1);

    /* spcol */
    int aux_x_need = ndof_dir_x * ncpx;
    int aux_y_need = ndof_dir_y * ncpy;
    QI_ENSURE(S->aux_x, S->aux_x_cap, aux_x_need, sizeof(double));
    QI_ENSURE(S->aux_y, S->aux_y_cap, aux_y_need, sizeof(double));
    qi_bs_spcol(Ux, nUx, dx, col_x, ncpx, S->aux_x);
    qi_bs_spcol(Uy, nUy, dy, col_y, ncpy, S->aux_y);

    /* bbasis (nl2d x npts_grid) */
    QI_ENSURE(S->bbasis, S->bbasis_cap, nl2d * npts_grid, sizeof(double));
    QI_ENSURE(S->poly_emat, S->poly_emat_cap, QI_MAX_NPB * npts_grid, sizeof(double));
    QI_ENSURE(S->bbasis_T, S->bbasis_T_cap, npts_grid * nl2d, sizeof(double));
    QI_ENSURE(S->rhs, S->rhs_cap, npts_grid * ncomp, sizeof(double));
    QI_ENSURE(S->piv, S->piv_cap, npts_grid, sizeof(int));

    for (int h_loc = 0; h_loc < nl_x; ++h_loc) {
        int h = mux + h_loc;
        for (int q_loc = 0; q_loc < nl_y; ++q_loc) {
            int q = muy + q_loc;
            int row = h_loc * nl_y + q_loc;
            for (int jj = 0; jj < ncpy; ++jj)
                for (int ii = 0; ii < ncpx; ++ii) {
                    int col = jj * ncpx + ii;
                    S->bbasis[row * npts_grid + col]
                        = S->aux_x[h * ncpx + ii] * S->aux_y[q * ncpy + jj];
                }
        }
    }

    /* degree descent loop */
    int dx_arr[QI_MAX_NPB], dy_arr[QI_MAX_NPB];
    int    chosen_pd       = -1;
    double chosen_sigmamin = 0.0;

    for (int try_pd = pd; try_pd >= 0; --try_pd) {
        int npb = enumerate_monomials(try_pd, dx_arr, dy_arr);
        if (nloc < npb) continue;

        /* Build A (nloc x npb) */
        QI_ENSURE(S->A, S->A_cap, nloc * npb, sizeof(double));
        for (int i = 0; i < nloc; ++i) {
            double u = S->xloc[i];
            double v = S->yloc[i];
            double upow[QI_MAX_DEGREE+1];
            double vpow[QI_MAX_DEGREE+1];
            upow[0] = 1.0; vpow[0] = 1.0;
            for (int p = 1; p <= try_pd; ++p) {
                upow[p] = upow[p-1] * u;
                vpow[p] = vpow[p-1] * v;
            }
            for (int p = 0; p < npb; ++p)
                S->A[i*npb + p] = upow[dx_arr[p]] * vpow[dy_arr[p]];
        }

        QI_ENSURE(S->coef_p, S->coef_p_cap, npb * ncomp, sizeof(double));
        QI_ENSURE(S->coef_p_xy, S->coef_p_xy_cap, npb * ncomp, sizeof(double));

        int wneed = 2*npb*npb + npb + 2*npb*ncomp;
        QI_ENSURE(S->work, S->work_cap, wneed, sizeof(double));

        double sigma_min = 0.0;
        int rc = qi_la_pinv_solve(nloc, npb, ncomp, S->A, S->floc,
                                  S->coef_p, &sigma_min,
                                  S->work, S->work_cap);
        if (rc != 0) continue;

        double MSV2 = (sigma_min > 0.0) ? (1.0 / sigma_min) : 1e300;
        if (MSV2 >= sigma_threshold) continue;

        /* ---- Basis change with power tables (optimisation #4) ---- */
        /* Pre-compute Dx^a, Dy^b, (-x_mu)^k, (-y_mu)^k tables */
        double Dx_pow[QI_MAX_DEGREE+1], Dy_pow[QI_MAX_DEGREE+1];
        double neg_xmu_pow[QI_MAX_DEGREE+1], neg_ymu_pow[QI_MAX_DEGREE+1];
        Dx_pow[0] = 1.0;  Dy_pow[0] = 1.0;
        neg_xmu_pow[0] = 1.0;  neg_ymu_pow[0] = 1.0;
        for (int p = 1; p <= try_pd; ++p) {
            Dx_pow[p]      = Dx_pow[p-1] * Dx;
            Dy_pow[p]      = Dy_pow[p-1] * Dy;
            neg_xmu_pow[p] = neg_xmu_pow[p-1] * (-x_mu);
            neg_ymu_pow[p] = neg_ymu_pow[p-1] * (-y_mu);
        }

        /* Pre-compute binomial table C(a,r) for a,r <= try_pd */
        double binom_tbl[QI_MAX_DEGREE+1][QI_MAX_DEGREE+1];
        for (int a = 0; a <= try_pd; ++a)
            for (int r2 = 0; r2 <= a; ++r2)
                binom_tbl[a][r2] = (double)binomial_int(a, r2);

        for (int idx_xy = 0; idx_xy < npb; ++idx_xy) {
            int r2 = dx_arr[idx_xy];
            int s2 = dy_arr[idx_xy];
            for (int c = 0; c < ncomp; ++c) {
                double sum = 0.0;
                for (int idx_uv = 0; idx_uv < npb; ++idx_uv) {
                    int a = dx_arr[idx_uv];
                    int b = dy_arr[idx_uv];
                    if (a < r2 || b < s2) continue;
                    double t = S->coef_p[idx_uv*ncomp + c];
                    t /= (Dx_pow[a] * Dy_pow[b]);
                    sum += t * binom_tbl[a][r2] * binom_tbl[b][s2]
                             * neg_xmu_pow[a - r2] * neg_ymu_pow[b - s2];
                }
                S->coef_p_xy[idx_xy*ncomp + c] = sum;
            }
        }

        /* ---- poly_emat with pre-computed powers ---- */
        /* Pre-compute col_x[ii]^a and col_y[jj]^b */
        double col_x_pow[64][QI_MAX_DEGREE+1];
        double col_y_pow[64][QI_MAX_DEGREE+1];
        for (int ii = 0; ii < ncpx; ++ii) {
            col_x_pow[ii][0] = 1.0;
            for (int p = 1; p <= try_pd; ++p)
                col_x_pow[ii][p] = col_x_pow[ii][p-1] * col_x[ii];
        }
        for (int jj = 0; jj < ncpy; ++jj) {
            col_y_pow[jj][0] = 1.0;
            for (int p = 1; p <= try_pd; ++p)
                col_y_pow[jj][p] = col_y_pow[jj][p-1] * col_y[jj];
        }

        for (int k_basis = 0; k_basis < npb; ++k_basis) {
            int a = dx_arr[k_basis];
            int b = dy_arr[k_basis];
            for (int jj = 0; jj < ncpy; ++jj) {
                double yp = col_y_pow[jj][b];
                for (int ii = 0; ii < ncpx; ++ii) {
                    int col = jj * ncpx + ii;
                    S->poly_emat[k_basis * npts_grid + col]
                        = col_x_pow[ii][a] * yp;
                }
            }
        }

        /* rhs = poly_emat^T * coef_p_xy */
        for (int i = 0; i < npts_grid; ++i) {
            for (int c = 0; c < ncomp; ++c) {
                double sum = 0.0;
                for (int p = 0; p < npb; ++p)
                    sum += S->poly_emat[p*npts_grid + i] * S->coef_p_xy[p*ncomp + c];
                S->rhs[i*ncomp + c] = sum;
            }
        }

        /* bbasis^T */
        for (int i = 0; i < npts_grid; ++i)
            for (int j = 0; j < nl2d; ++j)
                S->bbasis_T[i*nl2d + j] = S->bbasis[j*npts_grid + i];

        /* solve bbasis^T * coeffs_BB = rhs */
        rc = qi_la_lu_solve(npts_grid, ncomp, S->bbasis_T, S->rhs, S->piv);
        if (rc != 0) continue;

        int ind_x = k_x - mux;
        int ind_y = k_y - muy;
        int row_idx = ind_x * nl_y + ind_y;
        for (int c = 0; c < ncomp; ++c)
            out_coeff[c] = S->rhs[row_idx*ncomp + c];

        chosen_pd = try_pd;
        chosen_sigmamin = sigma_min;
        break;
    }

    if (chosen_pd < 0) return QI_ERR_NO_ADMISSIBLE_DEG;
    if (deg_used) *deg_used = chosen_pd;
    if (cond_out) *cond_out = (chosen_sigmamin > 0.0) ? (1.0/chosen_sigmamin) : 0.0;
    return QI_OK;
}

/* ------------------------------------------------------------------------ */
/* Main entry point                                                         */
/* ------------------------------------------------------------------------ */
qi_status_t qi_compute_coeffs(const qi_hspace_t  *hs,
                              const qi_data_t    *dat,
                              const qi_options_t *opts_in,
                              double             *coeffs,
                              qi_diag_t          *diag)
{
    if (!hs || !dat || !coeffs) return QI_ERR_INVALID_ARG;
    if (hs->nlevels < 1)        return QI_ERR_INVALID_ARG;
    if (hs->ndof < 1)           return QI_ERR_INVALID_ARG;
    if (dat->npoints < 1 || dat->ncomp < 1) return QI_ERR_INVALID_ARG;
    if (hs->degree_x < 0 || hs->degree_x > QI_MAX_DEGREE) return QI_ERR_INVALID_ARG;
    if (hs->degree_y < 0 || hs->degree_y > QI_MAX_DEGREE) return QI_ERR_INVALID_ARG;

    qi_options_t opts;
    qi_options_defaults(&opts);
    if (opts_in) opts = *opts_in;
    int max_degree = (opts.max_degree > 0)
                     ? opts.max_degree : hs->degree_x;
    if (max_degree > QI_MAX_DEGREE) max_degree = QI_MAX_DEGREE;
    double sigma_threshold = (opts.sigma_threshold > 0.0)
                             ? opts.sigma_threshold : 1.0e2;

    const int *level_offset = hs->active_offset;

#ifdef _OPENMP
    if (opts.n_threads > 0) omp_set_num_threads(opts.n_threads);
#endif

    /* ---- Optimisation #1: build spatial grid index ---- */
    qi_spatial_t spatial;
    {
        int finest_nel = 1;
        for (int lv = 0; lv < hs->nlevels; ++lv) {
            int n = hs->nel_x[lv] > hs->nel_y[lv] ? hs->nel_x[lv] : hs->nel_y[lv];
            if (n > finest_nel) finest_nel = n;
        }
        qi_spatial_build_2d(&spatial, dat->npoints, dat->x, dat->y, finest_nel);
    }

    /* ---- Optimisation #3: pre-compute knot vectors per level ---- */
    qi_level_knots_t *lev_knots = (qi_level_knots_t *)malloc(
        (size_t)hs->nlevels * sizeof(qi_level_knots_t));
    if (!lev_knots) { qi_spatial_free(&spatial); return QI_ERR_OUT_OF_MEMORY; }
    for (int lv = 0; lv < hs->nlevels; ++lv) {
        int nkx = hs->nel_x[lv] + 2*hs->degree_x + 1;
        int nky = hs->nel_y[lv] + 2*hs->degree_y + 1;
        lev_knots[lv].Ux = (double *)malloc((size_t)nkx * sizeof(double));
        lev_knots[lv].Uy = (double *)malloc((size_t)nky * sizeof(double));
        qi_make_clamped_knots(hs, lv, 0, lev_knots[lv].Ux, &lev_knots[lv].nUx);
        qi_make_clamped_knots(hs, lv, 1, lev_knots[lv].Uy, &lev_knots[lv].nUy);
    }

    int err_count = 0;

    #pragma omp parallel reduction(+:err_count)
    {
        qi_scratch_t S;
        qi_scratch_init(&S);

        int k;
        #pragma omp for schedule(dynamic, 16)
        for (k = 0; k < hs->ndof; ++k) {
            double  c_buf[64];
            double *c_use = (dat->ncomp <= 64) ? c_buf
                            : (double*)malloc((size_t)dat->ncomp*sizeof(double));
            int deg_used_k = -1, nloc = 0, level = -1;
            double cond  = 0.0;
            int rc = qi_solve_one_basis(hs, dat, &opts, level_offset,
                                        &spatial, lev_knots,
                                        k, max_degree, sigma_threshold,
                                        &S, c_use, &deg_used_k,
                                        &nloc, &cond, &level);
            if (rc != QI_OK) {
                err_count += 1;
                for (int c = 0; c < dat->ncomp; ++c) c_use[c] = 0.0;
            }
            for (int c = 0; c < dat->ncomp; ++c)
                coeffs[k*dat->ncomp + c] = c_use[c];
            if (diag) {
                if (diag->deg_used)      diag->deg_used[k]      = deg_used_k;
                if (diag->nlocal)        diag->nlocal[k]        = nloc;
                if (diag->cond_estimate) diag->cond_estimate[k] = cond;
                if (diag->level)         diag->level[k]         = level;
            }
            if (c_use != c_buf) free(c_use);
        }
        qi_scratch_free(&S);
    }

    /* cleanup pre-computed data */
    for (int lv = 0; lv < hs->nlevels; ++lv) {
        free(lev_knots[lv].Ux);
        free(lev_knots[lv].Uy);
    }
    free(lev_knots);
    qi_spatial_free(&spatial);

    if (err_count > 0 && opts.verbose) {
        fprintf(stderr, "[libqi] warning: %d basis function(s) failed\n", err_count);
    }
    return (err_count == 0) ? QI_OK : QI_ERR_NO_ADMISSIBLE_DEG;
}

/* ------------------------------------------------------------------------ */
/* Single-level evaluation (optimisation #7: factor Nx*Ny)                  */
/* ------------------------------------------------------------------------ */
qi_status_t qi_eval_tp(const qi_hspace_t *hs, int level,
                       const double *coeffs_lev, int ncomp,
                       int npts, const double *xs, const double *ys,
                       double *out)
{
    if (!hs || !coeffs_lev || !xs || !ys || !out) return QI_ERR_INVALID_ARG;
    if (level < 0 || level >= hs->nlevels) return QI_ERR_INVALID_ARG;

    int dx = hs->degree_x;
    int dy = hs->degree_y;
    int nx_el = hs->nel_x[level];
    int ny_el = hs->nel_y[level];
    int ndof_dir_x = nx_el + dx;
    int ndof_dir_y = ny_el + dy;

    double *Ux = (double *)malloc((size_t)(nx_el + 2*dx + 1)*sizeof(double));
    double *Uy = (double *)malloc((size_t)(ny_el + 2*dy + 1)*sizeof(double));
    if (!Ux || !Uy) { free(Ux); free(Uy); return QI_ERR_OUT_OF_MEMORY; }
    int nUx, nUy;
    qi_make_clamped_knots(hs, level, 0, Ux, &nUx);
    qi_make_clamped_knots(hs, level, 1, Uy, &nUy);

    /* max (p+1)^2 non-zero basis products per point */
    int max_nz = (dx + 1) * (dy + 1);

    for (int i = 0; i < npts; ++i) {
        double u = xs[i], v = ys[i];
        double Nx[16], Ny[16];
        int spx = qi_bs_findspan(ndof_dir_x - 1, dx, u, Ux);
        int spy = qi_bs_findspan(ndof_dir_y - 1, dy, v, Uy);
        qi_bs_basisfun(spx, u, dx, Ux, Nx);
        qi_bs_basisfun(spy, v, dy, Uy, Ny);

        /* Pre-compute NxNy products and coefficient offsets (opt #7) */
        double NxNy[81];       /* max (4+1)^2 = 25, padded */
        int    coeff_off[81];
        int nz = 0;
        for (int rx = 0; rx <= dx; ++rx) {
            int kx = spx - dx + rx;
            if (kx < 0 || kx >= ndof_dir_x) continue;
            for (int ry = 0; ry <= dy; ++ry) {
                int ky = spy - dy + ry;
                if (ky < 0 || ky >= ndof_dir_y) continue;
                NxNy[nz]      = Nx[rx] * Ny[ry];
                coeff_off[nz]  = (ky * ndof_dir_x + kx) * ncomp;
                nz++;
            }
        }

        for (int c = 0; c < ncomp; ++c) {
            double s = 0.0;
            for (int j = 0; j < nz; ++j)
                s += coeffs_lev[coeff_off[j] + c] * NxNy[j];
            out[i*ncomp + c] = s;
        }
    }
    free(Ux); free(Uy);
    return QI_OK;
}

/* ------------------------------------------------------------------------ */
/* Subdivision (convert hierarchical coeffs to finest level)                */
/* ------------------------------------------------------------------------ */
static int dyadic_log2(int x)
{
    if (x <= 0) return -1;
    int k = 0;
    while ((1 << k) < x) ++k;
    return ((1 << k) == x) ? k : -1;
}

static int qi_subdiv_matrix_1D(int p, int nel_coarse, int nel_fine, double *M)
{
    int ncoarse = nel_coarse + p;
    int nfine   = nel_fine + p;
    double *Uc = (double *)malloc((size_t)(nel_coarse + 2*p + 1)*sizeof(double));
    double *Uf = (double *)malloc((size_t)(nel_fine   + 2*p + 1)*sizeof(double));
    if (!Uc || !Uf) { free(Uc); free(Uf); return -1; }
    for (int j = 0; j <= p; ++j) Uc[j] = 0.0;
    for (int j = 1; j < nel_coarse; ++j) Uc[p + j] = (double)j / (double)nel_coarse;
    for (int j = 0; j <= p; ++j) Uc[p + nel_coarse + j] = 1.0;
    for (int j = 0; j <= p; ++j) Uf[j] = 0.0;
    for (int j = 1; j < nel_fine; ++j) Uf[p + j] = (double)j / (double)nel_fine;
    for (int j = 0; j <= p; ++j) Uf[p + nel_fine + j] = 1.0;

    double *gp = (double *)malloc((size_t)nfine * sizeof(double));
    for (int i = 0; i < nfine; ++i) {
        double s = 0.0;
        for (int k = 1; k <= p; ++k) s += Uf[i + k];
        gp[i] = s / (double)p;
    }

    for (int i = 0; i < nfine*ncoarse; ++i) M[i] = 0.0;
    double Nbuf[16];
    for (int i = 0; i < nfine; ++i) {
        double u = gp[i];
        int span = qi_bs_findspan(ncoarse - 1, p, u, Uc);
        qi_bs_basisfun(span, u, p, Uc, Nbuf);
        for (int r = 0; r <= p; ++r) {
            int j = span - p + r;
            if (j >= 0 && j < ncoarse)
                M[i*ncoarse + j] = Nbuf[r];
        }
    }
    free(Uc); free(Uf); free(gp);
    return 0;
}

qi_status_t qi_to_finest(const qi_hspace_t *hs, const double *coeffs,
                         int ncomp, double *coeffs_finest)
{
    if (!hs || !coeffs || !coeffs_finest) return QI_ERR_INVALID_ARG;

    int M = hs->nlevels;
    int dx = hs->degree_x, dy = hs->degree_y;
    int finest_nx = hs->nel_x[M-1];
    int finest_ny = hs->nel_y[M-1];
    int ndof_x_fin = finest_nx + dx;
    int ndof_y_fin = finest_ny + dy;
    int ndof_fin   = ndof_x_fin * ndof_y_fin;

    for (int i = 0; i < ndof_fin*ncomp; ++i) coeffs_finest[i] = 0.0;

    for (int lev = 0; lev < M; ++lev) {
        int nx_l = hs->nel_x[lev], ny_l = hs->nel_y[lev];
        int ndof_x_l = nx_l + dx, ndof_y_l = ny_l + dy;

        int rx = (lev == M-1) ? 1 : (finest_nx / nx_l);
        int ry = (lev == M-1) ? 1 : (finest_ny / ny_l);
        if (lev != M-1 && (rx*nx_l != finest_nx || ry*ny_l != finest_ny))
            return QI_ERR_INVALID_ARG;

        double *Mx = (double *)malloc((size_t)ndof_x_fin*ndof_x_l*sizeof(double));
        double *My = (double *)malloc((size_t)ndof_y_fin*ndof_y_l*sizeof(double));
        if (!Mx || !My) { free(Mx); free(My); return QI_ERR_OUT_OF_MEMORY; }
        if (qi_subdiv_matrix_1D(dx, nx_l, finest_nx, Mx) != 0)
            { free(Mx); free(My); return QI_ERR_OUT_OF_MEMORY; }
        if (qi_subdiv_matrix_1D(dy, ny_l, finest_ny, My) != 0)
            { free(Mx); free(My); return QI_ERR_OUT_OF_MEMORY; }

        int n_active = hs->active_offset[lev+1] - hs->active_offset[lev];
        const int *active = hs->active_indices + hs->active_offset[lev];
        for (int a = 0; a < n_active; ++a) {
            int idx_l = active[a];
            int kx = idx_l % ndof_x_l;
            int ky = idx_l / ndof_x_l;
            int k_global = hs->active_offset[lev] + a;
            const double *coef = coeffs + k_global*ncomp;

            for (int Kx = 0; Kx < ndof_x_fin; ++Kx) {
                double w_x = Mx[Kx*ndof_x_l + kx];
                if (w_x == 0.0) continue;
                for (int Ky = 0; Ky < ndof_y_fin; ++Ky) {
                    double w_y = My[Ky*ndof_y_l + ky];
                    if (w_y == 0.0) continue;
                    double w = w_x * w_y;
                    int F = Ky*ndof_x_fin + Kx;
                    for (int c = 0; c < ncomp; ++c)
                        coeffs_finest[F*ncomp + c] += w * coef[c];
                }
            }
        }
        free(Mx); free(My);
    }
    return QI_OK;
}

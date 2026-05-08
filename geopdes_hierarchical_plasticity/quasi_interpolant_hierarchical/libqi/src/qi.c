/*
 * qi.c — Main implementation of the THB-spline quasi-interpolant.
 *
 * Contains the parallel loop over active basis functions and the per-basis
 * computation (support box, local data selection, polynomial LS,
 * polynomial-to-B-spline basis change, coefficient extraction).
 */
#include "../include/qi.h"
#include "qi_linalg.h"
#include "qi_bspline.h"

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
    return "libqi 1.0.0  (Hierarchical Quasi-Interpolant; OpenMP)";
}

void qi_options_defaults(qi_options_t *opts)
{
    opts->max_degree      = 0;       /* meaning "use degree_x" */
    opts->sigma_threshold = 0.0;     /* meaning "use 1e2"      */
    opts->n_threads       = 0;
    opts->verbose         = 0;
}

/* ------------------------------------------------------------------------ */
/* helpers                                                                  */
/* ------------------------------------------------------------------------ */

/* Locate the level lev (0-based) of the k-th THB-spline (0-based) using
   ndof_per_level via a precomputed cumulative offset table.
   Returns the local index within the level via *kloc_out.            */
static inline int qi_locate_level(const int *ndof_offset, int nlevels,
                                  int k, int *kloc_out)
{
    /* binary search on monotone ndof_offset[0..nlevels] (== active_offset) */
    int lo = 0, hi = nlevels - 1;
    while (lo < hi) {
        int mid = (lo + hi + 1) / 2;
        if (ndof_offset[mid] <= k) lo = mid;
        else                       hi = mid - 1;
    }
    *kloc_out = k - ndof_offset[lo];
    return lo;
}

/* breaks_x[lev][i] for i in 0..nel_x[lev], possibly uniform */
static inline double qi_break_x(const qi_hspace_t *h, int lev, int i)
{
    if (h->breaks_x_data) {
        return h->breaks_x_data[h->breaks_x_offset[lev] + i];
    }
    return (double)i / (double)h->nel_x[lev];
}
static inline double qi_break_y(const qi_hspace_t *h, int lev, int i)
{
    if (h->breaks_y_data) {
        return h->breaks_y_data[h->breaks_y_offset[lev] + i];
    }
    return (double)i / (double)h->nel_y[lev];
}

/* Construct the open-uniform clamped knot vector at level `lev` for direction
   x (axis==0) or y (axis==1).  Length = nel + 2*p + 1, contains 0 (p+1 times),
   then breakpoints, then 1 (p+1 times).  If non-uniform breaks are provided
   they are used; otherwise i/nel are used. Caller provides buffer.        */
static void qi_make_clamped_knots(const qi_hspace_t *h, int lev, int axis,
                                  double *U_out, int *nknots_out)
{
    int p   = (axis == 0) ? h->degree_x : h->degree_y;
    int nel = (axis == 0) ? h->nel_x[lev] : h->nel_y[lev];
    int nb  = nel + 1; /* number of breakpoints */
    /* total knots = (p+1) at left + (nb - 2 interior + 2 endpoints with mult 1)
       ... wait, easier: clamped knot vector = each end with mult p+1, others
       with mult 1, total = nb + 2*p. */
    int nknots = nb + 2*p;
    int idx = 0;
    /* left end, p+1 times, value = first break (= 0) */
    double left = (axis == 0) ? qi_break_x(h, lev, 0) : qi_break_y(h, lev, 0);
    for (int j = 0; j <= p; ++j) U_out[idx++] = left;
    /* middle breaks */
    for (int j = 1; j < nb - 1; ++j) {
        U_out[idx++] = (axis == 0) ? qi_break_x(h, lev, j) : qi_break_y(h, lev, j);
    }
    /* right end, p+1 times */
    double right = (axis == 0) ? qi_break_x(h, lev, nb-1) : qi_break_y(h, lev, nb-1);
    for (int j = 0; j <= p; ++j) U_out[idx++] = right;
    *nknots_out = nknots;
}

/* enumerate (deg_x, deg_y) pairs with deg_x + deg_y <= pd in lex order
 * (deg_x outer, deg_y inner).  Returns number of pairs.                   */
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

/* binomial coefficient C(n,k); n,k small */
static long long binomial_int(int n, int k)
{
    if (k < 0 || k > n) return 0;
    if (k == 0 || k == n) return 1;
    if (k > n - k) k = n - k;
    long long r = 1;
    for (int i = 1; i <= k; ++i) {
        r = r * (n - i + 1) / i;
    }
    return r;
}

/* ------------------------------------------------------------------------ */
/* Per-basis-function computation                                           */
/* ------------------------------------------------------------------------ */

/* Constants for buffer sizing */
#define QI_MAX_DEGREE 4
/* (max_pd+1)(max_pd+2)/2 for max_pd=4 = 15 */
#define QI_MAX_NPB    ((QI_MAX_DEGREE+1)*(QI_MAX_DEGREE+2)/2)
/* local 1D B-spline count <= 2*p+1 */
#define QI_MAX_NL1D   (2*QI_MAX_DEGREE + 1)
/* local 2D basis count = local-1D-x * local-1D-y */
#define QI_MAX_NL2D   (QI_MAX_NL1D*QI_MAX_NL1D)

/* per-thread scratch buffer */
typedef struct {
    /* full clamped knot vectors (buffer sized to be safe) */
    double *Ux;   int Ux_cap;
    double *Uy;   int Uy_cap;

    /* indices of local data inside the disk */
    int    *loc_idx;  int loc_idx_cap;

    /* local data extracted */
    double *xloc;     int xloc_cap;
    double *yloc;     int yloc_cap;
    double *floc;     int floc_cap;       /* counts doubles, == nloc*ncomp */
    int     loc_cap;  /* legacy; kept zero */

    /* design matrix and pinv-solve scratch */
    double *A;        int A_cap;
    double *coef_p;   int coef_p_cap;        /* npb x ncomp */
    double *coef_p_xy;int coef_p_xy_cap;     /* coeffs in (x,y)-monomial basis */

    /* B-spline 1D collocation matrices (full, ndof_dir x ncpx)  */
    double *aux_x; int aux_x_cap;
    double *aux_y; int aux_y_cap;

    /* basis change matrices */
    double *bbasis;     int bbasis_cap;     /* nl2d x npts_grid */
    double *poly_emat;  int poly_emat_cap;  /* npb  x npts_grid */
    double *bbasis_T;   int bbasis_T_cap;   /* npts x nl2d   */
    double *rhs;        int rhs_cap;        /* npts x ncomp  */
    int    *piv;        int piv_cap;

    /* pinv work */
    double *work;       int work_cap;
} qi_scratch_t;

static void qi_scratch_init(qi_scratch_t *s)
{
    memset(s, 0, sizeof(*s));
}
static void *qi_xrealloc(void *p, int cap_old, int cap_new, size_t elem)
{
    if (cap_new <= cap_old) return p;
    void *q = realloc(p, (size_t)cap_new * elem);
    return q;
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
    free(s->Ux); free(s->Uy);
    free(s->loc_idx);
    free(s->xloc); free(s->yloc); free(s->floc);
    free(s->A);
    free(s->coef_p);
    free(s->coef_p_xy);
    free(s->aux_x); free(s->aux_y);
    free(s->bbasis); free(s->poly_emat);
    free(s->bbasis_T);
    free(s->rhs);
    free(s->piv);
    free(s->work);
}

/* ------------------------------------------------------------------------ */
/* Solve the QI for one basis function k.                                   */
/* ------------------------------------------------------------------------ */
static int qi_solve_one_basis(const qi_hspace_t  *hs,
                              const qi_data_t    *dat,
                              const qi_options_t *opts,
                              const int          *level_offset,
                              int                 k,
                              int                 max_degree,
                              double              sigma_threshold,
                              qi_scratch_t       *S,
                              double             *out_coeff,   /* ncomp */
                              int                *deg_used,
                              int                *nlocal_out,
                              double             *cond_out,
                              int                *level_out)
{
    const int dx = hs->degree_x;
    const int dy = hs->degree_y;
    const int ncomp = dat->ncomp;

    int kloc;
    int lev = qi_locate_level(level_offset, hs->nlevels, k, &kloc);
    if (level_out) *level_out = lev;

    /* active linear index in the level-lev tensor product space */
    int act_idx = hs->active_indices[level_offset[lev] + kloc];
    int ndof_dir_x = hs->nel_x[lev] + dx;
    int ndof_dir_y = hs->nel_y[lev] + dy;
    int k_x = act_idx % ndof_dir_x;          /* MATLAB column-major */
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

    /* select local data inside disk of radius r0 around centre.
       If too few, enlarge by integer multiples up to ratio_max.            */
    double cx = 0.5 * (x_mu + x_nu);
    double cy = 0.5 * (y_mu + y_nu);
    double r0 = 0.5 * sqrt(Dx*Dx + Dy*Dy);

    /* ratioray (max enlarging factor) — same expression as MATLAB code */
    int ref_lev = (lev >= 1) ? (lev - 1) : lev;
    double ref_nx = (double)hs->nel_x[ref_lev];
    double ref_ny = (double)hs->nel_y[ref_lev];
    double ratio_d = sqrt(((double)(dx+1)/ref_nx)*((double)(dx+1)/ref_nx)
                        + ((double)(dy+1)/ref_ny)*((double)(dy+1)/ref_ny))
                     / (2.0 * r0);
    int ratio_max = (int)ceil(ratio_d);
    if (ratio_max < 1) ratio_max = 1;
    int max_iter = 2*ratio_max + 1;

    /* Target initial degree */
    int pd = max_degree;
    int npb_target = (pd+1)*(pd+2)/2;

    /* allocate / grow loc_idx */
    QI_ENSURE(S->loc_idx, S->loc_idx_cap, dat->npoints, sizeof(int));

    int nloc = 0;
    int iter = 1;
    double r = r0;
    while (1) {
        /* recompute neighbour set with current r */
        nloc = 0;
        double r2 = r * r;
        for (int i = 0; i < dat->npoints; ++i) {
            double ddx = dat->x[i] - cx;
            double ddy = dat->y[i] - cy;
            if (ddx*ddx + ddy*ddy <= r2) S->loc_idx[nloc++] = i;
        }
        /* reduce pd if not enough data, but never below 0 */
        while (nloc < (pd+1)*(pd+2)/2 && pd > 0) --pd;
        npb_target = (pd+1)*(pd+2)/2;
        if (nloc >= npb_target) break;
        if (iter > max_iter) break;
        ++iter;
        r = (double)iter * r0;
    }

    if (nloc < npb_target) {
        /* try one final fallback: take all data points for this basis */
        nloc = dat->npoints;
        for (int i = 0; i < nloc; ++i) S->loc_idx[i] = i;
        while (nloc < (pd+1)*(pd+2)/2 && pd > 0) --pd;
        npb_target = (pd+1)*(pd+2)/2;
        if (nloc < npb_target) return QI_ERR_NOT_ENOUGH_DATA;
    }

    if (nlocal_out) *nlocal_out = nloc;

    /* ensure local-data buffers (each grown independently and consistently) */
    QI_ENSURE(S->xloc, S->xloc_cap, nloc,        sizeof(double));
    QI_ENSURE(S->yloc, S->yloc_cap, nloc,        sizeof(double));
    QI_ENSURE(S->floc, S->floc_cap, nloc*ncomp,  sizeof(double));

    /* extract local data, scaled to [0,1]^2 (datascal in MATLAB code) */
    for (int i = 0; i < nloc; ++i) {
        int idx = S->loc_idx[i];
        S->xloc[i] = (dat->x[idx] - x_mu) / Dx;
        S->yloc[i] = (dat->y[idx] - y_mu) / Dy;
        for (int c = 0; c < ncomp; ++c)
            S->floc[i*ncomp + c] = dat->f[idx*ncomp + c];
    }

    /* knot vectors at level lev */
    int Ux_need = nx + 2*dx + 1;
    int Uy_need = ny + 2*dy + 1;
    QI_ENSURE(S->Ux, S->Ux_cap, Ux_need, sizeof(double));
    QI_ENSURE(S->Uy, S->Uy_cap, Uy_need, sizeof(double));
    int nUx, nUy;
    qi_make_clamped_knots(hs, lev, 0, S->Ux, &nUx);
    qi_make_clamped_knots(hs, lev, 1, S->Uy, &nUy);

    /* Local 2D B-spline tile */
    int nl_x = nux - mux + dx;     /* number of 1D B-splines in x */
    int nl_y = nuy - muy + dy;     /* number of 1D B-splines in y */
    int ncpx = nl_x;               /* collocation points = number of 1D splines */
    int ncpy = nl_y;
    int npts_grid = ncpx * ncpy;
    int nl2d  = nl_x * nl_y;

    /* The local 2D basis spans the same number of points as collocation
       points, so the basis-change matrix is square: nl2d x nl2d.            */

    /* Build collocation points: linspace(x_mu, x_nu, ncpx) etc.            */
    double col_x[64], col_y[64];
    if (ncpx >= 64 || ncpy >= 64) return QI_ERR_INVALID_ARG;
    if (ncpx == 1) col_x[0] = 0.5*(x_mu + x_nu);
    else for (int i = 0; i < ncpx; ++i)
        col_x[i] = x_mu + (x_nu - x_mu) * (double)i / (double)(ncpx - 1);
    if (ncpy == 1) col_y[0] = 0.5*(y_mu + y_nu);
    else for (int i = 0; i < ncpy; ++i)
        col_y[i] = y_mu + (y_nu - y_mu) * (double)i / (double)(ncpy - 1);

    /* spcol of the full level-knot in x at col_x  -> aux_x has shape
       (ndof_dir_x, ncpx), row-major.  Then we extract rows mux..mux+nl_x-1.*/
    int aux_x_need = ndof_dir_x * ncpx;
    int aux_y_need = ndof_dir_y * ncpy;
    QI_ENSURE(S->aux_x, S->aux_x_cap, aux_x_need, sizeof(double));
    QI_ENSURE(S->aux_y, S->aux_y_cap, aux_y_need, sizeof(double));
    qi_bs_spcol(S->Ux, nUx, dx, col_x, ncpx, S->aux_x);
    qi_bs_spcol(S->Uy, nUy, dy, col_y, ncpy, S->aux_y);

    /* bbasis matrix (nl2d x npts_grid).  Layout (h_loc * nl_y + q_loc, j*ncpx+i).
       h_loc indexes 1D-x basis (h = mux+h_loc), q_loc indexes 1D-y basis.
       Column j*ncpx+i corresponds to col_y[j], col_x[i] (matching MATLAB
       column-major flatten of aux = aux_x * aux_y^T).                      */
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
                    /* column index: jj * ncpx + ii  (MATLAB column-major) */
                    int col = jj * ncpx + ii;
                    S->bbasis[row * npts_grid + col]
                        = S->aux_x[h * ncpx + ii] * S->aux_y[q * ncpy + jj];
                }
        }
    }

    /* The solving loop over polynomial degree (descending) — we accept the
       FIRST admissible pd, since that is the maximum admissible degree.    */
    int dx_arr[QI_MAX_NPB], dy_arr[QI_MAX_NPB];

    int    chosen_pd       = -1;
    double chosen_sigmamin = 0.0;

    /* `out_coeff` accumulates the chosen coefficients */

    for (int try_pd = pd; try_pd >= 0; --try_pd) {
        int npb = enumerate_monomials(try_pd, dx_arr, dy_arr);

        if (nloc < npb) continue;

        /* Build A (nloc x npb) where A[i,k] = u_i^{dx_arr[k]} * v_i^{dy_arr[k]} */
        QI_ENSURE(S->A, S->A_cap, nloc * npb, sizeof(double));
        for (int i = 0; i < nloc; ++i) {
            double u = S->xloc[i];
            double v = S->yloc[i];
            /* precompute powers of u, v */
            double upow[QI_MAX_NPB];
            double vpow[QI_MAX_NPB];
            upow[0] = 1.0; vpow[0] = 1.0;
            for (int p = 1; p <= try_pd; ++p) {
                upow[p] = upow[p-1] * u;
                vpow[p] = vpow[p-1] * v;
            }
            for (int p = 0; p < npb; ++p) {
                S->A[i*npb + p] = upow[dx_arr[p]] * vpow[dy_arr[p]];
            }
        }

        /* allocate coef_p */
        QI_ENSURE(S->coef_p, S->coef_p_cap, npb * ncomp, sizeof(double));
        QI_ENSURE(S->coef_p_xy, S->coef_p_xy_cap, npb * ncomp, sizeof(double));

        /* pinv-solve work buffer */
        int wneed = 2*npb*npb + npb + 2*npb*ncomp;
        QI_ENSURE(S->work, S->work_cap, wneed, sizeof(double));

        /* preserve A for reuse below — qi_la_pinv_solve doesn't destroy A */
        double sigma_min = 0.0;
        int rc = qi_la_pinv_solve(nloc, npb, ncomp, S->A, S->floc,
                                  S->coef_p, &sigma_min,
                                  S->work, S->work_cap);
        if (rc != 0) continue;

        /* MSV check: MSV2 = 1/sigma_min(A); admissible iff MSV2 < threshold */
        double MSV2 = (sigma_min > 0.0) ? (1.0 / sigma_min) : 1e300;
        if (MSV2 >= sigma_threshold) continue;

        /* Convert polynomial coefficients from scaled (u,v)-monomial basis
           to (x,y)-monomial basis.   p(u,v) = sum_k c_k u^{a_k} v^{b_k}.
           u = (x - x_mu) / Dx, v = (y - y_muy) / Dy.
                 (x - x_mu)^a = sum_{r=0..a} C(a,r) (-x_mu)^(a-r) x^r
           So p(x,y) = sum_k c_k / (Dx^a_k Dy^b_k)
                       sum_{r=0..a_k} C(a_k,r) (-x_mu)^(a_k-r) x^r
                       sum_{s=0..b_k} C(b_k,s) (-y_mu)^(b_k-s) y^s.

           Coefficient of x^r y^s in p(x,y) for (r,s) with r+s <= try_pd is:
              d_{rs} = sum_{(a,b): a>=r,b>=s,a+b<=try_pd}
                          c_{(a,b)} / (Dx^a Dy^b)
                          C(a,r) (-x_mu)^(a-r) C(b,s) (-y_mu)^(b-s)
        */
        for (int idx_xy = 0; idx_xy < npb; ++idx_xy) {
            int r = dx_arr[idx_xy];
            int s = dy_arr[idx_xy];
            for (int c = 0; c < ncomp; ++c) {
                double sum = 0.0;
                for (int idx_uv = 0; idx_uv < npb; ++idx_uv) {
                    int a = dx_arr[idx_uv];
                    int b = dy_arr[idx_uv];
                    if (a < r || b < s) continue;
                    /* term = c_{ab} / (Dx^a Dy^b) * C(a,r) C(b,s)
                                   * (-x_mu)^(a-r) * (-y_mu)^(b-s) */
                    double t = S->coef_p[idx_uv*ncomp + c];
                    /* divide */
                    double dxa = 1.0, dyb = 1.0;
                    for (int p = 0; p < a; ++p) dxa *= Dx;
                    for (int p = 0; p < b; ++p) dyb *= Dy;
                    t /= (dxa * dyb);
                    /* binomials and powers */
                    long long cab = binomial_int(a, r);
                    long long cbs = binomial_int(b, s);
                    double xpow = 1.0, ypow = 1.0;
                    for (int p = 0; p < a - r; ++p) xpow *= -x_mu;
                    for (int p = 0; p < b - s; ++p) ypow *= -y_mu;
                    sum += t * (double)cab * (double)cbs * xpow * ypow;
                }
                S->coef_p_xy[idx_xy*ncomp + c] = sum;
            }
        }

        /* Build poly_emat (npb x npts_grid): row k is x^{dx_arr[k]} y^{dy_arr[k]}
           evaluated on the (i,j) collocation grid, column = j*ncpx + i.    */
        for (int k_basis = 0; k_basis < npb; ++k_basis) {
            int a = dx_arr[k_basis];
            int b = dy_arr[k_basis];
            for (int jj = 0; jj < ncpy; ++jj) {
                double yv = col_y[jj];
                double yp = 1.0;
                for (int p = 0; p < b; ++p) yp *= yv;
                for (int ii = 0; ii < ncpx; ++ii) {
                    double xv = col_x[ii];
                    double xp = 1.0;
                    for (int p = 0; p < a; ++p) xp *= xv;
                    int col = jj * ncpx + ii;
                    S->poly_emat[k_basis * npts_grid + col] = xp * yp;
                }
            }
        }

        /* Compute rhs = poly_emat^T * coef_p_xy   (npts_grid x ncomp) */
        for (int i = 0; i < npts_grid; ++i) {
            for (int c = 0; c < ncomp; ++c) {
                double sum = 0.0;
                for (int p = 0; p < npb; ++p)
                    sum += S->poly_emat[p*npts_grid + i] * S->coef_p_xy[p*ncomp + c];
                S->rhs[i*ncomp + c] = sum;
            }
        }

        /* Build bbasis^T (npts_grid x nl2d) */
        for (int i = 0; i < npts_grid; ++i)
            for (int j = 0; j < nl2d; ++j)
                S->bbasis_T[i*nl2d + j] = S->bbasis[j*npts_grid + i];

        /* solve bbasis^T * coeffs_BB = rhs  (square, npts_grid == nl2d) */
        rc = qi_la_lu_solve(npts_grid, ncomp, S->bbasis_T, S->rhs, S->piv);
        if (rc != 0) continue;

        /* extract row corresponding to (k_x, k_y) */
        int ind_x = k_x - mux;     /* 0-based local */
        int ind_y = k_y - muy;
        int row_idx = ind_x * nl_y + ind_y;
        for (int c = 0; c < ncomp; ++c)
            out_coeff[c] = S->rhs[row_idx*ncomp + c];

        chosen_pd = try_pd;
        chosen_sigmamin = sigma_min;
        break;   /* first admissible == max admissible */
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

    /* level offsets (== active_offset, but we recompute for safety) */
    /* MATLAB convention is that active_offset[lev+1] = active_offset[lev]
       + ndof_per_level[lev], but we accept whichever the caller passes.  */
    const int *level_offset = hs->active_offset;

#ifdef _OPENMP
    if (opts.n_threads > 0) omp_set_num_threads(opts.n_threads);
#endif

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
            int deg_used = -1, nloc = 0, level = -1;
            double cond  = 0.0;
            int rc = qi_solve_one_basis(hs, dat, &opts, level_offset,
                                        k, max_degree, sigma_threshold,
                                        &S, c_use, &deg_used,
                                        &nloc, &cond, &level);
            if (rc != QI_OK) {
                err_count += 1;
                /* on failure, fill with zeros */
                for (int c = 0; c < dat->ncomp; ++c) c_use[c] = 0.0;
            }
            for (int c = 0; c < dat->ncomp; ++c)
                coeffs[k*dat->ncomp + c] = c_use[c];
            if (diag) {
                if (diag->deg_used)      diag->deg_used[k]      = deg_used;
                if (diag->nlocal)        diag->nlocal[k]        = nloc;
                if (diag->cond_estimate) diag->cond_estimate[k] = cond;
                if (diag->level)         diag->level[k]         = level;
            }
            if (c_use != c_buf) free(c_use);
        }
        qi_scratch_free(&S);
    }

    if (err_count > 0 && opts.verbose) {
        fprintf(stderr, "[libqi] warning: %d basis function(s) failed\n", err_count);
    }
    return (err_count == 0) ? QI_OK : QI_ERR_NO_ADMISSIBLE_DEG;
}

/* ------------------------------------------------------------------------ */
/* Single-level evaluation                                                  */
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
    int nx = hs->nel_x[level];
    int ny = hs->nel_y[level];
    int ndof_dir_x = nx + dx;
    int ndof_dir_y = ny + dy;

    /* knot vectors */
    double *Ux = malloc((size_t)(nx + 2*dx + 1)*sizeof(double));
    double *Uy = malloc((size_t)(ny + 2*dy + 1)*sizeof(double));
    if (!Ux || !Uy) { free(Ux); free(Uy); return QI_ERR_OUT_OF_MEMORY; }
    int nUx, nUy;
    qi_make_clamped_knots(hs, level, 0, Ux, &nUx);
    qi_make_clamped_knots(hs, level, 1, Uy, &nUy);

    for (int i = 0; i < npts; ++i) {
        double u = xs[i], v = ys[i];
        double Nx[16], Ny[16];
        int spx = qi_bs_findspan(ndof_dir_x - 1, dx, u, Ux);
        int spy = qi_bs_findspan(ndof_dir_y - 1, dy, v, Uy);
        qi_bs_basisfun(spx, u, dx, Ux, Nx);
        qi_bs_basisfun(spy, v, dy, Uy, Ny);

        for (int c = 0; c < ncomp; ++c) {
            double s = 0.0;
            for (int rx = 0; rx <= dx; ++rx) {
                int kx = spx - dx + rx;
                if (kx < 0 || kx >= ndof_dir_x) continue;
                for (int ry = 0; ry <= dy; ++ry) {
                    int ky = spy - dy + ry;
                    if (ky < 0 || ky >= ndof_dir_y) continue;
                    int idx = ky * ndof_dir_x + kx;  /* col-major */
                    s += coeffs_lev[idx*ncomp + c] * Nx[rx] * Ny[ry];
                }
            }
            out[i*ncomp + c] = s;
        }
    }
    free(Ux); free(Uy);
    return QI_OK;
}

/* ------------------------------------------------------------------------ */
/* Subdivision (used to convert hierarchical coeffs to finest level)        */
/* ------------------------------------------------------------------------ */
/*
 *  We compute the matrix C : R^ndof -> R^ndof_finest such that
 *    sum_k λ_k T^{ℓ_k}_{J_k}(x) = sum_K c_K B^{M-1}_K(x)
 *  where T are THB-splines and B are tensor-product B-splines on the finest
 *  level.  We use the standard preservation-of-coefficients property: the
 *  truncation merely zeros out coefficients of finer-level B-splines whose
 *  support is entirely inside Ω^{ℓ+1}; for the simple case where the
 *  caller already keeps track of only active basis functions, the fine-level
 *  coefficient of a mother B-spline B^ℓ_J is its image under the standard
 *  knot-insertion matrix between level ℓ and the finest level (M-1).
 *
 *  For tensor-product knot insertion between two clamped uniform meshes
 *  whose numbers of cells satisfy nel_fine = r * nel_coarse, the coefficient
 *  vector transforms by a tensor-product subdivision matrix (Boehm/Lyche).
 *  We implement only the case nel_fine = 2^k * nel_coarse for simplicity —
 *  i.e. dyadic refinement, which is the typical use case for THB-splines.
 *
 *  If the meshes are not dyadic refinements the function returns
 *  QI_ERR_INVALID_ARG.  Use the GeoPDEs hspace_subdivision_matrix in MATLAB
 *  as a fallback.
 */
static int dyadic_log2(int x)
{
    if (x <= 0) return -1;
    int k = 0;
    while ((1 << k) < x) ++k;
    return ((1 << k) == x) ? k : -1;
}

/* Knot insertion: given (n+1) coefficients of a degree-p B-spline expansion
 * on a clamped uniform knot vector with nel cells, produce (2*nel+p)
 * coefficients on the dyadically refined mesh. We use the Lane-Riesenfeld
 * (averaging mask) formulation.
 *
 *   refined coefficients c_new[i] = sum_j  M[i, j] c_old[j]
 *
 * Implementation: explicit knot-insertion mask for uniform B-splines of
 * degree p, halving the spacing.  The mask is the up-sampling+convolution
 * with the "subdivision mask" b_p.
 *
 * For B-spline of degree p, dyadic subdivision mask is:
 *     b_p(k) = 2^{-p} * C(p+1, k)
 * applied as c_new[2i + r] = sum_k b_p(2i + r - 2j) c_old[j]   ... etc.
 *
 * For clamped (open uniform) knot vectors, boundary corrections are needed;
 * the simplest approach is to reuse the Oslo / discrete B-spline algorithm.
 * To keep this self-contained we simply build the global subdivision matrix
 * by *evaluating* the coarse B-splines on the fine knots — equivalent to
 * the Oslo algorithm.
 *
 * That's exactly what we do below.
 */
static int qi_subdiv_matrix_1D(int p, int nel_coarse, int nel_fine,
                               double *M)
{
    /* M has shape (nel_fine + p) x (nel_coarse + p), row-major. */
    int ncoarse = nel_coarse + p;
    int nfine   = nel_fine + p;

    /* coarse and fine knot vectors */
    double *Uc = malloc((size_t)(nel_coarse + 2*p + 1)*sizeof(double));
    double *Uf = malloc((size_t)(nel_fine   + 2*p + 1)*sizeof(double));
    if (!Uc || !Uf) { free(Uc); free(Uf); return -1; }
    /* uniform clamped */
    for (int j = 0; j <= p; ++j) Uc[j] = 0.0;
    for (int j = 1; j < nel_coarse; ++j) Uc[p + j] = (double)j / (double)nel_coarse;
    for (int j = 0; j <= p; ++j) Uc[p + nel_coarse + j] = 1.0;
    for (int j = 0; j <= p; ++j) Uf[j] = 0.0;
    for (int j = 1; j < nel_fine; ++j) Uf[p + j] = (double)j / (double)nel_fine;
    for (int j = 0; j <= p; ++j) Uf[p + nel_fine + j] = 1.0;

    /* Use the Greville abscissae of the fine basis as collocation points
       (gives a square invertible system on the fine side, and we then
       interpolate the coarse basis at those points — gives the
       subdivision rule).                                                  */
    double *gp = malloc((size_t)nfine * sizeof(double));
    for (int i = 0; i < nfine; ++i) {
        double s = 0.0;
        for (int k = 1; k <= p; ++k) s += Uf[i + k];
        gp[i] = s / (double)p;
    }

    /* M[i, j] = N_j^{coarse}(gp[i])  (since interp on Greville exactly
       reproduces the coefficients for fine spline of same degree).        */
    /* zero M */
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

    /* zero output */
    for (int i = 0; i < ndof_fin*ncomp; ++i) coeffs_finest[i] = 0.0;

    /* For each level, lift its contributions to finest using subdivision.    */
    for (int lev = 0; lev < M; ++lev) {
        int nx = hs->nel_x[lev], ny = hs->nel_y[lev];
        int ndof_x_l = nx + dx, ndof_y_l = ny + dy;

        /* check dyadic refinement */
        int rx = (lev == M-1) ? 1 : (finest_nx / nx);
        int ry = (lev == M-1) ? 1 : (finest_ny / ny);
        if (lev != M-1 && (rx*nx != finest_nx || ry*ny != finest_ny)) {
            return QI_ERR_INVALID_ARG;  /* non-dyadic */
        }

        /* build 1D subdivision matrices for x and y */
        double *Mx = malloc((size_t)ndof_x_fin*ndof_x_l*sizeof(double));
        double *My = malloc((size_t)ndof_y_fin*ndof_y_l*sizeof(double));
        if (!Mx || !My) { free(Mx); free(My); return QI_ERR_OUT_OF_MEMORY; }
        if (qi_subdiv_matrix_1D(dx, nx, finest_nx, Mx) != 0)
            { free(Mx); free(My); return QI_ERR_OUT_OF_MEMORY; }
        if (qi_subdiv_matrix_1D(dy, ny, finest_ny, My) != 0)
            { free(Mx); free(My); return QI_ERR_OUT_OF_MEMORY; }

        /* iterate over active basis at this level */
        int n_active = hs->active_offset[lev+1] - hs->active_offset[lev];
        const int *active = hs->active_indices + hs->active_offset[lev];
        for (int a = 0; a < n_active; ++a) {
            int idx_l = active[a];
            int kx = idx_l % ndof_x_l;
            int ky = idx_l / ndof_x_l;
            /* the coefficient vector for this basis function */
            int k_global = hs->active_offset[lev] + a;
            const double *coef = coeffs + k_global*ncomp;

            /* contribute to every fine basis (Kx, Ky):
               coeffs_finest[Ky*ndof_x_fin+Kx] += Mx[Kx,kx] * My[Ky,ky] * coef */
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

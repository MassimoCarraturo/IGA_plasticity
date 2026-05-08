/*
 * qi_localls.c — Parallel implementation of getcoeff_localLS_Bspl.
 *
 * Algorithm (per active hierarchical THB-spline T^lev_kl):
 *   1. Locate the parametric support box [breaks[a_d], breaks[b_d+1]] in each
 *      direction from the cell box [a_d, b_d] = [k_d - p_d, k_d] (clipped).
 *   2. Gather data points whose coordinates lie in that box.  If too few
 *      (count < prod((b-a+1)+deg)), enlarge the box by 1 cell on each side
 *      until sufficient.
 *   3. Determine the set of TP basis functions nonzero on the (expanded) box;
 *      this gives the local row/col index set ind_active_on_supp.
 *   4. Form the local penalised normal-equations system
 *           S = C_loc * diag(w) * C_loc^T + lambda * M_loc
 *      where C_loc and M_loc are submatrices of the level's full collocation
 *      and gradgrad-mass matrices, indexed by ind_active_on_supp (rows /
 *      cols of M_loc; rows of C_loc; cols of C_loc are local data points).
 *   5. Build the right-hand side rhs = C_loc * diag(w) * f_loc  (n_active x ncomp)
 *      and solve S * X = rhs (LU with partial pivoting).
 *   6. Output coefficient row = X[loc_kl, :] where loc_kl is the position of
 *      kl in ind_active_on_supp.
 *
 * The outer loop over k is parallelised with OpenMP; each thread holds its
 * own scratch.  Matrices are stored row-major to match qi_linalg.[ch].
 */

#include "../include/qi.h"
#include "../include/qi_localls.h"
#include "qi_linalg.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef _OPENMP
#include <omp.h>
#endif

/* ---------- Options ------------------------------------------------------ */
void qi_ls_options_defaults(qi_ls_options_t *opts)
{
    if (!opts) return;
    opts->lambda    = 0.0;
    opts->n_threads = 0;
    opts->verbose   = 0;
}

/* ---------- Helpers ------------------------------------------------------ */
static inline int locate_level(const int *act_off, int nlev, int k, int *kloc)
{
    int lo = 0, hi = nlev - 1;
    while (lo < hi) {
        int mid = (lo + hi + 1) >> 1;
        if (act_off[mid] <= k) lo = mid;
        else                    hi = mid - 1;
    }
    *kloc = k - act_off[lo];
    return lo;
}

/* Column-major linear index -> subscripts:
 * idx = kx[0] + kx[1]*nd[0] + kx[2]*nd[0]*nd[1] + ... */
static inline void ind2sub_cm(int idx, int dim, const int *nd, int *kx)
{
    int rem = idx;
    for (int d = 0; d < dim; ++d) {
        kx[d] = rem % nd[d];
        rem  /= nd[d];
    }
}

/* ---------- Per-thread scratch with growable buffers --------------------- */
typedef struct {
    /* TP-index -> local position (used to restrict sparse rows/cols).
     * Allocated to ndof_lev_max once; entries -1 when not in current active set. */
    int    *pos_in_active;     int pos_size;

    int    *active;            int active_cap;       /* current active TP indices */
    int    *ind_loc;           int ind_loc_cap;      /* indices of local data points */
    double *C_loc;             int C_cap;            /* n_active x n_loc, row-major */
    double *Cw;                int Cw_cap;           /* n_active x n_loc, row-major */
    double *S;                 int S_cap;            /* n_active x n_active           */
    double *rhs;               int rhs_cap;          /* n_active x ncomp              */
    int    *piv;               int piv_cap;          /* [n_active]                    */
} scratch_t;

static void scratch_init(scratch_t *s) { memset(s, 0, sizeof(*s)); }

static void scratch_free(scratch_t *s)
{
    free(s->pos_in_active);
    free(s->active);
    free(s->ind_loc);
    free(s->C_loc);
    free(s->Cw);
    free(s->S);
    free(s->rhs);
    free(s->piv);
    memset(s, 0, sizeof(*s));
}

static int ensure_int(int **buf, int *cap, int needed)
{
    if (needed <= 0) return 0;
    if (*cap >= needed) return 0;
    int n = (*cap > 0) ? *cap : 16;
    while (n < needed) n *= 2;
    int *p = (int*)realloc(*buf, (size_t)n * sizeof(int));
    if (!p) return -1;
    *buf = p; *cap = n;
    return 0;
}

static int ensure_double(double **buf, int *cap, int needed)
{
    if (needed <= 0) return 0;
    if (*cap >= needed) return 0;
    int n = (*cap > 0) ? *cap : 16;
    while (n < needed) n *= 2;
    double *p = (double*)realloc(*buf, (size_t)n * sizeof(double));
    if (!p) return -1;
    *buf = p; *cap = n;
    return 0;
}

static int ensure_pos(scratch_t *s, int needed)
{
    if (s->pos_size >= needed) return 0;
    free(s->pos_in_active);
    s->pos_in_active = (int*)malloc((size_t)needed * sizeof(int));
    if (!s->pos_in_active) return -1;
    for (int i = 0; i < needed; ++i) s->pos_in_active[i] = -1;
    s->pos_size = needed;
    return 0;
}

/* ---------- Per-DOF solve ------------------------------------------------ */
static qi_status_t solve_one(int                          k_global,
                             int                          lev,
                             int                          kl,
                             const qi_ls_hspace_t        *hs,
                             const qi_ls_data_t          *dat,
                             const qi_ls_level_precomp_t *plev,
                             double                       lambda,
                             scratch_t                   *S,
                             double                      *coeffs_row)
{
    const int dim = hs->par_dim;
    const int *deg     = hs->degree;
    const int *ndof_d  = hs->ndof_dir + lev * dim;
    const int *nel_d   = hs->nel_dir  + lev * dim;
    const int  ncomp   = dat->ncomp;

    /* Subscripts of kl in level's TP space */
    int kx[3] = {0, 0, 0};
    ind2sub_cm(kl, dim, ndof_d, kx);

    /* Initial support cell box (0-based inclusive).
     * Use the per-cell connectivity (cell c's 1D basis are conn_first[c]..conn_first[c]+p)
     * so this works for any regularity (not just the default p-1).
     * Basis j is supported on cells c where conn_first[c] <= j <= conn_first[c] + p. */
    int a[3] = {0, 0, 0}, b[3] = {0, 0, 0};
    for (int d = 0; d < dim; ++d) {
        int p     = deg[d];
        int j     = kx[d];
        int nel_l = nel_d[d];
        int off   = hs->conn_first_offset[lev * dim + d];
        const int *cf = hs->conn_first_data + off;
        int a_d = -1, b_d = -1;
        for (int c = 0; c < nel_l; ++c) {
            if (cf[c] + p >= j && a_d < 0) a_d = c;
            if (cf[c] <= j)                b_d = c;
        }
        if (a_d < 0) a_d = 0;
        if (b_d < 0) b_d = nel_l - 1;
        a[d] = a_d;
        b[d] = b_d;
    }

    /* Adaptive expansion until enough local points */
    if (ensure_int(&S->ind_loc, &S->ind_loc_cap, dat->npoints) != 0) return QI_ERR_OUT_OF_MEMORY;
    int n_loc = 0, min_req = 0;
    int max_iters = 1;
    for (int d = 0; d < dim; ++d) max_iters += nel_d[d];

    for (int it = 0; it < max_iters; ++it) {
        double pmin[3], pmax[3];
        for (int d = 0; d < dim; ++d) {
            int off = hs->breaks_offset[lev * dim + d];
            pmin[d] = hs->breaks_data[off + a[d]];
            pmax[d] = hs->breaks_data[off + b[d] + 1];
        }
        int prod = 1;
        for (int d = 0; d < dim; ++d) prod *= (b[d] - a[d] + 1) + deg[d];
        min_req = prod;

        n_loc = 0;
        for (int p = 0; p < dat->npoints; ++p) {
            int inside = 1;
            for (int d = 0; d < dim; ++d) {
                double v = dat->coords[(size_t)p * dim + d];
                if (v < pmin[d] || v > pmax[d]) { inside = 0; break; }
            }
            if (inside) S->ind_loc[n_loc++] = p;
        }
        if (n_loc >= min_req || n_loc >= dat->npoints) break;

        int changed = 0;
        for (int d = 0; d < dim; ++d) {
            if (a[d] > 0)            { a[d]--; changed = 1; }
            if (b[d] < nel_d[d] - 1) { b[d]++; changed = 1; }
        }
        if (!changed) break;
    }
    if (n_loc < min_req && n_loc < dat->npoints) {
        return QI_ERR_NOT_ENOUGH_DATA;
    }
    if (n_loc == 0) return QI_ERR_NOT_ENOUGH_DATA;

    /* Active basis on cell box [a, b]:
     * 1D range = [conn_first[a[d]], conn_first[b[d]] + deg[d]] (clipped). */
    int rng_lo[3] = {0,0,0}, rng_hi[3] = {0,0,0}, n_active = 1;
    for (int d = 0; d < dim; ++d) {
        int off = hs->conn_first_offset[lev * dim + d];
        const int *cf = hs->conn_first_data + off;
        rng_lo[d] = cf[a[d]];
        rng_hi[d] = cf[b[d]] + deg[d];
        if (rng_hi[d] > ndof_d[d] - 1) rng_hi[d] = ndof_d[d] - 1;
        if (rng_lo[d] > rng_hi[d])     return QI_ERR_INVALID_ARG;
        n_active *= (rng_hi[d] - rng_lo[d] + 1);
    }
    if (ensure_int(&S->active, &S->active_cap, n_active) != 0) return QI_ERR_OUT_OF_MEMORY;
    {
        int idx = 0;
        if (dim == 2) {
            for (int j = rng_lo[1]; j <= rng_hi[1]; ++j)
                for (int i = rng_lo[0]; i <= rng_hi[0]; ++i)
                    S->active[idx++] = i + j * ndof_d[0];
        } else {  /* dim == 3 */
            int nx = ndof_d[0], nxy = ndof_d[0] * ndof_d[1];
            for (int kk = rng_lo[2]; kk <= rng_hi[2]; ++kk)
                for (int j = rng_lo[1]; j <= rng_hi[1]; ++j)
                    for (int i = rng_lo[0]; i <= rng_hi[0]; ++i)
                        S->active[idx++] = i + j * nx + kk * nxy;
        }
    }

    /* Position of kl in active */
    int loc_kl = -1;
    for (int i = 0; i < n_active; ++i) {
        if (S->active[i] == kl) { loc_kl = i; break; }
    }
    if (loc_kl < 0) return QI_ERR_INVALID_ARG;

    /* Set pos_in_active map for this DOF (row-restricted sparse extraction) */
    for (int i = 0; i < n_active; ++i) S->pos_in_active[S->active[i]] = i;

    /* Working buffers */
    if (ensure_double(&S->C_loc, &S->C_cap,   n_active * n_loc)         != 0) goto memfail;
    if (ensure_double(&S->Cw,    &S->Cw_cap,  n_active * n_loc)         != 0) goto memfail;
    if (ensure_double(&S->S,     &S->S_cap,   n_active * n_active)      != 0) goto memfail;
    if (ensure_double(&S->rhs,   &S->rhs_cap, n_active * ncomp)         != 0) goto memfail;
    if (ensure_int   (&S->piv,   &S->piv_cap, n_active)                 != 0) goto memfail;

    /* Extract C_loc (row-major, n_active rows x n_loc cols).
     * C_loc[r, c] = col_full[active[r], ind_loc[c]]. */
    memset(S->C_loc, 0, (size_t)n_active * n_loc * sizeof(double));
    {
        const qi_csc_t *cm = &plev->col;
        for (int c = 0; c < n_loc; ++c) {
            int j = S->ind_loc[c];
            int p0 = cm->col_ptr[j];
            int p1 = cm->col_ptr[j + 1];
            for (int q = p0; q < p1; ++q) {
                int row = cm->row_ind[q];
                int rr  = (row >= 0 && row < S->pos_size) ? S->pos_in_active[row] : -1;
                if (rr >= 0) {
                    S->C_loc[(size_t)rr * n_loc + c] = cm->values[q];
                }
            }
        }
    }

    /* Cw[r, c] = C_loc[r, c] * w[ind_loc[c]] */
    {
        for (int c = 0; c < n_loc; ++c) {
            double wc = (dat->weight) ? dat->weight[S->ind_loc[c]] : 1.0;
            for (int r = 0; r < n_active; ++r) {
                S->Cw[(size_t)r * n_loc + c] = S->C_loc[(size_t)r * n_loc + c] * wc;
            }
        }
    }

    /* S = Cw * C_loc^T  (n_active x n_active), row-major.
     * S[r, s] = sum_c Cw[r, c] * C_loc[s, c] */
    {
        memset(S->S, 0, (size_t)n_active * n_active * sizeof(double));
        for (int r = 0; r < n_active; ++r) {
            const double *cwr = S->Cw    + (size_t)r * n_loc;
            for (int s = 0; s < n_active; ++s) {
                const double *cls = S->C_loc + (size_t)s * n_loc;
                double acc = 0.0;
                for (int c = 0; c < n_loc; ++c) acc += cwr[c] * cls[c];
                S->S[(size_t)r * n_active + s] = acc;
            }
        }
    }

    /* Add lambda * M_loc */
    if (lambda != 0.0) {
        const qi_csc_t *mm = &plev->mass;
        /* M_loc[r, s] = mass_full[active[r], active[s]] */
        for (int s = 0; s < n_active; ++s) {
            int j  = S->active[s];
            int p0 = mm->col_ptr[j];
            int p1 = mm->col_ptr[j + 1];
            for (int q = p0; q < p1; ++q) {
                int row = mm->row_ind[q];
                int rr  = (row >= 0 && row < S->pos_size) ? S->pos_in_active[row] : -1;
                if (rr >= 0) {
                    S->S[(size_t)rr * n_active + s] += lambda * mm->values[q];
                }
            }
        }
    }

    /* rhs[r, k] = sum_c Cw[r, c] * f[ind_loc[c]*ncomp + k]   (row-major: rhs[r*ncomp + k]) */
    {
        memset(S->rhs, 0, (size_t)n_active * ncomp * sizeof(double));
        for (int c = 0; c < n_loc; ++c) {
            int p = S->ind_loc[c];
            const double *fp = dat->f + (size_t)p * ncomp;
            for (int r = 0; r < n_active; ++r) {
                double cw = S->Cw[(size_t)r * n_loc + c];
                if (cw == 0.0) continue;
                double *rs = S->rhs + (size_t)r * ncomp;
                for (int kk = 0; kk < ncomp; ++kk) rs[kk] += cw * fp[kk];
            }
        }
    }

    /* Solve S * X = rhs in place (rhs is overwritten with X) */
    {
        int rc = qi_la_lu_solve(n_active, ncomp, S->S, S->rhs, S->piv);
        if (rc != 0) {
            /* Singular: fall back to pinv (A^T A eigendecomposition).  We need
             * scratch ~ n_active * (n_active + 1 + 4 + ncomp) doubles + n_active
             * ints; reuse C_loc / Cw arrays which have been freed of meaning. */
            int lwork = n_active * n_active + n_active + 4 * n_active + n_active * ncomp + 16;
            if (ensure_double(&S->C_loc, &S->C_cap, lwork) != 0) goto memfail;
            /* A copy of S is needed because qi_la_pinv_solve overwrites it; but
             * S has just been destroyed by the failed LU.  Recompute. */
            /* Quick-and-dirty: re-form S, then solve. */
            memset(S->S, 0, (size_t)n_active * n_active * sizeof(double));
            for (int r = 0; r < n_active; ++r) {
                const double *cwr = S->Cw    + (size_t)r * n_loc;
                for (int s = 0; s < n_active; ++s) {
                    const double *cls = S->C_loc + (size_t)s * n_loc;
                    double acc = 0.0;
                    for (int c = 0; c < n_loc; ++c) acc += cwr[c] * cls[c];
                    S->S[(size_t)r * n_active + s] = acc;
                }
            }
            if (lambda != 0.0) {
                const qi_csc_t *mm = &plev->mass;
                for (int s = 0; s < n_active; ++s) {
                    int j  = S->active[s];
                    int p0 = mm->col_ptr[j];
                    int p1 = mm->col_ptr[j + 1];
                    for (int q = p0; q < p1; ++q) {
                        int row = mm->row_ind[q];
                        int rr  = (row >= 0 && row < S->pos_size) ? S->pos_in_active[row] : -1;
                        if (rr >= 0)
                            S->S[(size_t)rr * n_active + s] += lambda * mm->values[q];
                    }
                }
            }
            /* rhs already correct (qi_la_lu_solve only partly overwrote it on failure;
             * but to be safe, recompute rhs into the freshly-allocated C_loc buffer). */
            double *rhs_tmp = S->C_loc;  /* reuse */
            memset(rhs_tmp, 0, (size_t)n_active * ncomp * sizeof(double));
            for (int c = 0; c < n_loc; ++c) {
                int p = S->ind_loc[c];
                const double *fp = dat->f + (size_t)p * ncomp;
                for (int r = 0; r < n_active; ++r) {
                    double cw = S->Cw[(size_t)r * n_loc + c];
                    if (cw == 0.0) continue;
                    double *rs = rhs_tmp + (size_t)r * ncomp;
                    for (int kk = 0; kk < ncomp; ++kk) rs[kk] += cw * fp[kk];
                }
            }
            double sigma_min = 0.0;
            int rc2 = qi_la_pinv_solve(n_active, n_active, ncomp,
                                       S->S, rhs_tmp, S->rhs, &sigma_min,
                                       S->Cw, S->Cw_cap);
            if (rc2 != 0) {
                /* Reset pos_in_active before bailing out */
                for (int i = 0; i < n_active; ++i) S->pos_in_active[S->active[i]] = -1;
                return QI_ERR_SINGULAR_MATRIX;
            }
        }
    }

    /* Output: row loc_kl of the solution */
    {
        const double *row = S->rhs + (size_t)loc_kl * ncomp;
        for (int kk = 0; kk < ncomp; ++kk) coeffs_row[kk] = row[kk];
    }

    /* Reset pos_in_active for the visited entries */
    for (int i = 0; i < n_active; ++i) S->pos_in_active[S->active[i]] = -1;
    return QI_OK;

memfail:
    for (int i = 0; i < n_active; ++i) S->pos_in_active[S->active[i]] = -1;
    return QI_ERR_OUT_OF_MEMORY;
}

/* ---------- Driver ------------------------------------------------------- */
qi_status_t qi_localLS_compute(const qi_ls_hspace_t        *hs,
                               const qi_ls_data_t          *dat,
                               const qi_ls_level_precomp_t *level_pre,
                               const qi_ls_options_t       *opts,
                               double                      *coeffs)
{
    if (!hs || !dat || !level_pre || !coeffs) return QI_ERR_INVALID_ARG;
    if (hs->par_dim != 2 && hs->par_dim != 3)  return QI_ERR_INVALID_ARG;
    if (dat->par_dim != hs->par_dim)           return QI_ERR_INVALID_ARG;

    qi_ls_options_t local_opts;
    qi_ls_options_defaults(&local_opts);
    if (opts) local_opts = *opts;

#ifdef _OPENMP
    if (local_opts.n_threads > 0) omp_set_num_threads(local_opts.n_threads);
#endif

    /* Maximum ndof per level (used to pre-size each thread's pos_in_active map) */
    int max_ndof_lev = 1;
    for (int lv = 0; lv < hs->nlevels; ++lv) {
        int n = 1;
        for (int d = 0; d < hs->par_dim; ++d) n *= hs->ndof_dir[lv * hs->par_dim + d];
        if (n > max_ndof_lev) max_ndof_lev = n;
    }

    /* Pre-compute level offsets for binary search in locate_level */
    int err_count = 0;

    #pragma omp parallel reduction(+:err_count)
    {
        scratch_t S;
        scratch_init(&S);
        if (ensure_pos(&S, max_ndof_lev) != 0) {
            err_count += 1;
        } else {
            int k;
            #pragma omp for schedule(dynamic, 16)
            for (k = 0; k < hs->ndof; ++k) {
                int kloc = 0;
                int lev  = locate_level(hs->active_offset, hs->nlevels, k, &kloc);
                int kl   = hs->active_indices[hs->active_offset[lev] + kloc];

                qi_status_t st = solve_one(k, lev, kl, hs, dat,
                                           &level_pre[lev], local_opts.lambda,
                                           &S, coeffs + (size_t)k * dat->ncomp);
                if (st != QI_OK) {
                    err_count += 1;
                    /* Fill with zeros on failure */
                    for (int c = 0; c < dat->ncomp; ++c) coeffs[(size_t)k * dat->ncomp + c] = 0.0;
                }
            }
        }
        scratch_free(&S);
    }

    if (err_count > 0) {
        if (local_opts.verbose) {
            fprintf(stderr, "qi_localLS_compute: %d basis functions failed\n", err_count);
        }
    }
    return QI_OK;
}

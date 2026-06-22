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
 *
 * v1.1 optimisations:
 *   - Spatial grid index for O(N_local) data point queries (was O(N_total)).
 *   - Symmetric S formation: only upper triangle computed then mirrored.
 *   - Save/restore S and rhs before LU to avoid recomputation on fallback.
 *   - restrict qualifiers throughout for auto-vectorisation.
 */

#include "../include/qi.h"
#include "../include/qi_localls.h"
#include "qi_linalg.h"
#include "qi_spatial.h"

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
    double *S_bak;             int S_bak_cap;        /* backup of S before LU         */
    double *rhs;               int rhs_cap;          /* n_active x ncomp              */
    double *rhs_bak;           int rhs_bak_cap;      /* backup of rhs before LU       */
    int    *piv;               int piv_cap;          /* [n_active]                    */
    /* Diagnostics */
    int     n_active_last;                           /* local system size of last solve */
    int     pinv_used_last;                          /* 1 if last solve used pinv fallback */
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
    free(s->S_bak);
    free(s->rhs);
    free(s->rhs_bak);
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
                             const qi_spatial_t          *spatial,
                             const double *const         *soa_coords,
                             double                       lambda,
                             scratch_t                   *SC,
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
    if (ensure_int(&SC->ind_loc, &SC->ind_loc_cap, dat->npoints) != 0) return QI_ERR_OUT_OF_MEMORY;
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

        /* --- Optimisation #1: spatial grid query instead of brute-force --- */
        n_loc = qi_spatial_query_box_nd(spatial, pmin, pmax, soa_coords,
                                        SC->ind_loc, SC->ind_loc_cap);
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
    if (ensure_int(&SC->active, &SC->active_cap, n_active) != 0) return QI_ERR_OUT_OF_MEMORY;
    {
        int idx = 0;
        if (dim == 2) {
            for (int j = rng_lo[1]; j <= rng_hi[1]; ++j)
                for (int i = rng_lo[0]; i <= rng_hi[0]; ++i)
                    SC->active[idx++] = i + j * ndof_d[0];
        } else {  /* dim == 3 */
            int nx = ndof_d[0], nxy = ndof_d[0] * ndof_d[1];
            for (int kk = rng_lo[2]; kk <= rng_hi[2]; ++kk)
                for (int j = rng_lo[1]; j <= rng_hi[1]; ++j)
                    for (int i = rng_lo[0]; i <= rng_hi[0]; ++i)
                        SC->active[idx++] = i + j * nx + kk * nxy;
        }
    }

    /* Position of kl in active */
    int loc_kl = -1;
    for (int i = 0; i < n_active; ++i) {
        if (SC->active[i] == kl) { loc_kl = i; break; }
    }
    if (loc_kl < 0) return QI_ERR_INVALID_ARG;

    /* Set pos_in_active map for this DOF (row-restricted sparse extraction) */
    for (int i = 0; i < n_active; ++i) SC->pos_in_active[SC->active[i]] = i;

    /* Working buffers */
    if (ensure_double(&SC->C_loc,   &SC->C_cap,       n_active * n_loc)    != 0) goto memfail;
    if (ensure_double(&SC->Cw,      &SC->Cw_cap,      n_active * n_loc)    != 0) goto memfail;
    if (ensure_double(&SC->S,       &SC->S_cap,       n_active * n_active) != 0) goto memfail;
    if (ensure_double(&SC->S_bak,   &SC->S_bak_cap,   n_active * n_active) != 0) goto memfail;
    if (ensure_double(&SC->rhs,     &SC->rhs_cap,     n_active * ncomp)    != 0) goto memfail;
    if (ensure_double(&SC->rhs_bak, &SC->rhs_bak_cap, n_active * ncomp)    != 0) goto memfail;
    if (ensure_int   (&SC->piv,     &SC->piv_cap,     n_active)            != 0) goto memfail;

    /* Extract C_loc (row-major, n_active rows x n_loc cols).
     * C_loc[r, c] = col_full[active[r], ind_loc[c]]. */
    memset(SC->C_loc, 0, (size_t)n_active * n_loc * sizeof(double));
    {
        const qi_csc_t *cm = &plev->col;
        for (int c = 0; c < n_loc; ++c) {
            int j = SC->ind_loc[c];
            int p0 = cm->col_ptr[j];
            int p1 = cm->col_ptr[j + 1];
            for (int q = p0; q < p1; ++q) {
                int row = cm->row_ind[q];
                int rr  = (row >= 0 && row < SC->pos_size) ? SC->pos_in_active[row] : -1;
                if (rr >= 0) {
                    SC->C_loc[(size_t)rr * n_loc + c] = cm->values[q];
                }
            }
        }
    }

    /* Cw[r, c] = C_loc[r, c] * w[ind_loc[c]] */
    {
        for (int c = 0; c < n_loc; ++c) {
            double wc = (dat->weight) ? dat->weight[SC->ind_loc[c]] : 1.0;
            for (int r = 0; r < n_active; ++r) {
                SC->Cw[(size_t)r * n_loc + c] = SC->C_loc[(size_t)r * n_loc + c] * wc;
            }
        }
    }

    /* --- Optimisation #2: Symmetric S formation ---
     * S = Cw * C_loc^T is symmetric (since S = C*diag(w)*C^T with w >= 0).
     * Only compute upper triangle (r <= s) and mirror. */
    {
        memset(SC->S, 0, (size_t)n_active * n_active * sizeof(double));
        for (int r = 0; r < n_active; ++r) {
            const double *QI_RESTRICT cwr = SC->Cw    + (size_t)r * n_loc;
            for (int s = r; s < n_active; ++s) {
                const double *QI_RESTRICT cls = SC->C_loc + (size_t)s * n_loc;
                double acc = 0.0;
                for (int c = 0; c < n_loc; ++c) acc += cwr[c] * cls[c];
                SC->S[(size_t)r * n_active + s] = acc;
                if (s != r)
                    SC->S[(size_t)s * n_active + r] = acc;
            }
        }
    }

    /* Add lambda * M_loc */
    if (lambda != 0.0) {
        const qi_csc_t *mm = &plev->mass;
        /* M_loc[r, s] = mass_full[active[r], active[s]] */
        for (int s = 0; s < n_active; ++s) {
            int j  = SC->active[s];
            int p0 = mm->col_ptr[j];
            int p1 = mm->col_ptr[j + 1];
            for (int q = p0; q < p1; ++q) {
                int row = mm->row_ind[q];
                int rr  = (row >= 0 && row < SC->pos_size) ? SC->pos_in_active[row] : -1;
                if (rr >= 0) {
                    SC->S[(size_t)rr * n_active + s] += lambda * mm->values[q];
                }
            }
        }
    }

    /* rhs[r, k] = sum_c Cw[r, c] * f[ind_loc[c]*ncomp + k]   (row-major: rhs[r*ncomp + k]) */
    {
        memset(SC->rhs, 0, (size_t)n_active * ncomp * sizeof(double));
        for (int c = 0; c < n_loc; ++c) {
            int p = SC->ind_loc[c];
            const double *QI_RESTRICT fp = dat->f + (size_t)p * ncomp;
            for (int r = 0; r < n_active; ++r) {
                double cw = SC->Cw[(size_t)r * n_loc + c];
                if (cw == 0.0) continue;
                double *QI_RESTRICT rs = SC->rhs + (size_t)r * ncomp;
                for (int kk = 0; kk < ncomp; ++kk) rs[kk] += cw * fp[kk];
            }
        }
    }

    /* --- Optimisation #6: Save S and rhs before LU to avoid recomputation --- */
    {
        size_t s_bytes   = (size_t)n_active * n_active * sizeof(double);
        size_t rhs_bytes = (size_t)n_active * ncomp    * sizeof(double);
        memcpy(SC->S_bak,   SC->S,   s_bytes);
        memcpy(SC->rhs_bak, SC->rhs, rhs_bytes);
    }

    /* Diagnostics: record system size */
    SC->n_active_last = n_active;
    SC->pinv_used_last = 0;

    /* Solve S * X = rhs in place (rhs is overwritten with X) */
    {
        int rc = qi_la_lu_solve(n_active, ncomp, SC->S, SC->rhs, SC->piv);
        if (rc != 0) {
            /* Singular: fall back to pinv (A^T A eigendecomposition).
             * Restore S and rhs from backup instead of recomputing. */
            SC->pinv_used_last = 1;
            size_t s_bytes   = (size_t)n_active * n_active * sizeof(double);
            size_t rhs_bytes = (size_t)n_active * ncomp    * sizeof(double);
            memcpy(SC->S,   SC->S_bak,   s_bytes);
            memcpy(SC->rhs, SC->rhs_bak, rhs_bytes);

            int lwork = 2*n_active*n_active + n_active + 2*n_active*ncomp;
            if (ensure_double(&SC->C_loc, &SC->C_cap, lwork) != 0) goto memfail;
            double sigma_min = 0.0;
            int rc2 = qi_la_pinv_solve(n_active, n_active, ncomp,
                                       SC->S, SC->rhs, SC->rhs_bak, &sigma_min,
                                       SC->C_loc, SC->C_cap);
            if (rc2 != 0) {
                /* Reset pos_in_active before bailing out */
                for (int i = 0; i < n_active; ++i) SC->pos_in_active[SC->active[i]] = -1;
                return QI_ERR_SINGULAR_MATRIX;
            }
            /* pinv_solve wrote into rhs_bak as output x; copy result to rhs */
            memcpy(SC->rhs, SC->rhs_bak, rhs_bytes);
        }
    }

    /* Output: row loc_kl of the solution */
    {
        const double *row = SC->rhs + (size_t)loc_kl * ncomp;
        for (int kk = 0; kk < ncomp; ++kk) coeffs_row[kk] = row[kk];
    }

    /* Reset pos_in_active for the visited entries */
    for (int i = 0; i < n_active; ++i) SC->pos_in_active[SC->active[i]] = -1;
    return QI_OK;

memfail:
    for (int i = 0; i < n_active; ++i) SC->pos_in_active[SC->active[i]] = -1;
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

    const int dim = hs->par_dim;

    /* Maximum ndof per level (used to pre-size each thread's pos_in_active map) */
    int max_ndof_lev = 1;
    for (int lv = 0; lv < hs->nlevels; ++lv) {
        int n = 1;
        for (int d = 0; d < hs->par_dim; ++d) n *= hs->ndof_dir[lv * hs->par_dim + d];
        if (n > max_ndof_lev) max_ndof_lev = n;
    }

    /* --- Optimisation #1: build spatial grid index from data points --- */
    /* Convert AoS coords (dat->coords[p*dim+d]) to SoA for qi_spatial */
    double *soa_buf = (double*)malloc((size_t)dat->npoints * dim * sizeof(double));
    if (!soa_buf) return QI_ERR_OUT_OF_MEMORY;
    const double **soa_ptrs = (const double**)malloc((size_t)dim * sizeof(double*));
    if (!soa_ptrs) { free(soa_buf); return QI_ERR_OUT_OF_MEMORY; }
    for (int d = 0; d < dim; ++d) {
        double *col = soa_buf + (size_t)d * dat->npoints;
        for (int p = 0; p < dat->npoints; ++p)
            col[p] = dat->coords[(size_t)p * dim + d];
        soa_ptrs[d] = col;
    }

    qi_spatial_t spatial;
    {
        /* Pick bin count from the finest-level element count */
        int finest_nel = 1;
        for (int lv = 0; lv < hs->nlevels; ++lv) {
            for (int d = 0; d < dim; ++d) {
                if (hs->nel_dir[lv * dim + d] > finest_nel)
                    finest_nel = hs->nel_dir[lv * dim + d];
            }
        }
        qi_spatial_build_nd(&spatial, dim, dat->npoints, soa_ptrs, finest_nel);
    }

    /* Pre-compute level offsets for binary search in locate_level */
    int err_count = 0;
    int pinv_count = 0;       /* diagnostic: DOFs that fell back to pinv */
    int max_n_active = 0;     /* diagnostic: largest local system size */

    #pragma omp parallel reduction(+:err_count) reduction(+:pinv_count)
    {
        scratch_t SC;
        scratch_init(&SC);
        if (ensure_pos(&SC, max_ndof_lev) != 0) {
            err_count += 1;
        } else {
            int k;
            #pragma omp for schedule(dynamic, 16)
            for (k = 0; k < hs->ndof; ++k) {
                int kloc = 0;
                int lev  = locate_level(hs->active_offset, hs->nlevels, k, &kloc);
                int kl   = hs->active_indices[hs->active_offset[lev] + kloc];

                /* Track local system size before solve */
                int n_active_before = SC.n_active_last;

                qi_status_t st = solve_one(k, lev, kl, hs, dat,
                                           &level_pre[lev],
                                           &spatial, soa_ptrs,
                                           local_opts.lambda,
                                           &SC, coeffs + (size_t)k * dat->ncomp);
                if (SC.n_active_last > max_n_active) {
                    #pragma omp critical
                    if (SC.n_active_last > max_n_active)
                        max_n_active = SC.n_active_last;
                }
                if (SC.pinv_used_last)
                    pinv_count += 1;
                if (st != QI_OK) {
                    err_count += 1;
                    /* Fill with zeros on failure */
                    for (int c = 0; c < dat->ncomp; ++c) coeffs[(size_t)k * dat->ncomp + c] = 0.0;
                }
            }
        }
        scratch_free(&SC);
    }

    /* Cleanup spatial index and SoA buffers */
    qi_spatial_free(&spatial);
    free((void*)soa_ptrs);
    free(soa_buf);

    /* Diagnostic output: always print summary when pinv fallback used */
    if (pinv_count > 0 || err_count > 0) {
        fprintf(stderr, "qi_localLS_compute: ndof=%d, lambda=%.2e, "
                "pinv_fallback=%d/%d (%.1f%%), errors=%d, max_sys_size=%d\n",
                hs->ndof, local_opts.lambda,
                pinv_count, hs->ndof,
                100.0 * pinv_count / (hs->ndof > 0 ? hs->ndof : 1),
                err_count, max_n_active);
    } else {
        fprintf(stderr, "qi_localLS_compute: ndof=%d, lambda=%.2e, "
                "all LU ok, max_sys_size=%d\n",
                hs->ndof, local_opts.lambda, max_n_active);
    }
    return QI_OK;
}

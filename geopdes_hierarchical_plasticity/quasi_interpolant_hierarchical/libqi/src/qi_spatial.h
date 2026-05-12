/*
 * qi_spatial.h — Uniform-grid spatial index for fast point-in-box/disk queries.
 *
 * Thread-safe once built: all query functions are read-only.
 * Built via counting sort in O(N) time; each query touches only the bins
 * overlapping the query region, giving O(N_local) cost instead of O(N_total).
 */
#ifndef QI_SPATIAL_H
#define QI_SPATIAL_H

#define QI_SPATIAL_MAX_DIM 3

typedef struct {
    int    ndim;
    int    nbins[QI_SPATIAL_MAX_DIM];
    double lo[QI_SPATIAL_MAX_DIM];
    double hi[QI_SPATIAL_MAX_DIM];
    double inv_bin[QI_SPATIAL_MAX_DIM];   /* nbins[d] / (hi[d] - lo[d]) */
    int    total_bins;
    int   *bin_start;     /* [total_bins + 1]  CSR offsets  */
    int   *bin_sorted;    /* [npoints]  indices sorted by bin */
    int    npoints;
} qi_spatial_t;

/* Build a 2D index from separate x,y arrays (SoA).
 * nbins_hint: bins per dimension (0 = auto ~sqrt(npoints)). */
void qi_spatial_build_2d(qi_spatial_t *idx, int npoints,
                         const double *x, const double *y,
                         int nbins_hint);

/* Build an n-dimensional index from SoA coordinate arrays. */
void qi_spatial_build_nd(qi_spatial_t *idx, int ndim, int npoints,
                         const double *const *coords,
                         int nbins_hint);

/* Free index memory. Safe to call on a zeroed struct. */
void qi_spatial_free(qi_spatial_t *idx);

/* Box query (n-dim): return indices of points inside [lo_q, hi_q].
 * coords must be the same SoA arrays used to build the index.
 * Returns count of points found; writes at most result_cap indices. */
int qi_spatial_query_box_nd(const qi_spatial_t *idx,
                            const double *lo_q, const double *hi_q,
                            const double *const *coords,
                            int *result, int result_cap);

/* Disk query (2D only): return indices of points within radius r of (cx,cy).
 * x,y must be the same arrays used to build the index. */
int qi_spatial_query_disk_2d(const qi_spatial_t *idx,
                             double cx, double cy, double r,
                             const double *x, const double *y,
                             int *result, int result_cap);

#endif /* QI_SPATIAL_H */

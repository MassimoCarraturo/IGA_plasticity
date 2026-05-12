/*
 * qi_spatial.c — Uniform-grid spatial index (counting-sort based).
 *
 * Build: O(N) two-pass counting sort.
 * Query: O(N_local) — only bins overlapping the query region are visited.
 */
#include "qi_spatial.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* ------------------------------------------------------------------ */
/* Internal builder                                                    */
/* ------------------------------------------------------------------ */
static void build_internal(qi_spatial_t *idx, int ndim, int npoints,
                           const double *const *coords, int nbins_hint)
{
    memset(idx, 0, sizeof(*idx));
    idx->ndim    = ndim;
    idx->npoints = npoints;

    if (npoints <= 0) {
        idx->total_bins = 1;
        idx->bin_start  = (int *)calloc(2, sizeof(int));
        idx->bin_sorted = NULL;
        return;
    }

    /* ---- bounding box ---- */
    for (int d = 0; d < ndim; ++d) {
        idx->lo[d] = coords[d][0];
        idx->hi[d] = coords[d][0];
    }
    for (int i = 1; i < npoints; ++i)
        for (int d = 0; d < ndim; ++d) {
            double v = coords[d][i];
            if (v < idx->lo[d]) idx->lo[d] = v;
            if (v > idx->hi[d]) idx->hi[d] = v;
        }

    /* ---- bin count ---- */
    int nb = (nbins_hint > 0) ? nbins_hint
                              : (int)ceil(sqrt((double)npoints));
    if (nb < 1)    nb = 1;
    if (nb > 1024) nb = 1024;

    idx->total_bins = 1;
    for (int d = 0; d < ndim; ++d) {
        idx->nbins[d] = nb;
        double range = idx->hi[d] - idx->lo[d];
        if (range < 1e-15) range = 1.0;
        idx->hi[d]     = idx->lo[d] + range * 1.0000001;   /* epsilon guard */
        idx->inv_bin[d] = (double)nb / (idx->hi[d] - idx->lo[d]);
        idx->total_bins *= nb;
    }

    /* ---- counting sort (two passes) ---- */
    idx->bin_start  = (int *)calloc((size_t)(idx->total_bins + 1), sizeof(int));
    idx->bin_sorted = (int *)malloc((size_t)npoints * sizeof(int));

    /* helper: compute bin for point i */
    #define BIN_OF(i) do {                                      \
        int _bin = 0, _stride = 1;                              \
        for (int _d = 0; _d < ndim; ++_d) {                    \
            int _b = (int)((coords[_d][(i)] - idx->lo[_d])     \
                           * idx->inv_bin[_d]);                 \
            if (_b < 0) _b = 0;                                 \
            if (_b >= idx->nbins[_d]) _b = idx->nbins[_d] - 1; \
            _bin += _b * _stride;                               \
            _stride *= idx->nbins[_d];                          \
        }                                                       \
        bin = _bin;                                             \
    } while (0)

    /* pass 1: count */
    for (int i = 0; i < npoints; ++i) {
        int bin;
        BIN_OF(i);
        idx->bin_start[bin + 1]++;
    }

    /* prefix sum */
    for (int b = 0; b < idx->total_bins; ++b)
        idx->bin_start[b + 1] += idx->bin_start[b];

    /* pass 2: place */
    int *wpos = (int *)malloc((size_t)idx->total_bins * sizeof(int));
    memcpy(wpos, idx->bin_start, (size_t)idx->total_bins * sizeof(int));
    for (int i = 0; i < npoints; ++i) {
        int bin;
        BIN_OF(i);
        idx->bin_sorted[wpos[bin]++] = i;
    }
    free(wpos);
    #undef BIN_OF
}

/* ------------------------------------------------------------------ */
/* Public builders                                                     */
/* ------------------------------------------------------------------ */
void qi_spatial_build_2d(qi_spatial_t *idx, int npoints,
                         const double *x, const double *y,
                         int nbins_hint)
{
    const double *c[2] = { x, y };
    build_internal(idx, 2, npoints, c, nbins_hint);
}

void qi_spatial_build_nd(qi_spatial_t *idx, int ndim, int npoints,
                         const double *const *coords, int nbins_hint)
{
    build_internal(idx, ndim, npoints, coords, nbins_hint);
}

void qi_spatial_free(qi_spatial_t *idx)
{
    free(idx->bin_start);
    free(idx->bin_sorted);
    memset(idx, 0, sizeof(*idx));
}

/* ------------------------------------------------------------------ */
/* Queries                                                             */
/* ------------------------------------------------------------------ */
int qi_spatial_query_box_nd(const qi_spatial_t *idx,
                            const double *lo_q, const double *hi_q,
                            const double *const *coords,
                            int *result, int result_cap)
{
    if (!idx->bin_start || idx->npoints <= 0) return 0;

    const int ndim = idx->ndim;
    int blo[QI_SPATIAL_MAX_DIM], bhi[QI_SPATIAL_MAX_DIM];
    for (int d = 0; d < ndim; ++d) {
        blo[d] = (int)((lo_q[d] - idx->lo[d]) * idx->inv_bin[d]);
        bhi[d] = (int)((hi_q[d] - idx->lo[d]) * idx->inv_bin[d]);
        if (blo[d] < 0) blo[d] = 0;
        if (bhi[d] >= idx->nbins[d]) bhi[d] = idx->nbins[d] - 1;
        if (blo[d] > bhi[d]) return 0;
    }

    int nfound = 0;

    if (ndim == 2) {
        const int nx = idx->nbins[0];
        for (int by = blo[1]; by <= bhi[1]; ++by) {
            for (int bx = blo[0]; bx <= bhi[0]; ++bx) {
                int bin = bx + by * nx;
                int p0 = idx->bin_start[bin];
                int p1 = idx->bin_start[bin + 1];
                for (int q = p0; q < p1; ++q) {
                    int i = idx->bin_sorted[q];
                    if (coords[0][i] >= lo_q[0] && coords[0][i] <= hi_q[0] &&
                        coords[1][i] >= lo_q[1] && coords[1][i] <= hi_q[1]) {
                        if (nfound < result_cap) result[nfound] = i;
                        nfound++;
                    }
                }
            }
        }
    } else {   /* ndim == 3 */
        const int nx = idx->nbins[0];
        const int nxy = idx->nbins[0] * idx->nbins[1];
        for (int bz = blo[2]; bz <= bhi[2]; ++bz) {
            for (int by = blo[1]; by <= bhi[1]; ++by) {
                for (int bx = blo[0]; bx <= bhi[0]; ++bx) {
                    int bin = bx + by * nx + bz * nxy;
                    int p0 = idx->bin_start[bin];
                    int p1 = idx->bin_start[bin + 1];
                    for (int q = p0; q < p1; ++q) {
                        int i = idx->bin_sorted[q];
                        int ok = 1;
                        for (int d = 0; d < 3; ++d)
                            if (coords[d][i] < lo_q[d] ||
                                coords[d][i] > hi_q[d]) { ok = 0; break; }
                        if (ok) {
                            if (nfound < result_cap) result[nfound] = i;
                            nfound++;
                        }
                    }
                }
            }
        }
    }
    return nfound;
}

int qi_spatial_query_disk_2d(const qi_spatial_t *idx,
                             double cx, double cy, double r,
                             const double *x, const double *y,
                             int *result, int result_cap)
{
    if (!idx->bin_start || idx->npoints <= 0) return 0;

    int bxlo = (int)((cx - r - idx->lo[0]) * idx->inv_bin[0]);
    int bxhi = (int)((cx + r - idx->lo[0]) * idx->inv_bin[0]);
    int bylo = (int)((cy - r - idx->lo[1]) * idx->inv_bin[1]);
    int byhi = (int)((cy + r - idx->lo[1]) * idx->inv_bin[1]);
    if (bxlo < 0) bxlo = 0;
    if (bxhi >= idx->nbins[0]) bxhi = idx->nbins[0] - 1;
    if (bylo < 0) bylo = 0;
    if (byhi >= idx->nbins[1]) byhi = idx->nbins[1] - 1;

    const double r2 = r * r;
    int nfound = 0;

    for (int by = bylo; by <= byhi; ++by) {
        for (int bx = bxlo; bx <= bxhi; ++bx) {
            int bin = bx + by * idx->nbins[0];
            int p0 = idx->bin_start[bin];
            int p1 = idx->bin_start[bin + 1];
            for (int q = p0; q < p1; ++q) {
                int i  = idx->bin_sorted[q];
                double dx = x[i] - cx;
                double dy = y[i] - cy;
                if (dx * dx + dy * dy <= r2) {
                    if (nfound < result_cap) result[nfound] = i;
                    nfound++;
                }
            }
        }
    }
    return nfound;
}

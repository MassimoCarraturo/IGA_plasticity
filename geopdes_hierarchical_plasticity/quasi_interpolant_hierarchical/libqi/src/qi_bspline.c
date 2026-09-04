/*
 * qi_bspline.c — univariate B-spline routines (NURBS Book Algs 2.1 & 2.2).
 *
 * v1.1: restrict qualifiers for auto-vectorisation.
 */
#include "qi_bspline.h"
#include <string.h>

int qi_bs_findspan(int n, int p, double u, const double *QI_RESTRICT U)
{
    /* clamp u to last span if at the right boundary */
    if (u >= U[n + 1]) return n;
    if (u <= U[p])     return p;

    int low  = p;
    int high = n + 1;
    int mid  = (low + high) / 2;
    while (u < U[mid] || u >= U[mid + 1]) {
        if (u < U[mid]) high = mid;
        else            low  = mid;
        mid = (low + high) / 2;
    }
    return mid;
}

void qi_bs_basisfun(int i, double u, int p, const double *QI_RESTRICT U,
                    double *QI_RESTRICT N)
{
    /* N[0..p] -- standard NURBS Book Alg 2.2.
       Uses local left[1..p], right[1..p].  Bound p <= 8 in this library. */
    double left [16];
    double right[16];

    N[0] = 1.0;
    for (int j = 1; j <= p; ++j) {
        left [j] = u - U[i + 1 - j];
        right[j] = U[i + j] - u;
        double saved = 0.0;
        for (int r = 0; r < j; ++r) {
            double denom = right[r + 1] + left[j - r];
            double tmp   = (denom != 0.0) ? N[r] / denom : 0.0;
            N[r] = saved + right[r + 1] * tmp;
            saved = left[j - r] * tmp;
        }
        N[j] = saved;
    }
}

void qi_bs_spcol(const double *QI_RESTRICT U, int nknots, int p,
                 const double *QI_RESTRICT pts, int npts,
                 double *QI_RESTRICT M)
{
    int ndof = nknots - p - 1;
    /* zero M */
    memset(M, 0, (size_t)ndof * npts * sizeof(double));

    double Nbuf[16];
    for (int j = 0; j < npts; ++j) {
        int    span = qi_bs_findspan(ndof - 1, p, pts[j], U);
        qi_bs_basisfun(span, pts[j], p, U, Nbuf);
        /* the p+1 non-zero basis functions are at indices span-p .. span */
        for (int r = 0; r <= p; ++r) {
            int row = span - p + r;
            if (row >= 0 && row < ndof)
                M[row*npts + j] = Nbuf[r];
        }
    }
}

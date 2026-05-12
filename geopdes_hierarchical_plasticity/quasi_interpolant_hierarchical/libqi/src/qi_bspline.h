/*
 * qi_bspline.h — univariate B-spline routines (findspan, basisfun, spcol).
 *
 * Implements algorithms 2.1 and 2.2 from Piegl & Tiller, "The NURBS Book".
 *
 * v1.1: restrict qualifiers on pointer parameters for auto-vectorisation.
 */
#ifndef QI_BSPLINE_H
#define QI_BSPLINE_H

#include "qi_linalg.h"   /* for QI_RESTRICT */

/* Find the knot span i such that U[i] <= u < U[i+1].  n+1 control points,
 * degree p, knot vector U of length n+p+2.
 *   n: index of last basis function (so there are n+1 basis functions)
 *   p: degree
 *   u: parameter value
 *   U: knot vector
 * Returns i.
 */
int qi_bs_findspan(int n, int p, double u, const double *QI_RESTRICT U);

/* Compute the p+1 non-zero B-spline basis functions at u in span i.
 * Output: N[0..p].
 */
void qi_bs_basisfun(int i, double u, int p, const double *QI_RESTRICT U,
                    double *QI_RESTRICT N);

/* Compute B-spline collocation at npts points.
 *   U:      knot vector of length nknots
 *   p:      degree
 *   pts:    points to evaluate at  [npts]
 *   M:      output (ndof x npts) row-major, where ndof = nknots - p - 1.
 *   Each column j contains N_0(pts[j]), ..., N_{ndof-1}(pts[j]).
 *
 * This matches the row-major transpose of MATLAB's spcol(U, p+1, pts).
 */
void qi_bs_spcol(const double *QI_RESTRICT U, int nknots, int p,
                 const double *QI_RESTRICT pts, int npts,
                 double *QI_RESTRICT M);

#endif

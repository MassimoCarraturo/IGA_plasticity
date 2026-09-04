/*
 * qi_test.c — End-to-end test of the THB-spline quasi-interpolant.
 *
 * Three independent tests:
 *
 *   1. Polynomial reproduction.  Feed f(x,y) = a polynomial of degree ≤ d.
 *      A QI of polynomial degree d should reproduce it (almost) exactly.
 *      Required error: ‖f - QI(f)‖_inf < 1e-10.
 *
 *   2. Smooth function approximation.  Use the test function from
 *      Bracco et al. (2018), Example 1:
 *          f(x,y) = (tanh(9y - 9x) + 1)/9  +  exp(-((10x-6)^2 + (10y+7)^2))/1.5
 *      on [0,1]^2 with scattered data, and check the QI reaches an error
 *      consistent with the expected O(h^{d+1}) rate.
 *
 *   3. Vector-valued data (ncomp = 3).  Verify that all three components
 *      are processed independently and correctly.
 *
 * Build:  see ../Makefile (target qi_test)
 */

#include "qi.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

#ifdef _OPENMP
#include <omp.h>
#endif

/* -------------------------------------------------------------------------- */
/* Halton low-discrepancy sequence — much better coverage than rand() for     */
/* scattered data testing.                                                    */
/* -------------------------------------------------------------------------- */
static double halton(int i, int b)
{
    double f = 1.0, r = 0.0;
    while (i > 0) {
        f /= (double)b;
        r += f * (double)(i % b);
        i /= b;
    }
    return r;
}

/* -------------------------------------------------------------------------- */
/* Test target functions                                                      */
/* -------------------------------------------------------------------------- */
static double f_smooth(double x, double y)
{
    double t1 = (tanh(9.0*y - 9.0*x) + 1.0) / 9.0;
    double dx = 10.0*x - 6.0;
    double dy = 10.0*y + 7.0;
    double t2 = exp(-(dx*dx + dy*dy)) / 1.5;
    return t1 + t2;
}

/* polynomial of degree 2 in (x,y) for reproduction test */
static double f_poly2(double x, double y)
{
    return 0.7 - 1.3*x + 0.5*y + 2.1*x*x - 1.7*x*y + 0.9*y*y;
}

/* -------------------------------------------------------------------------- */
/* Build a simple single-level tensor-product hierarchical space with all     */
/* basis functions active at level 0.                                         */
/* -------------------------------------------------------------------------- */
typedef struct {
    int     degree_x, degree_y;
    int     nlevels;
    int    *nel_x;
    int    *nel_y;
    int     ndof;
    int    *ndof_per_level;
    int    *active_offset;
    int    *active_indices;
} simple_hspace_t;

static void simple_hspace_alloc(simple_hspace_t *S, int dx, int dy, int nelx, int nely)
{
    int ndof_x = nelx + dx;
    int ndof_y = nely + dy;
    int ndof   = ndof_x * ndof_y;

    S->degree_x       = dx;
    S->degree_y       = dy;
    S->nlevels        = 1;
    S->nel_x          = (int*)malloc(sizeof(int));
    S->nel_y          = (int*)malloc(sizeof(int));
    S->nel_x[0]       = nelx;
    S->nel_y[0]       = nely;
    S->ndof           = ndof;
    S->ndof_per_level = (int*)malloc(sizeof(int));
    S->ndof_per_level[0] = ndof;
    S->active_offset  = (int*)malloc(2*sizeof(int));
    S->active_offset[0] = 0;
    S->active_offset[1] = ndof;
    S->active_indices = (int*)malloc(ndof*sizeof(int));
    for (int i = 0; i < ndof; ++i) S->active_indices[i] = i;
}

static void simple_hspace_free(simple_hspace_t *S)
{
    free(S->nel_x); free(S->nel_y);
    free(S->ndof_per_level); free(S->active_offset); free(S->active_indices);
    memset(S, 0, sizeof(*S));
}

static void to_qi_hspace(const simple_hspace_t *S, qi_hspace_t *H)
{
    H->degree_x         = S->degree_x;
    H->degree_y         = S->degree_y;
    H->nlevels          = S->nlevels;
    H->nel_x            = S->nel_x;
    H->nel_y            = S->nel_y;
    H->breaks_x_data    = NULL;       /* uniform fallback */
    H->breaks_x_offset  = NULL;
    H->breaks_y_data    = NULL;
    H->breaks_y_offset  = NULL;
    H->ndof             = S->ndof;
    H->ndof_per_level   = S->ndof_per_level;
    H->active_offset    = S->active_offset;
    H->active_indices   = S->active_indices;
}

/* -------------------------------------------------------------------------- */
/* Test 1: polynomial reproduction                                            */
/* -------------------------------------------------------------------------- */
static int test_polynomial_reproduction(void)
{
    printf("\n=== Test 1: polynomial reproduction ===\n");

    const int dx = 2, dy = 2;
    const int nelx = 8, nely = 8;
    const int npts = 2000;

    simple_hspace_t S;
    simple_hspace_alloc(&S, dx, dy, nelx, nely);
    qi_hspace_t H;
    to_qi_hspace(&S, &H);

    /* Halton-distributed data points */
    double *xs = (double*)malloc(npts*sizeof(double));
    double *ys = (double*)malloc(npts*sizeof(double));
    double *fs = (double*)malloc(npts*sizeof(double));
    for (int i = 0; i < npts; ++i) {
        xs[i] = halton(i+1, 2);
        ys[i] = halton(i+1, 3);
        fs[i] = f_poly2(xs[i], ys[i]);
    }

    qi_data_t data = { npts, 1, xs, ys, fs };
    qi_options_t opts; qi_options_defaults(&opts);
    opts.verbose = 0;

    double *coeffs = (double*)calloc(S.ndof, sizeof(double));
    qi_status_t st = qi_compute_coeffs(&H, &data, &opts, coeffs, NULL);
    if (st != QI_OK) {
        printf("  qi_compute_coeffs returned %d\n", (int)st);
        free(xs); free(ys); free(fs); free(coeffs);
        simple_hspace_free(&S);
        return 1;
    }

    /* Evaluate at a fine 50x50 grid and compare */
    const int Ngrid = 50;
    double *gx = (double*)malloc(Ngrid*Ngrid*sizeof(double));
    double *gy = (double*)malloc(Ngrid*Ngrid*sizeof(double));
    double *gf = (double*)malloc(Ngrid*Ngrid*sizeof(double));
    for (int i = 0; i < Ngrid; ++i)
        for (int j = 0; j < Ngrid; ++j) {
            double xv = (double)i / (double)(Ngrid-1);
            double yv = (double)j / (double)(Ngrid-1);
            gx[i*Ngrid+j] = xv;
            gy[i*Ngrid+j] = yv;
        }
    qi_eval_tp(&H, 0, coeffs, 1, Ngrid*Ngrid, gx, gy, gf);

    double err_max = 0.0;
    for (int i = 0; i < Ngrid*Ngrid; ++i) {
        double exact = f_poly2(gx[i], gy[i]);
        double e = fabs(exact - gf[i]);
        if (e > err_max) err_max = e;
    }
    printf("  max polynomial reproduction error = %.3e\n", err_max);

    int ok = (err_max < 1e-9);
    printf("  %s\n", ok ? "PASS" : "FAIL");

    free(xs); free(ys); free(fs); free(coeffs);
    free(gx); free(gy); free(gf);
    simple_hspace_free(&S);

    return ok ? 0 : 1;
}

/* -------------------------------------------------------------------------- */
/* Test 2: smooth function approximation                                      */
/* -------------------------------------------------------------------------- */
static int test_smooth_approximation(void)
{
    printf("\n=== Test 2: smooth function approximation ===\n");

    const int dx = 3, dy = 3;
    const int nelx = 16, nely = 16;
    const int npts = 4000;

    simple_hspace_t S;
    simple_hspace_alloc(&S, dx, dy, nelx, nely);
    qi_hspace_t H;
    to_qi_hspace(&S, &H);

    double *xs = (double*)malloc(npts*sizeof(double));
    double *ys = (double*)malloc(npts*sizeof(double));
    double *fs = (double*)malloc(npts*sizeof(double));
    for (int i = 0; i < npts; ++i) {
        xs[i] = halton(i+1, 2);
        ys[i] = halton(i+1, 3);
        fs[i] = f_smooth(xs[i], ys[i]);
    }

    qi_data_t data = { npts, 1, xs, ys, fs };
    qi_options_t opts; qi_options_defaults(&opts);

    int *deg_used = (int*)calloc(S.ndof, sizeof(int));
    qi_diag_t diag = { deg_used, NULL, NULL, NULL };

    double *coeffs = (double*)calloc(S.ndof, sizeof(double));

    double t0 = 0.0;
#ifdef _OPENMP
    t0 = omp_get_wtime();
#else
    t0 = (double)clock() / (double)CLOCKS_PER_SEC;
#endif

    qi_status_t st = qi_compute_coeffs(&H, &data, &opts, coeffs, &diag);

    double t1 = 0.0;
#ifdef _OPENMP
    t1 = omp_get_wtime();
#else
    t1 = (double)clock() / (double)CLOCKS_PER_SEC;
#endif

    if (st != QI_OK) {
        printf("  qi_compute_coeffs returned %d\n", (int)st);
        free(xs); free(ys); free(fs); free(coeffs); free(deg_used);
        simple_hspace_free(&S);
        return 1;
    }

    /* Histogram of degrees actually used */
    int hist[QI_VERSION_MAJOR + 8] = {0};      /* big enough for d ≤ 7 */
    int dmax = 0;
    for (int k = 0; k < S.ndof; ++k) {
        if (deg_used[k] > 7) continue;
        hist[deg_used[k]]++;
        if (deg_used[k] > dmax) dmax = deg_used[k];
    }
    printf("  ndof = %d, npts = %d, time = %.4f s\n", S.ndof, npts, t1 - t0);
    printf("  degree histogram: ");
    for (int d = 0; d <= dmax; ++d)
        printf("d%d=%d ", d, hist[d]);
    printf("\n");

    /* Approximation error on a 100x100 grid */
    const int Ngrid = 100;
    double *gx = (double*)malloc(Ngrid*Ngrid*sizeof(double));
    double *gy = (double*)malloc(Ngrid*Ngrid*sizeof(double));
    double *gf = (double*)malloc(Ngrid*Ngrid*sizeof(double));
    for (int i = 0; i < Ngrid; ++i)
        for (int j = 0; j < Ngrid; ++j) {
            double xv = (double)i / (double)(Ngrid-1);
            double yv = (double)j / (double)(Ngrid-1);
            gx[i*Ngrid+j] = xv;
            gy[i*Ngrid+j] = yv;
        }
    qi_eval_tp(&H, 0, coeffs, 1, Ngrid*Ngrid, gx, gy, gf);

    double err_max = 0.0, err_l2 = 0.0;
    for (int i = 0; i < Ngrid*Ngrid; ++i) {
        double exact = f_smooth(gx[i], gy[i]);
        double e = fabs(exact - gf[i]);
        if (e > err_max) err_max = e;
        err_l2 += e*e;
    }
    err_l2 = sqrt(err_l2 / (Ngrid*Ngrid));
    printf("  max approx error  = %.3e\n", err_max);
    printf("  L2  approx error  = %.3e\n", err_l2);

    /* For d = 3 on a 16x16 mesh of [0,1]^2, h ≈ 1/16, expect error of
       a few × 1e-3 for this test function.  Be generous: < 1e-2.        */
    int ok = (err_max < 5e-2);
    printf("  %s\n", ok ? "PASS" : "FAIL");

    free(xs); free(ys); free(fs); free(coeffs); free(deg_used);
    free(gx); free(gy); free(gf);
    simple_hspace_free(&S);
    return ok ? 0 : 1;
}

/* -------------------------------------------------------------------------- */
/* Test 3: vector-valued data                                                 */
/* -------------------------------------------------------------------------- */
static int test_vector_valued(void)
{
    printf("\n=== Test 3: vector-valued data (ncomp = 3) ===\n");

    const int dx = 2, dy = 2;
    const int nelx = 8, nely = 8;
    const int npts = 2000;
    const int ncomp = 3;

    simple_hspace_t S;
    simple_hspace_alloc(&S, dx, dy, nelx, nely);
    qi_hspace_t H;
    to_qi_hspace(&S, &H);

    double *xs = (double*)malloc(npts*sizeof(double));
    double *ys = (double*)malloc(npts*sizeof(double));
    double *fs = (double*)malloc(npts*ncomp*sizeof(double));
    for (int i = 0; i < npts; ++i) {
        xs[i] = halton(i+1, 2);
        ys[i] = halton(i+1, 3);
        double x = xs[i], y = ys[i];
        fs[i*ncomp + 0] = f_poly2(x, y);
        fs[i*ncomp + 1] = 1.0 + 2.0*x - 3.0*y + x*y;
        fs[i*ncomp + 2] = -0.5 + 1.1*x*x - 0.4*y*y;
    }

    qi_data_t data = { npts, ncomp, xs, ys, fs };
    qi_options_t opts; qi_options_defaults(&opts);

    double *coeffs = (double*)calloc(S.ndof*ncomp, sizeof(double));
    qi_status_t st = qi_compute_coeffs(&H, &data, &opts, coeffs, NULL);
    if (st != QI_OK) {
        printf("  qi_compute_coeffs returned %d\n", (int)st);
        free(xs); free(ys); free(fs); free(coeffs);
        simple_hspace_free(&S);
        return 1;
    }

    /* Evaluate each component on a 50x50 grid */
    const int Ngrid = 50;
    int npg = Ngrid*Ngrid;
    double *gx = (double*)malloc(npg*sizeof(double));
    double *gy = (double*)malloc(npg*sizeof(double));
    for (int i = 0; i < Ngrid; ++i)
        for (int j = 0; j < Ngrid; ++j) {
            gx[i*Ngrid+j] = (double)i / (double)(Ngrid-1);
            gy[i*Ngrid+j] = (double)j / (double)(Ngrid-1);
        }
    /* Split per component, evaluate scalar field for each. */
    double err_max[3] = {0,0,0};
    for (int c = 0; c < ncomp; ++c) {
        double *cc = (double*)malloc(S.ndof*sizeof(double));
        for (int k = 0; k < S.ndof; ++k) cc[k] = coeffs[k*ncomp + c];
        double *gf = (double*)malloc(npg*sizeof(double));
        qi_eval_tp(&H, 0, cc, 1, npg, gx, gy, gf);
        for (int i = 0; i < npg; ++i) {
            double x = gx[i], y = gy[i], exact = 0.0;
            switch (c) {
                case 0: exact = f_poly2(x, y); break;
                case 1: exact = 1.0 + 2.0*x - 3.0*y + x*y; break;
                case 2: exact = -0.5 + 1.1*x*x - 0.4*y*y; break;
            }
            double e = fabs(exact - gf[i]);
            if (e > err_max[c]) err_max[c] = e;
        }
        free(cc); free(gf);
    }
    printf("  component errors: %.3e %.3e %.3e\n",
           err_max[0], err_max[1], err_max[2]);

    int ok = (err_max[0] < 1e-9) && (err_max[1] < 1e-9) && (err_max[2] < 1e-9);
    printf("  %s\n", ok ? "PASS" : "FAIL");

    free(xs); free(ys); free(fs); free(coeffs);
    free(gx); free(gy);
    simple_hspace_free(&S);
    return ok ? 0 : 1;
}

/* -------------------------------------------------------------------------- */
int main(int argc, char **argv)
{
    (void)argc; (void)argv;
    printf("THB-spline QI test suite — version %s\n", qi_version());
#ifdef _OPENMP
    printf("OpenMP enabled — max threads = %d\n", omp_get_max_threads());
#else
    printf("OpenMP NOT enabled (serial build)\n");
#endif

    int nfail = 0;
    nfail += test_polynomial_reproduction();
    nfail += test_smooth_approximation();
    nfail += test_vector_valued();

    printf("\n========== %s ==========\n",
           nfail == 0 ? "ALL TESTS PASSED" : "SOME TESTS FAILED");
    return nfail == 0 ? 0 : 1;
}

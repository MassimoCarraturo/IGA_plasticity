/*
 * qi_bench.c — quick parallel-scaling benchmark.
 *
 * Builds a single-level d=3 space with a large mesh, ~50k scattered points,
 * and times qi_compute_coeffs at several thread counts.
 */
#include "qi.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#ifdef _OPENMP
#include <omp.h>
#endif

static double halton(int i, int b) {
    double f = 1.0, r = 0.0;
    while (i > 0) { f /= b; r += f * (i % b); i /= b; }
    return r;
}

int main(void)
{
    const int dx = 3, dy = 3, nelx = 64, nely = 64, npts = 50000;
    int ndof_x = nelx + dx, ndof_y = nely + dy, ndof = ndof_x * ndof_y;

    int nlevels = 1;
    int nel_x[1] = { nelx }, nel_y[1] = { nely };
    int ndof_per_level[1] = { ndof };
    int *active_offset = malloc(2*sizeof(int));
    active_offset[0] = 0; active_offset[1] = ndof;
    int *active = malloc(ndof*sizeof(int));
    for (int i = 0; i < ndof; ++i) active[i] = i;

    qi_hspace_t H = {0};
    H.degree_x = dx; H.degree_y = dy; H.nlevels = nlevels;
    H.nel_x = nel_x; H.nel_y = nel_y;
    H.ndof = ndof;
    H.ndof_per_level = ndof_per_level;
    H.active_offset = active_offset;
    H.active_indices = active;

    double *xs = malloc(npts*sizeof(double));
    double *ys = malloc(npts*sizeof(double));
    double *fs = malloc(npts*sizeof(double));
    for (int i = 0; i < npts; ++i) {
        xs[i] = halton(i+1, 2);
        ys[i] = halton(i+1, 3);
        double t1 = (tanh(9.0*ys[i] - 9.0*xs[i]) + 1.0)/9.0;
        double dx0 = 10.0*xs[i] - 6.0, dy0 = 10.0*ys[i] + 7.0;
        double t2 = exp(-(dx0*dx0 + dy0*dy0))/1.5;
        fs[i] = t1 + t2;
    }
    qi_data_t D = { npts, 1, xs, ys, fs };

    qi_options_t opts; qi_options_defaults(&opts);
    double *coeffs = calloc(ndof, sizeof(double));

    int thread_counts[] = {1, 2, 4, 8};
    int n_tc = sizeof(thread_counts)/sizeof(thread_counts[0]);

    printf("Benchmark: ndof=%d, npts=%d, d=%d\n", ndof, npts, dx);
    printf("threads      time(s)    speedup\n");
    double t1 = 0.0;
    for (int k = 0; k < n_tc; ++k) {
        int nt = thread_counts[k];
        opts.n_threads = nt;
        double t0 =
#ifdef _OPENMP
            omp_get_wtime();
#else
            (double)clock()/CLOCKS_PER_SEC;
#endif
        qi_compute_coeffs(&H, &D, &opts, coeffs, NULL);
        double te =
#ifdef _OPENMP
            omp_get_wtime();
#else
            (double)clock()/CLOCKS_PER_SEC;
#endif
        double dt = te - t0;
        if (k == 0) t1 = dt;
        printf("  %2d        %7.3f     %5.2fx\n", nt, dt, t1/dt);
    }

    free(coeffs);
    free(xs); free(ys); free(fs);
    free(active); free(active_offset);
    return 0;
}

/*
 * qi_linalg.h — small dense linear algebra used internally by libqi.
 *
 * The matrices treated here are small (typically <= 50x15 for the LS design
 * and <= 25x25 for the basis-change system), so we use:
 *   - Householder QR for least squares
 *   - Jacobi rotation eigendecomposition for symmetric matrices (used to get
 *     the smallest singular value of A via A^T A)
 *   - LU with partial pivoting for square solves
 *
 * All matrices are stored in row-major C order: M[i*nc + j] is row i col j.
 */
#ifndef QI_LINALG_H
#define QI_LINALG_H

#include <stddef.h>

/* Compute x = pinv(A) * b   for A (m x n, m >= n) — minimum norm LS solution.
 * Internally uses A^T A eigendecomposition (Jacobi).
 *   A:        [m * n]   row-major (overwritten)
 *   b:        [m * nrhs] row-major
 *   x:        [n * nrhs] row-major (output)
 *   sigma_min:  returns smallest singular value of A (>=0)
 *   work:     scratch space, must be at least (n*n + n + 4*n + m*nrhs) doubles
 *   iwork:    scratch space, must be at least n ints
 *   Returns 0 on success, non-zero if A is severely rank-deficient.
 */
int qi_la_pinv_solve(int m, int n, int nrhs,
                     double *A, const double *b,
                     double *x, double *sigma_min,
                     double *work, int lwork);

/* Solve a square system L*x = b, where L is n x n, in-place LU with partial
 * pivoting.  L is overwritten.
 *   L:    [n*n]
 *   b:    [n*nrhs]   becomes x in place  (row-major)
 *   piv:  [n]        scratch
 * Returns 0 on success, non-zero if singular.
 */
int qi_la_lu_solve(int n, int nrhs, double *L, double *b, int *piv);

/* Symmetric eigendecomposition of an n x n matrix S (small).  Uses Jacobi
 * rotations.  S is overwritten with eigenvalues on the diagonal during the
 * sweeps; on return, S is destroyed and `evals` contains eigenvalues sorted
 * in ascending order; `Q` (n x n) contains the corresponding eigenvectors as
 * columns.
 *
 *   S:     [n*n]  (input/destroyed)
 *   evals: [n]    (output)
 *   Q:     [n*n]  (output, eigenvectors as columns, row-major)
 * Returns 0 on success.
 */
int qi_la_sym_eig(int n, double *S, double *evals, double *Q);

/* GEMM helpers — naive triple loops, fine for small sizes. */
void qi_la_gemm_NN(int m, int k, int n,
                   const double *A, const double *B, double *C);
void qi_la_gemm_TN(int m, int k, int n,
                   const double *A, const double *B, double *C);

#endif

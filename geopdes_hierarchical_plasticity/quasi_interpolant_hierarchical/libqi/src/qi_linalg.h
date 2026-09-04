/*
 * qi_linalg.h — small dense linear algebra used internally by libqi.
 *
 * The matrices treated here are small (typically <= 50x15 for the LS design
 * and <= 25x25 for the basis-change system), so we use:
 *   - Jacobi rotation eigendecomposition for symmetric matrices (used to get
 *     the smallest singular value of A via A^T A)
 *   - LU with partial pivoting for square solves
 *
 * All matrices are stored in row-major C order: M[i*nc + j] is row i col j.
 *
 * Optimisation notes (v1.1):
 *   - restrict qualifiers on all pointer parameters for auto-vectorisation.
 *   - LU solve uses a stack buffer to avoid per-call malloc.
 */
#ifndef QI_LINALG_H
#define QI_LINALG_H

#include <stddef.h>

/* MSVC uses __restrict instead of restrict in C mode */
#if defined(_MSC_VER) && !defined(__cplusplus)
#define QI_RESTRICT __restrict
#elif defined(__STDC_VERSION__) && __STDC_VERSION__ >= 199901L
#define QI_RESTRICT restrict
#else
#define QI_RESTRICT
#endif

/* Compute x = pinv(A) * b  for tall (m >= n) A — minimum-norm LS solution.
 * Internally uses A^T A eigendecomposition (Jacobi).
 *   A:        [m * n]    row-major (read-only in practice)
 *   b:        [m * nrhs] row-major
 *   x:        [n * nrhs] row-major (output)
 *   sigma_min:  returns sqrt(smallest eigenvalue of A^T A)
 *   work:     scratch, >= 2*n*n + n + 2*n*nrhs doubles
 * Returns 0 on success.
 */
int qi_la_pinv_solve(int m, int n, int nrhs,
                     double *QI_RESTRICT A,
                     const double *QI_RESTRICT b,
                     double *QI_RESTRICT x,
                     double *sigma_min,
                     double *QI_RESTRICT work, int lwork);

/* Solve a square system L*x = b in-place (LU with partial pivoting).
 *   L:   [n*n]     row-major — overwritten with LU factors
 *   b:   [n*nrhs]  row-major — overwritten with solution x
 *   piv: [n]       scratch
 * Returns 0 on success, non-zero if singular.
 */
int qi_la_lu_solve(int n, int nrhs,
                   double *QI_RESTRICT L,
                   double *QI_RESTRICT b,
                   int *QI_RESTRICT piv);

/* Symmetric eigendecomposition (Jacobi rotations). n is small.
 *   S:      [n*n]   input, destroyed
 *   evals:  [n]     output, ascending
 *   Q:      [n*n]   output, eigenvectors as columns (row-major)
 */
int qi_la_sym_eig(int n, double *S, double *evals, double *Q);

/* C = A * B    (m,k)*(k,n) -> (m,n)   row-major */
void qi_la_gemm_NN(int m, int k, int n,
                   const double *QI_RESTRICT A,
                   const double *QI_RESTRICT B,
                   double *QI_RESTRICT C);

/* C = A^T * B  where A is (k,m)       row-major */
void qi_la_gemm_TN(int m, int k, int n,
                   const double *QI_RESTRICT A,
                   const double *QI_RESTRICT B,
                   double *QI_RESTRICT C);

#endif

/*
 * qi_linalg.c — small dense linear algebra: Jacobi eig, pinv-via-eig, LU.
 *
 * All routines work on matrices stored in row-major C order.
 *
 * v1.1 changes:
 *   - restrict qualifiers for auto-vectorisation.
 *   - LU solve uses a stack buffer (no malloc for n*nrhs <= 512).
 */
#include "qi_linalg.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>

/* ---------- naive matrix multiplications --------------------------------- */
void qi_la_gemm_NN(int m, int k, int n,
                   const double *QI_RESTRICT A,
                   const double *QI_RESTRICT B,
                   double *QI_RESTRICT C)
{
    /* C = A * B  with shapes (m,k)*(k,n) -> (m,n) */
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            double s = 0.0;
            for (int p = 0; p < k; ++p) s += A[i*k + p] * B[p*n + j];
            C[i*n + j] = s;
        }
    }
}

void qi_la_gemm_TN(int m, int k, int n,
                   const double *QI_RESTRICT A,
                   const double *QI_RESTRICT B,
                   double *QI_RESTRICT C)
{
    /* C = A^T * B  where A is (k,m), so C = (m,n).  B is (k,n). */
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            double s = 0.0;
            for (int p = 0; p < k; ++p) s += A[p*m + i] * B[p*n + j];
            C[i*n + j] = s;
        }
    }
}

/* ---------- Symmetric eigendecomposition via Jacobi ---------------------- *
 *  Classical cyclic Jacobi.  Robust for small n; converges quickly.
 */
int qi_la_sym_eig(int n, double *S, double *evals, double *Q)
{
    const int    max_sweeps = 60;
    const double tol = 1e-14;

    /* initialise Q to the identity */
    for (int i = 0; i < n*n; ++i) Q[i] = 0.0;
    for (int i = 0; i < n;   ++i) Q[i*n + i] = 1.0;

    for (int sweep = 0; sweep < max_sweeps; ++sweep) {
        /* off-diagonal Frobenius norm */
        double off = 0.0;
        for (int p = 0; p < n; ++p)
            for (int q = p + 1; q < n; ++q) off += S[p*n + q] * S[p*n + q];
        if (off < tol*tol) break;

        for (int p = 0; p < n - 1; ++p) {
            for (int q = p + 1; q < n; ++q) {
                double Spq = S[p*n + q];
                if (fabs(Spq) < 1e-300) continue;

                double Spp = S[p*n + p];
                double Sqq = S[q*n + q];
                double theta = (Sqq - Spp) / (2.0 * Spq);
                double t;
                if (fabs(theta) > 1e150) {
                    t = 0.5 / theta;
                } else {
                    double sgn = (theta >= 0.0) ? 1.0 : -1.0;
                    t = sgn / (fabs(theta) + sqrt(theta*theta + 1.0));
                }
                double c = 1.0 / sqrt(t*t + 1.0);
                double s = t * c;

                /* update S: rotate rows/cols p, q */
                S[p*n + p] = Spp - t*Spq;
                S[q*n + q] = Sqq + t*Spq;
                S[p*n + q] = 0.0;
                S[q*n + p] = 0.0;
                for (int i = 0; i < n; ++i) {
                    if (i == p || i == q) continue;
                    double Sip = S[i*n + p];
                    double Siq = S[i*n + q];
                    S[i*n + p] = c*Sip - s*Siq;
                    S[i*n + q] = s*Sip + c*Siq;
                    S[p*n + i] = S[i*n + p];
                    S[q*n + i] = S[i*n + q];
                }
                /* update Q (eigenvectors as columns): rotate cols p, q */
                for (int i = 0; i < n; ++i) {
                    double Qip = Q[i*n + p];
                    double Qiq = Q[i*n + q];
                    Q[i*n + p] = c*Qip - s*Qiq;
                    Q[i*n + q] = s*Qip + c*Qiq;
                }
            }
        }
    }

    /* extract diagonal as eigenvalues, then sort ascending (carry Q) */
    for (int i = 0; i < n; ++i) evals[i] = S[i*n + i];

    /* simple selection sort, n is tiny */
    for (int i = 0; i < n - 1; ++i) {
        int    imin = i;
        double vmin = evals[i];
        for (int j = i + 1; j < n; ++j)
            if (evals[j] < vmin) { vmin = evals[j]; imin = j; }
        if (imin != i) {
            double tmp = evals[i]; evals[i] = evals[imin]; evals[imin] = tmp;
            for (int r = 0; r < n; ++r) {
                tmp = Q[r*n + i]; Q[r*n + i] = Q[r*n + imin]; Q[r*n + imin] = tmp;
            }
        }
    }
    return 0;
}

/* ---------- Pseudo-inverse-times-vector via eigendecomposition ----------- */
int qi_la_pinv_solve(int m, int n, int nrhs,
                     double *QI_RESTRICT A,
                     const double *QI_RESTRICT b,
                     double *QI_RESTRICT x,
                     double *sigma_min,
                     double *QI_RESTRICT work, int lwork)
{
    int need = 2*n*n + n + 2*n*nrhs;
    if (lwork < need) return -1;

    double *S      = work;
    double *Q      = S + n*n;
    double *evals  = Q + n*n;
    double *AtB    = evals + n;
    double *VtAtB  = AtB + n*nrhs;

    /* S = A^T A */
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            double s = 0.0;
            for (int k = 0; k < m; ++k) s += A[k*n + i] * A[k*n + j];
            S[i*n + j] = s;
        }
    }
    /* AtB = A^T b */
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < nrhs; ++j) {
            double s = 0.0;
            for (int k = 0; k < m; ++k) s += A[k*n + i] * b[k*nrhs + j];
            AtB[i*nrhs + j] = s;
        }
    }

    qi_la_sym_eig(n, S, evals, Q);

    double lam_max = evals[n-1];
    double tol     = (lam_max > 0.0 ? lam_max : 1.0) * 1e-13 * (double)n;
    double lam_min = evals[0] > 0.0 ? evals[0] : 0.0;
    *sigma_min = sqrt(lam_min);

    qi_la_gemm_TN(n, n, nrhs, Q, AtB, VtAtB);

    for (int i = 0; i < n; ++i) {
        double lam = evals[i];
        double inv = (lam > tol) ? 1.0/lam : 0.0;
        for (int j = 0; j < nrhs; ++j) VtAtB[i*nrhs + j] *= inv;
    }

    qi_la_gemm_NN(n, n, nrhs, Q, VtAtB, x);
    return 0;
}

/* ---------- LU with partial pivoting ------------------------------------- */
int qi_la_lu_solve(int n, int nrhs,
                   double *QI_RESTRICT L,
                   double *QI_RESTRICT b,
                   int *QI_RESTRICT piv)
{
    /* in-place LU factorisation of L; row-major */
    for (int k = 0; k < n; ++k) piv[k] = k;

    for (int k = 0; k < n; ++k) {
        /* find pivot */
        double pmax = fabs(L[k*n + k]);
        int    irow = k;
        for (int i = k + 1; i < n; ++i) {
            double v = fabs(L[i*n + k]);
            if (v > pmax) { pmax = v; irow = i; }
        }
        if (pmax < 1e-300) return -1;
        if (irow != k) {
            for (int j = 0; j < n; ++j) {
                double tmp = L[k*n + j]; L[k*n + j] = L[irow*n + j]; L[irow*n + j] = tmp;
            }
            int t = piv[k]; piv[k] = piv[irow]; piv[irow] = t;
        }
        double pivot = L[k*n + k];
        for (int i = k + 1; i < n; ++i) {
            double f = L[i*n + k] / pivot;
            L[i*n + k] = f;
            for (int j = k + 1; j < n; ++j) L[i*n + j] -= f * L[k*n + j];
        }
    }

    /* apply permutation to b: b_perm[i] = b[piv[i]]   (row-major nrhs)
     * Use a stack buffer to avoid malloc in the hot parallel path. */
    double bp_stack[512];
    int    heap = (n * nrhs > 512);
    double *bp = heap ? (double*)malloc((size_t)n * nrhs * sizeof(double)) : bp_stack;
    if (!bp) return -2;

    for (int i = 0; i < n; ++i)
        for (int j = 0; j < nrhs; ++j)
            bp[i*nrhs + j] = b[piv[i]*nrhs + j];
    memcpy(b, bp, (size_t)n * nrhs * sizeof(double));

    if (heap) free(bp);

    /* forward substitution: L y = b   (unit diag) */
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < nrhs; ++j)
            for (int k = 0; k < i; ++k)
                b[i*nrhs + j] -= L[i*n + k] * b[k*nrhs + j];

    /* backward substitution: U x = y */
    for (int i = n - 1; i >= 0; --i) {
        double diag = L[i*n + i];
        for (int j = 0; j < nrhs; ++j) {
            double s = b[i*nrhs + j];
            for (int k = i + 1; k < n; ++k) s -= L[i*n + k] * b[k*nrhs + j];
            b[i*nrhs + j] = s / diag;
        }
    }
    return 0;
}

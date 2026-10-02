#include "linsys.h"
#include "def_solnp_lapack.h"

solnp_int SOLNP(chol)(
    solnp_int n,
    solnp_float *H)
{
    solnp_int info;
    char uploup = LAPACK_UPLOW_UP;

    dpotrf(&uploup, &n, H, &n, &info);
    // info = LAPACKE_dpotrf(LAPACK_COL_MAJOR, 'U', n, H, n);
    return info; // 0 if successful
}

solnp_int SOLNP(solve_lin_sys)(
    solnp_int n,
    solnp_int nrhs,
    solnp_float *L,
    solnp_float *b)
{ // This subroutine solve PSD linear system
    solnp_int info;
    char uploup = LAPACK_UPLOW_UP;
    dpotrs(&uploup, &n, &nrhs, L, &n, b, &n, &info);

    // info = LAPACKE_dpotrs(LAPACK_COL_MAJOR, 'U', n, nrhs, L, n, b, n);
    return info; // 0 if successful
}

void SOLNP(cond)(
    solnp_int n,
    solnp_float *a,
    solnp_float *cond)
{
    solnp_float norm = 0;
    for (solnp_int i = 0; i < n; i++)
    {
        for (solnp_int j = 0; i < n; i++)
        {
            norm += ABS(a[i * n + j]);
        }
    }
    // LAPACKE_dpotrf(LAPACK_COL_MAJOR, 'U', n, a, n);
    // LAPACKE_dpocon(LAPACK_COL_MAJOR, 'U', n, a, n, norm, cond);

    solnp_int info;
    char uploup = LAPACK_UPLOW_UP;
    solnp_float *work = (solnp_float *)solnp_malloc(3 * n * sizeof(solnp_float));
    solnp_int *iwork = (solnp_int *)solnp_malloc(n * sizeof(solnp_int));
    dpotrf(&uploup, &n, a, &n, &info);
    dpocon(&uploup, &n, a, &n, &norm, cond, work, iwork, &info);
    solnp_free(work);
    solnp_free(iwork);
}
// solve linear system of square matrix
void SOLNP(solve_general_lin_sys)(
    solnp_int n,
    solnp_float *a,
    solnp_float *b)
{
    solnp_int *ipiv = (solnp_int *)solnp_malloc(n * sizeof(solnp_int));
    // LAPACKE_dgetrf(LAPACK_COL_MAJOR, n, n, a, n, ipiv);
    // LAPACKE_dgetrs(LAPACK_COL_MAJOR, 'N', n, 1, a, n, ipiv, b, n);

    solnp_int info;
    char trans = LAPACK_NOTRANS;
    solnp_int nrhs = 1;

    dgetrf(&n, &n, a, &n, ipiv, &info);
    dgetrs(&trans, &n, &nrhs, a, &n, ipiv, b, &n, &info);

    solnp_free(ipiv);
    return;
}

solnp_int SOLNP(solve_sys_lin_sys)(
    solnp_int n,
    solnp_float *a,
    solnp_float *b,
    solnp_int max_neg_eig)
{
    // This subroutine solve general symmetric linear system
    solnp_int info, i, neg_eig;
    neg_eig = 0;
    solnp_int *ipiv = (solnp_int *)solnp_malloc(n * sizeof(solnp_int));
    // LAPACKE_dsytrf(LAPACK_COL_MAJOR, 'L', n, a, n, ipiv);

    char uplolow = LAPACK_UPLOW_LOW;
    solnp_int lwork = 8 * n;
    solnp_float *work = (solnp_float *)solnp_malloc(MAX(1, lwork) * sizeof(solnp_float));
    dsytrf(&uplolow, &n, a, &n, ipiv, work, &lwork, &info);
    solnp_free(work);

    for (i = 0; i < n; i++)
    {
        if (a[i + i * n] < 0)
        {
            neg_eig++;
        }
    }
    if (neg_eig > max_neg_eig)
    {
        return -1;
    }
    // info = LAPACKE_dsytrs(LAPACK_COL_MAJOR, 'L', n, 1, a, n, ipiv, b, n);
    dsytrs(&uplolow, &n, &n, a, &n, ipiv, b, &n, &info);

    solnp_free(ipiv);
    return info;
}

solnp_int SOLNP(least_square)(
    solnp_int m,
    solnp_int n,
    solnp_float *A,
    solnp_float *b)
{
    // Solve the least square problem using SVD decompoisiton
    solnp_int info, rank;
    solnp_float *s = (solnp_float *)solnp_malloc(n * sizeof(solnp_float));
    // info = LAPACKE_dgelss(LAPACK_COL_MAJOR, m, n, 1, A, m, b, MAX(m, n), s, 1e-8, &rank);

    solnp_int nrhs = 1;
    solnp_int ldb = MAX(m, n);
    solnp_float rcond = 1e-8;
    solnp_int min_lwork = 3 * MIN(m, n) + MAX(MAX(2 * MIN(m, n), MAX(m, n)), nrhs);
    solnp_int lwork = 8 * MAX(n, min_lwork);
    solnp_float *work = (solnp_float *)solnp_malloc(MAX(1, lwork) * sizeof(solnp_float));
    dgelss(&m, &n, &nrhs, A, &m, b, &ldb, s, &rcond, &rank, work, &lwork, &info);
    solnp_free(work);
    solnp_free(s);

    return info;
}

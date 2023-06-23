#include "linsys.h"

solnp_int SOLNP(chol)
(
    solnp_int n,
    solnp_float *H
)
{
    solnp_int info;

    info = LAPACKE_dpotrf(LAPACK_COL_MAJOR, 'U',  n, H, n);
    return info; // 0 if successful

}

solnp_int SOLNP(solve_lin_sys)
(
    solnp_int n,
    solnp_int nrhs,
    solnp_float *L,
    solnp_float *b
)
{   // This subroutine solve PSD linear system
    solnp_int info;

    info = LAPACKE_dpotrs(LAPACK_COL_MAJOR, 'U', n, nrhs, L, n, b, n);

    return info; // 0 if successful
}

void SOLNP(cond)
(
    solnp_int n,
    solnp_float* a,
    solnp_float* cond
    )
{
    solnp_float norm = 0;
    for (solnp_int i = 0; i < n; i++) {
        for (solnp_int j = 0; i < n; i++) {
            norm += ABS(a[i * n + j]);
        }
    }
    LAPACKE_dpotrf(LAPACK_COL_MAJOR, 'U', n, a, n);
    LAPACKE_dpocon(LAPACK_COL_MAJOR, 'U', n, a, n, norm, cond);
}
//solve linear system of square matrix
void SOLNP(solve_general_lin_sys)
(
    solnp_int n,
    solnp_float* a,
    solnp_float* b
){
    solnp_int* ipiv = (solnp_int*) solnp_malloc(n*sizeof(solnp_int));
    LAPACKE_dgetrf(LAPACK_COL_MAJOR, n, n, a, n, ipiv);
    LAPACKE_dgetrs(LAPACK_COL_MAJOR,'N', n, 1, a, n, ipiv, b, n);
    solnp_free(ipiv);
    return;
}

solnp_int SOLNP(solve_sys_lin_sys)
(
    solnp_int n,
    solnp_float* a,
    solnp_float* b,
    solnp_int max_neg_eig
 ){
    // This subroutine solve general symmetric linear system
    solnp_int state,i,neg_eig;
    neg_eig = 0;
    solnp_int* ipiv = (solnp_int*)solnp_malloc(n * sizeof(solnp_int));
    LAPACKE_dsytrf(LAPACK_COL_MAJOR,'L', n, a, n, ipiv);
    for (i = 0; i < n; i++) {
        if (a[i + i * n] < 0) {
            neg_eig++;
        }
    }
    if (neg_eig > max_neg_eig) {
        return -1;
    }
    state = LAPACKE_dsytrs(LAPACK_COL_MAJOR, 'L', n, 1, a, n, ipiv, b, n);
    solnp_free(ipiv);
    return state;
}
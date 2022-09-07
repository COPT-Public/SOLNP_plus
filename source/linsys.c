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
{
    solnp_int info;

    info = LAPACKE_dpotrs(LAPACK_COL_MAJOR, 'U', n, nrhs, L, n, b, n);

    return info; // 0 if successful
}
#include "solnp.h"
// #include "mkl_lapacke.h"

solnp_int SOLNP(chol)(
    solnp_int n,
    solnp_float *H);

solnp_int SOLNP(solve_lin_sys)(
    solnp_int n,
    solnp_int nrhs,
    solnp_float *L,
    solnp_float *b);
void SOLNP(cond)(
    solnp_int n,
    solnp_float *a,
    solnp_float *cond);
void SOLNP(solve_general_lin_sys)(
    solnp_int n,
    solnp_float *a,
    solnp_float *b);
solnp_int SOLNP(solve_sys_lin_sys)(
    solnp_int n,
    solnp_float *a,
    solnp_float *b,
    solnp_int max_neg_eig);
solnp_int SOLNP(least_square)(
    solnp_int m,
    solnp_int n,
    solnp_float *A,
    solnp_float *b);
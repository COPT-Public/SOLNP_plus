#include "solnp.h"
#include "solnp_util.h"
#include <stdio.h>

// void rosenbrock(solnp_float *x, solnp_float *result, solnp_int np, solnp_int nfeval)
void rosenbrock(solnp_float *x, solnp_float *result, solnp_int np)
{
    result[0] = (1.0 - x[0]) * (1.0 - x[0]) + 100.0 * (x[1] - x[0] * x[0]) * (x[1] - x[0] * x[0]);
}

int main(void)
{

    SOLNPSettings *stgs = (SOLNPSettings *)solnp_calloc(1, sizeof(SOLNPSettings));
    SOLNPSol *sol = (SOLNPSol *)solnp_calloc(1, sizeof(SOLNPSol));
    SOLNPInfo *info = (SOLNPInfo *)solnp_calloc(1, sizeof(SOLNPInfo));
    solnp_float *l = SOLNP_NULL;
    solnp_float *h = SOLNP_NULL;
    SOLNPProb *prob = SOLNP(init_prob)(2, 0, 0);

    SOLNP(set_default_settings)
    (stgs, prob->np);

    prob->pbl = (solnp_float *)solnp_malloc(prob->np * sizeof(solnp_float));
    prob->pbl[0] = -1.0;
    prob->pbl[1] = -1.0;

    prob->pbu = (solnp_float *)solnp_malloc(prob->np * sizeof(solnp_float));
    prob->pbu[0] = 2.0;
    prob->pbu[1] = 2.0;

    prob->p0 = (solnp_float *)solnp_malloc(prob->np * sizeof(solnp_float));
    memset(prob->p0, 0, prob->np * sizeof(solnp_float));

    SOLNP_PLUS(prob, stgs, sol, info, rosenbrock, SOLNP_NULL, SOLNP_NULL, l, h);

    printf("SOLNP_PLUS done in %d iterations, %f seconds\nThe objective is %f\n", sol->iter, info->total_time, sol->obj);

    // free the memory allocated to call solnp
    SOLNP(free_sol)
    (sol);
    SOLNP(free_info)
    (info);
    SOLNP(free_stgs)
    (stgs);
    SOLNP(free_prob)
    (prob);
    if (l)
    {
        solnp_free(l);
    }
    if (h)
    {
        solnp_free(h);
    }

    return 0;
}

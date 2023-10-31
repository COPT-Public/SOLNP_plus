#include "solnp_util.h"

// void cost(solnp_float *x, solnp_float *result, solnp_int np, solnp_int nfeval)
void cost(solnp_float *x, solnp_float *result, solnp_int np)
{
    result[0] = pow(x[0] - 5, 2) + pow(x[1], 2) - 25;
    result[1] = -pow(x[0], 2) + x[1];
}

int main(void)
{
    SOLNPSettings *stgs = (SOLNPSettings *)solnp_calloc(1, sizeof(SOLNPSettings));
    SOLNPSol *sol = (SOLNPSol *)solnp_calloc(1, sizeof(SOLNPSol));
    SOLNPInfo *info = (SOLNPInfo *)solnp_calloc(1, sizeof(SOLNPInfo));
    solnp_float *l = SOLNP_NULL;
    solnp_float *h = SOLNP_NULL;
    SOLNPProb *prob = SOLNP(init_prob)(2, 1, 0);

    SOLNP(set_default_settings)
    (stgs, prob->np);

    prob->p0 = (solnp_float *)solnp_malloc(prob->np * sizeof(solnp_float));
    prob->p0[0] = 4.9;
    prob->p0[1] = 0.1;

    prob->ib0 = (solnp_float *)solnp_malloc(prob->nic * sizeof(solnp_float));
    prob->ib0[0] = 1.0;

    prob->ibl = (solnp_float *)solnp_malloc(prob->nic * sizeof(solnp_float));
    prob->ibl[0] = 0.0;

    SOLNP_PLUS(prob, stgs, sol, info, cost, SOLNP_NULL, SOLNP_NULL, l, h);

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
#include "solnp.h"
#include "solnp_util.h"
#include <stdio.h>

// void p78(solnp_float *x, solnp_float *result, solnp_int np, solnp_int nfeval)
void p78(solnp_float *x, solnp_float *result, solnp_int np)
{
    result[0] = 1.0;
    for (int i = 0; i < np; ++i)
    {
        result[0] *= x[i];
    }
    result[1] = -10.0;
    for (int i = 0; i < np; ++i)
    {
        result[1] += x[i] * x[i];
    }
    result[2] = x[1] * x[2] - 5 * x[3] * x[4];
    result[3] = x[0] * x[0] * x[0] + x[1] * x[1] * x[1] + 1.0;
}

int main(void)
{

    SOLNPSettings *stgs = (SOLNPSettings *)solnp_calloc(1, sizeof(SOLNPSettings));
    SOLNPSol *sol = (SOLNPSol *)solnp_calloc(1, sizeof(SOLNPSol));
    SOLNPInfo *info = (SOLNPInfo *)solnp_calloc(1, sizeof(SOLNPInfo));
    solnp_float *l = SOLNP_NULL;
    solnp_float *h = SOLNP_NULL;
    SOLNPProb *prob = SOLNP(init_prob)(5, 0, 3);

    SOLNP(set_default_settings)
    (stgs, prob->np);

    prob->p0 = (solnp_float *)solnp_malloc(prob->np * sizeof(solnp_float));
    prob->p0[0] = -2.0;
    prob->p0[1] = 1.5;
    prob->p0[2] = 2.0;
    prob->p0[3] = -1.0;
    prob->p0[4] = -1.0;

    SOLNP_PLUS(prob, stgs, sol, info, p78, SOLNP_NULL, SOLNP_NULL, l, h);

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

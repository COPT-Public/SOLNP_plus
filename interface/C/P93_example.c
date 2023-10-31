#include "solnp.h"
#include "solnp_util.h"
#include <stdio.h>

// void p93(solnp_float *x, solnp_float *result, solnp_int np, solnp_int nfeval)
void p93(solnp_float *x, solnp_float *result, solnp_int np)
{
    result[0] = 0.0;
    result[0] += 0.0204 * x[0] * x[3] * (x[0] + x[1] + x[2]);
    result[0] += 0.0187 * x[1] * x[2] * (x[0] + 1.57 * x[1] + x[3]);
    result[0] += 0.0607 * x[0] * x[3] * x[4] * x[4] * (x[0] + x[1] + x[2]);
    result[0] += 0.0437 * x[1] * x[2] * x[5] * x[5] * (x[0] + 1.57 * x[1] + x[3]);
    result[1] = 0.001 * x[0] * x[1] * x[2] * x[3] * x[4] * x[5] - 2.07;
    result[2] = 1.0;
    result[2] -= 0.00062 * x[0] * x[3] * x[4] * x[4] * (x[0] + x[1] + x[2]);
    result[2] -= 0.00058 * x[1] * x[2] * x[5] * x[5] * (x[0] + 1.57 * x[1] + x[3]);
}

int main(void)
{

    SOLNPSettings *stgs = (SOLNPSettings *)solnp_calloc(1, sizeof(SOLNPSettings));
    SOLNPSol *sol = (SOLNPSol *)solnp_calloc(1, sizeof(SOLNPSol));
    SOLNPInfo *info = (SOLNPInfo *)solnp_calloc(1, sizeof(SOLNPInfo));
    solnp_float *l = SOLNP_NULL;
    solnp_float *h = SOLNP_NULL;
    SOLNPProb *prob = SOLNP(init_prob)(6, 2, 0);

    SOLNP(set_default_settings)
    (stgs, prob->np);

    prob->pbl = (solnp_float *)solnp_malloc(prob->np * sizeof(solnp_float));
    prob->ibl = (solnp_float *)solnp_malloc(prob->nic * sizeof(solnp_float));
    prob->ib0 = (solnp_float *)solnp_malloc(prob->nic * sizeof(solnp_float));
    prob->p0 = (solnp_float *)solnp_malloc(prob->np * sizeof(solnp_float));

    memset(prob->pbl, 0, prob->np * sizeof(solnp_float));
    memset(prob->ibl, 0, prob->nic * sizeof(solnp_float));
    prob->ib0[0] = 1.0;
    prob->ib0[1] = 1.0;
    prob->p0[0] = 5.54;
    prob->p0[1] = 4.4;
    prob->p0[2] = 12.02;
    prob->p0[3] = 11.82;
    prob->p0[4] = 0.702;
    prob->p0[5] = 0.852;

    SOLNP_PLUS(prob, stgs, sol, info, p93, SOLNP_NULL, SOLNP_NULL, l, h);

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

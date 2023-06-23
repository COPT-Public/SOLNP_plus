#include "solnp.h"
#include "linalg.h"
#include "solnp_util.h"

// solnp_float cost_time = 0; // global cost timer

typedef void cost_temple(double *p, double *result, int n);
typedef void g_temple(double *p, double *result);
typedef void h_temple(double *p, double *result);

cost_temple *py_cost;
g_temple *py_grad;
h_temple *py_hess;

void def_python_callback(cost_temple *cost, g_temple *grad, h_temple *hess)
{
    // Function called by Python once
    // Defines what "py_cost" is pointing to
    py_cost = cost;
    py_grad = grad;
    py_hess = hess;
}

void C_call_python_cost(SOLNPCost **c, solnp_float *p, solnp_int np, solnp_int nfeval)
{

    SOLNP(timer)
    cost_timer;
    SOLNP(tic)
    (&cost_timer);
    solnp_int i, j;

    double *result = (double *)solnp_malloc(sizeof(double) * (1 + c[0]->nic + c[0]->nec) * nfeval);
    py_cost(p, result, nfeval);

    for (j = 0; j < nfeval; j++)
    {
        c[j]->obj = (solnp_float)result[j * (1 + c[j]->nic + c[j]->nec)];

        for (i = 0; i < c[j]->nec; i++)
        {
            c[j]->ec[i] = (solnp_float)result[i + 1 + j * (1 + c[j]->nic + c[j]->nec)];
        }
        for (i = 0; i < c[j]->nic; i++)
        {
            c[j]->ic[i] = (solnp_float)result[i + 1 + c[j]->nec + j * (1 + c[j]->nic + c[j]->nec)];
        }
    }

    solnp_free(result);
    // cost_time += SOLNP(tocq)(&cost_timer) / 1e3;
}

void C_call_python_grad(solnp_float *g, solnp_float *p, solnp_int np, solnp_int ngeval)
{

    SOLNP(timer)
    cost_timer;
    SOLNP(tic)
    (&cost_timer);
    solnp_int i;

    double *result = (double *)solnp_malloc(sizeof(double) * np * ngeval);
    py_grad(p, result);

    for (i = 0; i < np * ngeval; i++)
    {
        g[i] = (solnp_float)result[i];
    }
    // cost_time += SOLNP(tocq)(&cost_timer) / 1e3;
}

void C_call_python_hess(solnp_float *h, solnp_float *p, solnp_int np, solnp_int nheval)
{

    SOLNP(timer)
    cost_timer;
    SOLNP(tic)
    (&cost_timer);
    solnp_int i;

    double *result = (double *)solnp_malloc(sizeof(double) * np * np * nheval);
    py_hess(p, result);

    for (i = 0; i < np * np * nheval; i++)
    {
        h[i] = (solnp_float)result[i];
    }
    // cost_time += SOLNP(tocq)(&cost_timer) / 1e3;
}

void SOLNP_C(
    // input
    solnp_float *ibl,
    solnp_float *ibu,
    solnp_float *pbl,
    solnp_float *pbu,
    solnp_float *Ipc,
    solnp_float *Ipb,
    solnp_float *ib0,
    solnp_float *p,
    solnp_float *op,
    solnp_float *l,
    solnp_float *h,
    solnp_int np,
    solnp_int nic,
    solnp_int nec,
    // output
    solnp_float *scalars,
    solnp_float *p_out,
    solnp_float *best_fea_p,
    solnp_float *ic,
    solnp_float *jh,
    solnp_float *ch,
    solnp_float *l_out,
    solnp_float *h_out)
{
    SOLNP(timer)
    total_timer;
    SOLNP(tic)
    (&total_timer);

    solnp_int index, i;
    /*--------------------- construct settings ------------------------*/
    SOLNPSettings *stgs = (SOLNPSettings *)solnp_malloc(sizeof(SOLNPSettings));
    set_default_settings(stgs);
    stgs->maxfev = 500 * np;
    index = 0;

    stgs->rho = (solnp_float)op[index++];
    stgs->pen_l1 = (solnp_float)op[index++];
    stgs->max_iter = (solnp_int)op[index++];
    stgs->min_iter = (solnp_int)op[index++];
    stgs->max_iter_rescue = (solnp_int)op[index++];
    stgs->min_iter_rescue = (solnp_int)op[index++];
    stgs->delta = (solnp_float)op[index++];
    stgs->tol = (solnp_float)op[index++];
    stgs->tol_con = (solnp_float)op[index++];
    stgs->ls_time = (solnp_int)op[index++];
    stgs->batchsize = (solnp_int)op[index++];
    stgs->tol_restart = (solnp_float)op[index++];
    stgs->re_time = (solnp_int)op[index++];
    stgs->delta_end = (solnp_float)op[index++];
    stgs->maxfev = (solnp_int)op[index++];
    stgs->noise = (solnp_int)op[index++];
    stgs->qpsolver = (solnp_int)op[index++];
    stgs->scale = (solnp_int)op[index++];
    stgs->bfgs = (solnp_int)op[index++];
    stgs->rs = (solnp_int)op[index++];
    stgs->grad = (solnp_int)op[index++];
    stgs->k_i = (solnp_float)op[index++];
    stgs->k_r = (solnp_float)op[index++];
    stgs->c_r = (solnp_float)op[index++];
    stgs->c_i = (solnp_float)op[index++];
    stgs->ls_way = (solnp_int)op[index++];
    stgs->rescue = (solnp_int)op[index++];
    stgs->drsom = (solnp_int)op[index++];

    /*--------------------- check function from python ------------------------*/
    if (py_cost == SOLNP_NULL)
    {
        solnp_printf("SOLNP+--> The user does not provided cost function in the fun structure. \n");
        solnp_printf("SOLNP+--> SOLNP stops.\n");
        exit(1);
    }

    if (py_grad != SOLNP_NULL)
    {
        stgs->grad = 1;
    }
    else
    {
        solnp_printf("SOLNP+--> The user does not provided gradient function of cost in the fun structure. \n");
        solnp_printf("SOLNP+--> SOLNP uses zero-order method instead.\n");
        if (stgs->rs)
        {
            solnp_printf("SOLNP+--> SOLNP uses random sampling method to estimate the gradient.\n");
        }
        else
        {
            stgs->grad = 0;
        }
    }

    if (py_hess != SOLNP_NULL)
    {
        stgs->hess = 1;
        solnp_printf("SOLNP+--> Using second-order method in optimization. \n");
    }
    else
    {
        if (stgs->grad && !stgs->rs)
        {
            solnp_printf("SOLNP+--> Using first-order method in optimization. \n");
        }
        stgs->hess = 0;
    }

    /*--------------------- construct constraint ------------------------*/
    SOLNPConstraint *constraint = malloc_constriant(np, nic);
    if (nic > 0)
    {
        memcpy(constraint->il, ibl, nic * sizeof(solnp_float));
        memcpy(constraint->iu, ibu, nic * sizeof(solnp_float));
    }

    if (np > 0)
    {
        memcpy(constraint->pl, pbl, np * sizeof(solnp_float));
        memcpy(constraint->pu, pbu, np * sizeof(solnp_float));
    }

    memcpy(constraint->Ipc, Ipc, 2 * sizeof(solnp_float));
    memcpy(constraint->Ipb, Ipb, 2 * sizeof(solnp_float));

    /*--------------------- assemble input ------------------------*/
    solnp_int nc = nec + nic;
    solnp_int h_len = (np + nic) * (np + nic);
    SOLNPIput *input = (SOLNPIput *)solnp_malloc(sizeof(SOLNPIput));
    input->cnstr = constraint;
    input->stgs = stgs;
    input->n = np;
    input->h = (solnp_float *)solnp_malloc(sizeof(solnp_float) * h_len);
    input->l = (solnp_float *)solnp_malloc(sizeof(solnp_float) * nc);
    memcpy(input->h, h, h_len * sizeof(solnp_float));
    memcpy(input->l, l, nc * sizeof(solnp_float));
    /*--------------------- construct cost ------------------------*/
    SOLNP(timer)
    cost_timer;
    SOLNP(tic)
    (&cost_timer);
    // first call python cost function
    double *result = (double *)solnp_malloc(sizeof(double) * (nec + nic + 1));
    py_cost(p, result, 1);
    // cost_time += SOLNP(tocq)(&cost_timer) / 1e3;

    solnp_int m = nec + nic + 1;
    solnp_float *ob = (solnp_float *)solnp_malloc(sizeof(solnp_float) * m);
    for (i = 0; i < m; i++)
    {
        ob[i] = (solnp_float)result[i];
    }

    SOLNPCost *cost = malloc_cost(nec, nic, &C_call_python_cost, &C_call_python_grad, &C_call_python_hess);
    cost->obj = ob[0];
    for (i = 0; i < nec; i++)
    {
        cost->ec[i] = ob[1 + i];
    }
    for (i = 0; i < nic; i++)
    {
        cost->ic[i] = ob[1 + nec + i];
    }
    solnp_free(ob);
    // solnp_free(result);
    /*--------------------- construct other params ------------------------*/
    SOLNPSol *sol = (SOLNPSol *)solnp_malloc(sizeof(SOLNPSol));
    SOLNPInfo *info = (SOLNPInfo *)solnp_malloc(sizeof(SOLNPInfo));

    solnp_float *ib0_p = (solnp_float *)solnp_malloc(sizeof(solnp_float) * (np + nic));

    if (nic > 0)
    {
        memcpy(ib0_p, ib0, nic * sizeof(solnp_float));
        // solnp_free(ib0);
    }

    memcpy(&ib0_p[nic], p, np * sizeof(solnp_float));

    // solnp_printf("Ipb[0]: %f\n", constraint->Ipb[0]);
    // solnp_printf("Ipb[1]: %f\n", constraint->Ipb[1]);
    // solnp_printf("Ipc[0]: %f\n", constraint->Ipc[0]);
    // solnp_printf("Ipc[1]: %f\n", constraint->Ipc[1]);

    /*--------------------- call solnp ------------------------*/
    solnp_int status = SOLNP(main)(input, cost, ib0_p, sol, info);
    /*--------------------- assemble output ------------------------*/
    index = 0;

    scalars[index++] = (solnp_float)sol->iter;
    scalars[index++] = (solnp_float)sol->count_cost;
    scalars[index++] = (solnp_float)sol->count_grad;
    scalars[index++] = (solnp_float)sol->count_hess;
    scalars[index++] = (solnp_float)sol->constraint;
    scalars[index++] = (solnp_float)sol->restart_time;
    scalars[index++] = (solnp_float)sol->obj;
    scalars[index++] = (solnp_float)sol->status;
    scalars[index++] = (solnp_float)info->total_time;

    memcpy(p_out, sol->p, np * sizeof(solnp_float));
    memcpy(best_fea_p, sol->best_fea_p, np * sizeof(solnp_float));
    memcpy(ic, sol->ic, MAX(nic, 1) * sizeof(solnp_float));
    memcpy(jh, sol->jh, (sol->iter + 1) * sizeof(solnp_float));
    memcpy(ch, sol->ch, (sol->iter + 1) * sizeof(solnp_float));
    memcpy(l_out, sol->l, MAX(nic + nec, 1) * sizeof(solnp_float));
    if (stgs->bfgs)
    {
        memcpy(h_out, sol->h, (np + nic) * (np + nic) * sizeof(solnp_float));
    }
    else
    {
        memcpy(h_out, sol->h, 1 * sizeof(solnp_float));
    }

    solnp_free(info);
    // free_sol(sol);
    solnp_printf("SOLNP+--> SOLNP finished.\n");
}
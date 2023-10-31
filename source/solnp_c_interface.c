#include "solnp.h"
#include "solnp_util.h"

cost_temple *cost_fun_c;
g_temple *grad_fun_c;
h_temple *hess_fun_c;

void def_c_callback(cost_temple *cost, g_temple *grad, h_temple *hess)
{
    cost_fun_c = cost;
    grad_fun_c = grad;
    hess_fun_c = hess;
}

void C_call_cost(SOLNPCost **c, solnp_float *p, solnp_int np, solnp_int nfeval)
{

    SOLNP(timer)
    cost_timer;
    SOLNP(tic)
    (&cost_timer);
    solnp_int i, j;

    double *result = (double *)solnp_malloc(sizeof(double) * (1 + c[0]->nic + c[0]->nec) * nfeval);
    // cost_fun_c(p, result, np, nfeval);
    cost_fun_c(p, result, np);

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

void C_call_grad(solnp_float *g, solnp_float *p, solnp_int np, solnp_int ngeval)
{

    SOLNP(timer)
    cost_timer;
    SOLNP(tic)
    (&cost_timer);
    solnp_int i;

    double *result = (double *)solnp_malloc(sizeof(double) * np * ngeval);
    grad_fun_c(p, result);

    for (i = 0; i < np * ngeval; i++)
    {
        g[i] = (solnp_float)result[i];
    }
    // cost_time += SOLNP(tocq)(&cost_timer) / 1e3;
}

void C_call_hess(solnp_float *h, solnp_float *p, solnp_int np, solnp_int nheval)
{

    SOLNP(timer)
    cost_timer;
    SOLNP(tic)
    (&cost_timer);
    solnp_int i;

    double *result = (double *)solnp_malloc(sizeof(double) * np * np * nheval);
    hess_fun_c(p, result);

    for (i = 0; i < np * np * nheval; i++)
    {
        h[i] = (solnp_float)result[i];
    }
    // cost_time += SOLNP(tocq)(&cost_timer) / 1e3;
}

void SOLNP_C_INTERFACE(
    // input
    SOLNPSettings *stgs,
    solnp_float *ibl,
    solnp_float *ibu,
    solnp_float *pbl,
    solnp_float *pbu,
    solnp_int *Ipc,
    solnp_int *Ipb,
    solnp_float *ib0,
    solnp_float *p,
    solnp_float *l,
    solnp_float *h,
    solnp_int np,
    solnp_int nic,
    solnp_int nec,
    // output
    SOLNPSol *sol,
    SOLNPInfo *info)
{

    SOLNP(timer)
    total_timer;
    SOLNP(tic)
    (&total_timer);

    solnp_int index, i;

    /*--------------------- check function ------------------------*/
    if (cost_fun_c == SOLNP_NULL)
    {
        solnp_printf("SOLNP+--> The user does not provided cost function in the fun structure. \n");
        solnp_printf("SOLNP+--> SOLNP stops.\n");
        exit(1);
    }

    if (grad_fun_c != SOLNP_NULL)
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

    if (hess_fun_c != SOLNP_NULL)
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

    memcpy(constraint->Ipc, Ipc, 2 * sizeof(solnp_int));
    memcpy(constraint->Ipb, Ipb, 2 * sizeof(solnp_int));

    /*--------------------- assemble input ------------------------*/
    solnp_int nc = nec + nic;
    solnp_int h_len = (np + nic) * (np + nic);
    SOLNPIput *input = (SOLNPIput *)solnp_malloc(sizeof(SOLNPIput));
    input->cnstr = constraint;
    input->stgs = (SOLNPSettings *)solnp_malloc(sizeof(SOLNPSettings));
    memcpy(input->stgs, stgs, sizeof(SOLNPSettings));
    // input->stgs = stgs;
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
    // first call cost function
    double *result = (double *)solnp_malloc(sizeof(double) * (nec + nic + 1));
    // cost_fun_c(p, result, np, 1);
    cost_fun_c(p, result, np);
    solnp_int m = nec + nic + 1;
    solnp_float *ob = (solnp_float *)solnp_malloc(sizeof(solnp_float) * m);
    for (i = 0; i < m; i++)
    {
        ob[i] = (solnp_float)result[i];
    }
    SOLNPCost *cost = malloc_cost(nec, nic, C_call_cost, C_call_grad, C_call_hess);

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

    solnp_float *ib0_p = (solnp_float *)solnp_malloc(sizeof(solnp_float) * (np + nic));

    if (nic > 0)
    {
        memcpy(ib0_p, ib0, nic * sizeof(solnp_float));
    }

    memcpy(&ib0_p[nic], p, np * sizeof(solnp_float));

    info->total_time = 0;
    /*--------------------- call solnp ------------------------*/
    solnp_int status = SOLNP(main)(input, cost, ib0_p, sol, info);
    /*--------------------- assemble output ------------------------*/
    info->total_time += SOLNP(tocq)(&total_timer) / 1e3;

    solnp_printf("SOLNP+--> SOLNP finished.\n");
}

solnp_int check_strict_in_bound(solnp_float *x, solnp_float *lb, solnp_float *ub, solnp_int len)
{
    /*
    check lb < x < ub
    */

    solnp_int inbound_status = 1;

    for (solnp_int i = 0; i < len; i++)
    {
        if (x[i] <= lb[i] || x[i] >= ub[i])
        {
            inbound_status = 0;
            break;
        }
    }

    return inbound_status;
}

void cal_avg_arr(solnp_float *x, solnp_float *y, solnp_float *x_y_2, solnp_int len)
{

    for (solnp_int i = 0; i < len; i++)
    {
        x_y_2[i] = (x[i] + y[i]) / 2;
    }
}

solnp_int check_var_bound(solnp_float *x, solnp_float *lb, solnp_float *ub, solnp_int len)
{
    /*
    check
        - x is not None
        - len(x) == len(lb) == len(ub)
        - lb < x < ub
    */

    solnp_int var_bound_status = 1;

    if (x == SOLNP_NULL)
    {
        var_bound_status = 0;
    }
    else
    {
        if (check_strict_in_bound(x, lb, ub, len) == 0)
        {
            var_bound_status = 0;
        }
    }
    return var_bound_status;
}

solnp_int check_bound(solnp_float *lb, solnp_float *ub, solnp_int len)
{
    /*
    check lb < ub
    */

    solnp_int bound_status = 1;

    for (solnp_int i = 0; i < len; i++)
    {
        if (lb[i] >= ub[i])
        {
            bound_status = 0;
            break;
        }
    }

    return bound_status;
}

void fill_array(solnp_float *arr, solnp_float value, solnp_int len)
{
    /*
    fill array with value
    */

    for (solnp_int i = 0; i < len; i++)
    {
        arr[i] = value;
    }
}

solnp_int check_prob(SOLNPProb *prob)
{
    /*
    check if prob is legal
    and init those None items
    */

    solnp_int prob_status = 1;

    // Ipc, Ipb
    if (prob->Ipc == SOLNP_NULL)
    {
        prob->Ipc = (solnp_int *)solnp_malloc(sizeof(solnp_int) * 2);
        prob->Ipc[0] = 0;
        prob->Ipc[1] = 0;
    }

    if (prob->Ipb == SOLNP_NULL)
    {
        prob->Ipb = (solnp_int *)solnp_malloc(sizeof(solnp_int) * 2);
        prob->Ipb[0] = 0;
        prob->Ipb[1] = 0;
    }

    if (prob->pbl)
    {
        prob->Ipc[0] = 1;
        prob->Ipb[0] = 1;
    }

    if (prob->pbu)
    {
        prob->Ipc[1] = 1;
        prob->Ipb[0] = 1;
    }

    if (prob->Ipb[0] + prob->nic > 0.5)
    {
        prob->Ipb[1] = 1;
    }

    // ibl, ibu and ib0
    if (prob->nic > 0)
    {
        // init ibl, ibu
        if (prob->ibl == SOLNP_NULL)
        {
            prob->ibl = (solnp_float *)solnp_malloc(sizeof(solnp_float) * prob->nic);
            fill_array(prob->ibl, -INFINITY, prob->nic);
        }
        if (prob->ibu == SOLNP_NULL)
        {
            prob->ibu = (solnp_float *)solnp_malloc(sizeof(solnp_float) * prob->nic);
            fill_array(prob->ibu, INFINITY, prob->nic);
        }
        // check and set ib0
        if (check_bound(prob->ibl, prob->ibu, prob->nic) == 0)
        {
            prob_status = 0;
            solnp_printf("Inequality bound error!");
            exit(1);
        }
        if (check_var_bound(prob->ib0, prob->ibl, prob->ibu, prob->nic) == 0)
        {
            if (prob->ib0 == SOLNP_NULL)
            {
                prob->ib0 = (solnp_float *)solnp_malloc(sizeof(solnp_float) * prob->nic);
            }
            cal_avg_arr(prob->ibl, prob->ibu, prob->ib0, prob->nic);
        }
    }
    // pbl, pbu and p0
    if (prob->np > 0)
    {
        // init pbl, pbu
        if (prob->pbl == SOLNP_NULL)
        {
            prob->pbl = (solnp_float *)solnp_malloc(sizeof(solnp_float) * prob->np);
            fill_array(prob->pbl, -INFINITY, prob->np);
        }
        if (prob->pbu == SOLNP_NULL)
        {
            prob->pbu = (solnp_float *)solnp_malloc(sizeof(solnp_float) * prob->np);
            fill_array(prob->pbu, INFINITY, prob->np);
        }
        // check and set p0
        if (check_bound(prob->pbl, prob->pbu, prob->np) == 0)
        {
            prob_status = 0;
            solnp_printf("Variable bound error!");
            exit(1);
        }
        if (check_var_bound(prob->p0, prob->pbl, prob->pbu, prob->np) == 0)
        {
            if (prob->p0 == SOLNP_NULL)
            {
                prob->p0 = (solnp_float *)solnp_malloc(sizeof(solnp_float) * prob->np);
            }
            cal_avg_arr(prob->pbl, prob->pbu, prob->p0, prob->np);
        }
    }
    else
    {
        prob_status = 0;
        solnp_printf("Dim of variables is 0!");
        exit(1);
    }

    return prob_status;
}

void SOLNP_PLUS(
    SOLNPProb *prob,
    SOLNPSettings *stgs,
    SOLNPSol *sol,
    SOLNPInfo *info,
    cost_temple *cost_fun,
    g_temple *grad_fun,
    h_temple *hess_fun,
    solnp_float *l,
    solnp_float *h)
{

    if (check_prob(prob) == 0)
    {
        solnp_printf("Illegal Problem!");
        exit(1);
    }

    solnp_int l_status = 1;

    if (l == SOLNP_NULL)
    {
        l_status = 0;
        l = (solnp_float *)solnp_malloc(sizeof(solnp_float) * prob->nc);
        fill_array(l, 0, prob->nc);
    }

    solnp_int h_status = 1;

    if (h == SOLNP_NULL)
    {
        h_status = 0;

        h = (solnp_float *)solnp_malloc(sizeof(solnp_float) * (prob->np + prob->nic) * (prob->np + prob->nic));
        fill_array(h, 0, (prob->np + prob->nic) * (prob->np + prob->nic));

        for (int diag = 0; diag < (prob->np + prob->nic); diag++)
        {
            h[diag * (prob->np + prob->nic) + diag] = 1;
        }
    }

    def_c_callback(cost_fun, grad_fun, hess_fun);

    SOLNP_C_INTERFACE(
        stgs,
        prob->ibl,
        prob->ibu,
        prob->pbl,
        prob->pbu,
        prob->Ipc,
        prob->Ipb,
        prob->ib0,
        prob->p0,
        l,
        h,
        prob->np,
        prob->nic,
        prob->nec,
        sol,
        info);

    if (l_status == 0)
    {
        solnp_free(l);
    }
    if (h_status == 0)
    {
        solnp_free(h);
    }
}
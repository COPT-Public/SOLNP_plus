#include "subnp.h"
#include "linalg.h"
#include "solnp_util.h"
#include "rescue.h"

SOLNPWork *init_work_RESCUE(
    SOLNPWork *w_old,
    SOLNPSettings *stgs)
{
    SOLNPWork *w = (SOLNPWork *)solnp_malloc(sizeof(SOLNPWork));
    w->n = w_old->n + w_old->nec;
    w->nec = w_old->nec;
    w->nic = w_old->nic;
    w->nc = w->nec + w->nic;
    w->ob = init_cost(w->nec, w->nic);
    copySOLNPCost(w->ob, w_old->best_fea_ob);
    w->ob->cost = w_old->ob->cost;
    w->pb = (SOLNPConstraint *)solnp_malloc(sizeof(SOLNPConstraint));

    w->pb->pu = (solnp_float *)solnp_malloc((w->n) * sizeof(solnp_float));
    memcpy(w->pb->pu, w_old->pb->pu, (w->n - w->nec) * sizeof(solnp_float));
    w->pb->pl = (solnp_float *)solnp_malloc((w->n) * sizeof(solnp_float));
    memcpy(w->pb->pl, w_old->pb->pl, (w->n - w->nec) * sizeof(solnp_float));
    for (solnp_int i = 0; i < w->nec; i++)
    {
        w->pb->pu[i + w->n - w->nec] = INFINITY;
        w->pb->pl[i + w->n - w->nec] = -INFINITY;
    }
    w->pb->n = w_old->pb->n;
    w->pb->nic = w_old->pb->nic;
    w->pb->Ipb = (solnp_int *)solnp_malloc(2 * sizeof(solnp_int));
    w->pb->Ipc = (solnp_int *)solnp_malloc(2 * sizeof(solnp_int));
    w->pb->il = (solnp_float *)solnp_malloc(w->nic * sizeof(solnp_float));
    w->pb->iu = (solnp_float *)solnp_malloc(w->nic * sizeof(solnp_float));
    memcpy(w->pb->il, w_old->pb->il, w->nic * sizeof(solnp_float));
    memcpy(w->pb->iu, w_old->pb->iu, w->nic * sizeof(solnp_float));
    memcpy(w->pb->Ipb, w_old->pb->Ipb, 2 * sizeof(solnp_int));
    memcpy(w->pb->Ipc, w_old->pb->Ipc, 2 * sizeof(solnp_int));

    w->rho = 1;

    w->constraint = (solnp_float *)solnp_malloc(w->nc * sizeof(solnp_float));
    w->p = (solnp_float *)solnp_malloc((w->n + w->nic) * sizeof(solnp_float));
    w->l = (solnp_float *)solnp_malloc(MAX(1, w->nc) * sizeof(solnp_float));
    w->bestl = (solnp_float *)solnp_malloc(MAX(1, w->nc) * sizeof(solnp_float));
    w->best_fea_l = (solnp_float *)solnp_malloc(MAX(1, w->nc) * sizeof(solnp_float));

    w->h = (solnp_float *)solnp_calloc((w->nic + w->n) * (w->nic + w->n), sizeof(solnp_float));
    w->jh = (solnp_float *)solnp_malloc((stgs->max_iter + 2) * sizeof(solnp_float));
    w->ch = (solnp_float *)solnp_malloc((stgs->max_iter + 2) * sizeof(solnp_float));
    w->bestp = (solnp_float *)solnp_malloc((w->n) * sizeof(solnp_float));
    w->best_fea_p = (solnp_float *)solnp_malloc((w->n) * sizeof(solnp_float));
    w->best_fea_ob = init_cost(w->nec, w->nic);
    w->bestob = init_cost(w->nec, w->nic);

    w->const_time = 0;
    w->count_cost = w_old->count_cost;
    w->count_grad = w_old->count_grad;
    w->count_hess = w_old->count_hess;
    w->exit = 0;

    memcpy(w->p, w_old->best_fea_p, (w_old->n) * sizeof(solnp_float));
    memcpy(w->best_fea_p, w_old->best_fea_p, (w_old->n) * sizeof(solnp_float));

    for (solnp_int i = 0; i < w->nec; i++)
    {
        w->p[i + w->n - w->nec] = 1;
    }
    memcpy(w->l, w_old->best_fea_l, w->nc * sizeof(solnp_float));
    SOLNP(add_scaled_array)
    (w->ob->ec, w->p + w->n + w->nic - w->nec, w->nec, -1.0);

    return w;
}

SOLNPWork *SOLNP_RESCUE_init(
    SOLNPWork *w_old,
    SOLNPSettings *stgs)
{
    SOLNPWork *w;
    w = init_work_RESCUE(w_old, stgs);

    return w;
}

solnp_int update_work_rescue(
    SOLNPWork *w,
    SOLNPSettings *stgs)
{
    solnp_int nec = w->nec;
    solnp_int nic = w->nic;
    solnp_int nc = w->nc;
    solnp_int i;

    w->bestobj = INFINITY;
    w->best_fea_con = INFINITY;
    if (nc > 0.5)
    {

        memcpy(w->constraint, w->ob->ec, nec * sizeof(solnp_float));
        memcpy(&w->constraint[nec], w->ob->ic, nic * sizeof(solnp_float));

        if (nic > 0.5)
        {
            for (i = 0; i < nic; i++)
            {
                if (w->constraint[nec + i] <= w->pb->il[i] || w->constraint[nec + i] >= w->pb->iu[i])
                    break;
            }
            if (i == nic)
                memcpy(w->p, &w->constraint[nec], nic * sizeof(solnp_float));

            SOLNP(add_scaled_array)
            (&w->constraint[nec], w->p, nic, -1.0);
        }

        w->cons_nm1 = SOLNP(norm)(w->constraint, w->nc);
        w->bestcon = w->cons_nm1;
    }
    else
    {
        w->l[0] = 0;
        w->constraint = SOLNP_NULL;
    }

    w->mu = w->n;
    w->best_fea_con = w->cons_nm1;
    w->j = SOLNP(norm_sq)(w->p + w->nic + w->n - w->nec, w->nec);
    w->ob->obj = w->j;
    w->jh[0] = w->j;
    // initiate w->h
    for (i = 0; i < w->nic + w->n; i++)
    {
        w->h[i + i * w->n] = 1;
    }

    return 0;
}

solnp_int SOLNP_RESCUE_solve(
    SOLNPWork *w,
    SOLNPSettings *stgs,
    SOLNPInfo *info)
{
    solnp_int i = 0;
    solnp_int j;

    solnp_int n = w->n;
    solnp_int nc = w->nc;
    solnp_int nec = w->nec;
    solnp_int nic = w->nic;
    solnp_int restart = 0;
    solnp_float delta0 = stgs->delta;
    update_work_rescue(w, stgs);

    for (i = 0; i < stgs->max_iter; i++)
    {
        subnp_qp(w, stgs, info);

        // w->ob->cost(w->ob, &w->p[w->nic], n, i);
        w->obj_dif = (w->j - w->ob->obj) / MAX(ABS(w->ob->obj), 1);
        w->j = w->ob->obj;

        if (nc > 0.5)
        {

            memcpy(w->constraint, w->ob->ec, nec * sizeof(solnp_float));
            memcpy(&w->constraint[nec], w->ob->ic, nic * sizeof(solnp_float));

            if (nic > 0.5)
            {

                for (j = 0; j < nic; j++)
                {
                    if (w->constraint[nec + j] <= w->pb->il[j] || w->constraint[nec + j] >= w->pb->iu[j])
                        break;
                }
                if (j == nic)
                    memcpy(w->p, &w->constraint[nec], nic * sizeof(solnp_float));

                SOLNP(add_scaled_array)
                (&w->constraint[nec], w->p, nic, -1.0);
            }

            w->cons_nm2 = SOLNP(norm)(w->constraint, nc);
            // previous 10
            if (w->cons_nm2 < 1 * stgs->tol_con)
            {
                w->rho = 0;
                w->mu = MIN(w->mu, stgs->tol);
            }
            // previously 5
            if (w->cons_nm2 < w->cons_nm1 && w->cons_nm2 < 10 * stgs->tol_con)
            {
                w->rho /= 5;
            }
            else if (w->cons_nm2 > 10 * w->cons_nm1)
            {
                w->rho = 5 * MAX(w->rho, SQRTF(stgs->tol));
            }
            if (w->exit == 0 && MAX(w->obj_dif, w->cons_nm1 - w->cons_nm2) <= 0 && restart < stgs->re_time && (stgs->delta <= MAX(3 * stgs->tol, stgs->delta_end) || stgs->grad))
            {
                SOLNP(scale_array)
                (w->l, 0, nc);
                restart++;
                if (stgs->noise)
                {
                    stgs->delta = delta0;
                }
                // h=diag(diag(h));
                if (stgs->noise)
                {
                    stgs->delta = delta0;
                }
                for (int col = 0; col < w->nic + w->n; col++)
                {
                    for (int row = 0; row < w->nic + w->n; row++)
                    {
                        if (col != row)
                        {
                            w->h[col * (w->nic + w->n) + row] = 0;
                        }
                    }
                }
            }

            w->cons_nm1 = w->cons_nm2;
        }

        w->jh[1 + i] = w->j;
        w->jh[2 + i] = INFINITY;
        w->ch[1 + i] = w->cons_nm1;

        if ((w->cons_nm1 <= stgs->tol_con && w->obj_dif <= stgs->tol && (stgs->delta <= MAX(stgs->tol, stgs->delta_end) || stgs->grad)) || w->exit)
        {
            if (restart < stgs->re_time && w->alm_crit > stgs->tol_restart)
            {
                restart++;
                if (stgs->noise)
                {
                    stgs->delta = delta0;
                }
                // h=diag(diag(h));
                for (int col = 0; col < w->nic + w->n; col++)
                {
                    for (int row = 0; row < w->nic + w->n; row++)
                    {
                        if (col != row)
                        {
                            w->h[col * (w->nic + w->n) + row] = 0;
                        }
                    }
                }
            }
            else
            {
                i++;
                break;
            }
        }
    }

    if (w->cons_nm1 <= stgs->tol_con && w->obj_dif <= stgs->tol && (stgs->delta <= MAX(stgs->tol, stgs->delta_end) || stgs->grad == 1))
    {
        printf("SOLNP+--> RESCUE process Success! Completed in %d iterations\n", i);
    }
    else
    {
        printf("SOLNP+--> Exiting after maximum number of function evaluation. Rescue Process Fails.\n");
    }
    return i;
}

SOLNPWork *SOLNP_RESCUE(
    SOLNPWork *w_old,
    SOLNPSettings *stgs,
    SOLNPInfo *info // record running time
)
{
    solnp_int status;
    SOLNPWork *w = SOLNP_RESCUE_init(w_old, stgs);
    if (w)
    {
        SOLNP_RESCUE_solve(w, stgs, info);
    }

    return w;
}

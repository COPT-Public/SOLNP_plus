#include "qp_solver.h"
#include "der_info.h"

void proj(
    solnp_float *p,
    SOLNPWork *w)
{
    solnp_int i;
    for (i = 0; i < w->nic; i++)
    {
        if (p[i] > w->pb->iu[i])
        {
            p[i] = w->pb->iu[i];
        }
        if (p[i] < w->pb->il[i])
        {
            p[i] = w->pb->il[i];
        }
    }
    for (i = 0; i < w->n; i++)
    {
        if (p[w->nic + i] > w->pb->pu[i])
        {
            p[w->nic + i] = w->pb->pu[i];
        }
        if (p[w->nic + i] < w->pb->pl[i])
        {
            p[w->nic + i] = w->pb->pl[i];
        }
    }
}

// solnp_int qpsolver(
solnp_float qpsolver(
    solnp_int n,
    solnp_float *H,
    solnp_float *g,
    solnp_int m1,
    solnp_int m2,
    solnp_float *A,
    solnp_float *b,
    solnp_float *lb,
    solnp_float *ub,
    solnp_float *p,
    solnp_float *lambda)
{
    /*
    This function solve the following QP(LP) problem:
           0.5 * x^T H x + g^T x
           subject to Ax = b,
               lb <= x <= ub,
    where H is n *n SPD matrix, A is m1 * m2. p is primal solution, lambda is the multiplier
    with repect to equality constraints.
    If lambda = SOLNP_NULL, it means we don't need a multiplier. If n = 0, it becomes a LP problem.
    */
    c_int m = m1;
    c_int n1 = n;
    c_int i;
    c_int P_nnz, A_nnz;
    c_float *P_x;
    c_int *P_i;
    c_int *P_p;

    SOLNP(timer)
    qpsolver_timer;
    SOLNP(tic)
    (&qpsolver_timer);

    /*Setup two matrice*/
    if (H == SOLNP_NULL)
    {
        P_nnz = 0;
        P_x = OSQP_NULL;
        P_i = OSQP_NULL;
        P_p = c_malloc((n + 1) * sizeof(c_int));
        for (i = 0; i < n + 1; i++)
        {
            P_p[i] = 0;
        }
    }
    else
    {
        P_nnz = countA_sys(n, n, H);
        P_p = (c_int *)c_malloc((n + 1) * sizeof(c_int));
        P_i = (c_int *)c_malloc(P_nnz * sizeof(c_int));
        P_x = (c_float *)c_malloc(P_nnz * sizeof(c_float));
        calculate_csc_sys(n, n, H, P_x, P_i, P_p);
    }
    A_nnz = countA(m, n, A);
    c_int *A_p = (c_int *)c_malloc((n + 1) * sizeof(c_int));
    c_int *A_i = (c_int *)c_malloc(A_nnz * sizeof(c_int));
    c_float *A_x = (c_float *)c_malloc(A_nnz * sizeof(c_float));
    calculate_csc(m, n, A, A_x, A_i, A_p);

    c_float *q = (c_float *)c_malloc(n * sizeof(c_float));
    memcpy(q, g, n * sizeof(c_float));

    c_float *l = (c_float *)c_malloc((m + n) * sizeof(c_float));
    memcpy(l, b, m * sizeof(c_float));
    memcpy(&l[m], lb, n * sizeof(c_float));

    c_float *u = (c_float *)c_malloc((m + n) * sizeof(c_float));
    memcpy(u, b, m * sizeof(c_float));
    memcpy(&u[m], ub, n * sizeof(c_float));

    // Exitflag
    c_int exitflag = 0;

    // Workspace structures
    OSQPWorkspace *work;
    OSQPSettings *settings = (OSQPSettings *)c_malloc(sizeof(OSQPSettings));
    OSQPData *data = (OSQPData *)c_malloc(sizeof(OSQPData));

    // Populate data
    if (data)
    {
        data->n = n1;
        data->m = m + n1;
        data->P = csc_matrix(data->n, data->n, P_nnz, P_x, P_i, P_p);
        data->q = q;
        data->A = csc_matrix(data->m, data->n, A_nnz, A_x, A_i, A_p);
        data->l = l;
        data->u = u;
    }

    // Define solver settings as default
    if (settings)
        osqp_set_default_settings(settings);
    settings->eps_abs = 1e-4;
    settings->eps_rel = 1e-4;
    // Setup workspace
    exitflag = osqp_setup(&work, data, settings);

    // Solve Problem
    osqp_solve(work);

    memcpy(p, work->solution->x, n * sizeof(c_float));
    if (lambda)
    {
        memcpy(lambda, work->solution->y, m1 * sizeof(c_float));
        SOLNP(scale_array)
        (lambda, -1.0, m1);
    }

    // Clean workspace
    osqp_cleanup(work);
    if (data)
    {
        if (data->A)
            c_free(data->A);
        if (data->P)
            c_free(data->P);
        c_free(q);
        c_free(l);
        c_free(u);
        c_free(data);
    }
    c_free(A_x);
    c_free(A_i);
    c_free(A_p);

    if (P_x)
    {
        c_free(P_x);
    }
    if (P_i)
    {
        c_free(P_i);
    }
    c_free(P_p);
    if (settings)
        c_free(settings);

    // qpsolver_time += SOLNP(tocq)(&qpsolver_timer) / 1e3;
    // return exitflag;

    solnp_float qpsolver_time = SOLNP(tocq)(&qpsolver_timer) / 1e3;
    return qpsolver_time;
}
solnp_int find_int_feas_sol_osqp(
    SUBNPWork *w_sub,
    SOLNPWork *w,
    SOLNPSettings *stgs,
    SOLNPInfo *info)
{
    w_sub->ch = -1;
    w_sub->alp[0] = stgs->tol - SOLNP(norm_inf)(&w_sub->J->a[w_sub->npic * w_sub->J->a_row], w_sub->J->a_row);

    if (w_sub->alp[0] <= 0)
    {
        w_sub->ch = 1;
        if (w->pb->Ipb[1] == 0)
        {
            // TODO: project onto hyperplane
            // need linear system solver
            // p0=p0-a'*((a*a')\constraint);

            // p0: w->p is a w_sub->npic vector
            // a:  w_sub->J->a is a w_sub->J->a_row * w_sub->npic matrix
            //     in col major form
            // constraint: &w_sub->J->a[w_sub->J->a_row * w_sub->npic] is a w_sub->J->a_row vector

            solnp_float *aT = (solnp_float *)solnp_malloc(w_sub->J->a_row * w_sub->npic * sizeof(solnp_float));
            SOLNP(transpose)
            (w_sub->J->a_row, w_sub->npic, w_sub->J->a, aT);
            solnp_float *aaT = (solnp_float *)solnp_malloc(w_sub->J->a_row * w_sub->J->a_row * sizeof(solnp_float));
            SOLNP(AB)
            (aaT, w_sub->J->a, aT, w_sub->J->a_row, w_sub->npic, w_sub->J->a_row);
            SOLNP(chol)
            (w_sub->J->a_row, aaT);
            solnp_float *inv_aaT_cons = (solnp_float *)solnp_malloc(w_sub->J->a_row * sizeof(solnp_float));
            memcpy(inv_aaT_cons, &w_sub->J->a[w_sub->J->a_row * w_sub->npic], w_sub->J->a_row * sizeof(solnp_float));
            SOLNP(solve_lin_sys)
            (w_sub->J->a_row, 1, aaT, inv_aaT_cons);
            solnp_float *aT_inv_aaT_cons = (solnp_float *)solnp_malloc(w_sub->npic * sizeof(solnp_float));
            SOLNP(Ax)
            (aT_inv_aaT_cons, aT, inv_aaT_cons, w_sub->npic, w_sub->J->a_row);
            SOLNP(add_scaled_array)
            (w->p, aT_inv_aaT_cons, w_sub->npic, -1.0);

            solnp_free(aaT);
            solnp_free(aT);
            solnp_free(inv_aaT_cons);
            solnp_free(aT_inv_aaT_cons);

            w_sub->alp[0] = 1;
        }
        else
        {
            // TODO: affine scaling method to solve LP
            solnp_float *p = (solnp_float *)solnp_calloc((w_sub->npic + 1), sizeof(solnp_float));
            solnp_float *c = (solnp_float *)solnp_calloc((w_sub->npic + 1), sizeof(solnp_float));
            solnp_float *y = (solnp_float *)solnp_calloc((w_sub->J->a_row), sizeof(solnp_float));
            solnp_float *lb = (solnp_float *)solnp_calloc((w_sub->npic + 1), sizeof(solnp_float));
            solnp_float *ub = (solnp_float *)solnp_calloc((w_sub->npic + 1), sizeof(solnp_float));
            // init p, c, dx

            c[w_sub->npic] = 1.;
            lb[w_sub->npic] = 0;
            ub[w_sub->npic] = INFINITY;
            memcpy(lb, w->pb->il, w->pb->nic * sizeof(solnp_float));
            memcpy(ub, w->pb->iu, w->pb->nic * sizeof(solnp_float));
            memcpy(&lb[w->pb->nic], w->pb->pl, w->pb->n * sizeof(solnp_float));
            memcpy(&ub[w->pb->nic], w->pb->pu, w->pb->n * sizeof(solnp_float));

            // last column of a should be **negative** constraint
            SOLNP(scale_array)
            (&w_sub->J->a[w_sub->J->a_row * w_sub->npic], -1.0, w_sub->J->a_row);

            info->qpsolver_time += qpsolver(w_sub->npic + 1, SOLNP_NULL, c, w->nc, w_sub->npic + 1, w_sub->J->a, w_sub->b, lb, ub, p, SOLNP_NULL);

            proj(p, w);
            memcpy(w->p, p, w_sub->npic * sizeof(solnp_float));
            SOLNP(Ax)
            (w_sub->b, w_sub->J->a, w->p, w_sub->J->a_row, w_sub->npic);

            // free pointers
            solnp_free(p);
            solnp_free(c);
            solnp_free(y);
            solnp_free(lb);
            solnp_free(ub);
        }
    }

    return 0;
}

// solnp_int find_int_feas_sol_aff(
solnp_float find_int_feas_sol_aff(
    SUBNPWork *w_sub,
    SOLNPWork *w,
    SOLNPSettings *stgs,
    SOLNPInfo *info)
{

    SOLNP(timer)
    qpsolver_timer;
    SOLNP(tic)
    (&qpsolver_timer);
    w_sub->ch = -1;
    w_sub->alp[0] = stgs->tol - SOLNP(norm_inf)(&w_sub->J->a[w_sub->npic * w_sub->J->a_row], w_sub->J->a_row);

    if (w_sub->alp[0] <= 0)
    {
        w_sub->ch = 1;
        if (w->pb->Ipb[1] == 0)
        {
            // TODO: project onto hyperplane
            // need linear system solver
            // p0=p0-a'*((a*a')\constraint);

            // p0: w->p is a w_sub->npic vector
            // a:  w_sub->J->a is a w_sub->J->a_row * w_sub->npic matrix
            //     in col major form
            // constraint: &w_sub->J->a[w_sub->J->a_row * w_sub->npic] is a w_sub->J->a_row vector

            /*solnp_float* aT = (solnp_float*)solnp_malloc(w_sub->J->a_row * w_sub->npic * sizeof(solnp_float));
            SOLNP(transpose)(w_sub->J->a_row, w_sub->npic, w_sub->J->a, aT);
            solnp_float* aaT = (solnp_float*)solnp_malloc(w_sub->J->a_row * w_sub->J->a_row * sizeof(solnp_float));
            SOLNP(AB)(aaT, w_sub->J->a, aT, w_sub->J->a_row, w_sub->npic, w_sub->J->a_row);
            SOLNP(chol)(w_sub->J->a_row, aaT);
            solnp_float* inv_aaT_cons = (solnp_float*)solnp_malloc(w_sub->J->a_row * sizeof(solnp_float));
            memcpy(inv_aaT_cons, &w_sub->J->a[w_sub->J->a_row * w_sub->npic], w_sub->J->a_row * sizeof(solnp_float));
            SOLNP(solve_lin_sys)(w_sub->J->a_row, 1, aaT, inv_aaT_cons);
            solnp_float* aT_inv_aaT_cons = (solnp_float*)solnp_malloc(w_sub->npic * sizeof(solnp_float));
            SOLNP(Ax)(aT_inv_aaT_cons, aT, inv_aaT_cons, w_sub->npic, w_sub->J->a_row);*/

            // Use SVD to solve the Least square problem

            solnp_float *aT_inv_aaT_cons = (solnp_float *)solnp_malloc(w_sub->npic * sizeof(solnp_float));
            memcpy(aT_inv_aaT_cons, &w_sub->J->a[w_sub->J->a_row * w_sub->npic], w_sub->J->a_row * sizeof(solnp_float));
            solnp_float *a_temp = (solnp_float *)solnp_malloc(w_sub->J->a_row * w_sub->npic * sizeof(solnp_float));
            memcpy(a_temp, w_sub->J->a, w_sub->J->a_row * w_sub->npic * sizeof(solnp_float));
            solnp_int info = SOLNP(least_square)(w_sub->J->a_row, w_sub->npic, a_temp, aT_inv_aaT_cons);

            SOLNP(add_scaled_array)
            (w->p, aT_inv_aaT_cons, w_sub->npic, -1.0);

            /*solnp_free(aaT);
            solnp_free(aT);
            solnp_free(inv_aaT_cons);*/
            solnp_free(a_temp);
            solnp_free(aT_inv_aaT_cons);

            w_sub->alp[0] = 1;
        }
        else
        {
            // TODO: affine scaling method to solve LP
            solnp_float *p = (solnp_float *)solnp_calloc((w_sub->npic + 1), sizeof(solnp_float));
            solnp_float *c = (solnp_float *)solnp_calloc((w_sub->npic + 1), sizeof(solnp_float));
            solnp_float *dx = (solnp_float *)solnp_calloc((w_sub->npic + 1), sizeof(solnp_float));
            solnp_float *y = (solnp_float *)solnp_calloc((w_sub->J->a_row), sizeof(solnp_float));
            solnp_float *v = (solnp_float *)solnp_calloc((w_sub->npic + 1), sizeof(solnp_float));
            solnp_int i, j, minit = 0;
            solnp_float go = 1.0, z;

            // init p, c, dx
            memcpy(p, w->p, w_sub->npic * sizeof(solnp_float));
            p[w_sub->npic] = 1.;
            c[w_sub->npic] = 1.;

            // last column of a should be **negative** constraint
            SOLNP(scale_array)
            (&w_sub->J->a[w_sub->J->a_row * w_sub->npic], -1.0, w_sub->J->a_row);

            for (i = 0; i <= w_sub->npic; i++)
            {
                dx[i] = 1.;
            }
            // POTENTIAL PROBLEM: the last column should be -constraint? -Right. This bug has been fixed.

            // affine scaling while loop
            while (go >= stgs->tol)
            {
                minit += 1;

                // form ellipse
                for (i = 0; i < w->nic; i++)
                {
                    dx[i] = MIN(MIN(p[i] - w->pb->il[i], w->pb->iu[i] - p[i]), 1e2);
                }
                dx[w_sub->npic] = p[w_sub->npic];
                if (w->pb->Ipb[0] == 0)
                { // x is free
                    for (i = 0; i < w->n; i++)
                    {
                        solnp_float dx_norm_inf = SOLNP(norm_inf)(dx, w->nic);
                        dx_norm_inf = MAX(dx_norm_inf, 100);
                        dx[w->nic + i] = dx_norm_inf;
                    }
                }
                else
                { // x is not free
                    for (i = 0; i < w->n; i++)
                    {
                        dx[w->nic + i] = MIN(MIN(p[w->nic + i] - w->pb->pl[i], w->pb->pu[i] - p[w->nic + i]), 1e1);
                    }
                }

                // TODO: solve for direction v
                // y=(a*diag(dx))'\(dx.*c');
                // v=dx.*(dx.*(c'-a'*y));

                // a: w_sub->J->a is a w_sub->J->a_row * (w_sub->npic+1) matrix
                //    in col major form
                // c': w_sub->npic+1 col vector (c is row vector in matlab)
                // dx: w_sub->npic+1 col vector
                // y: w_sub->J->a_row col vector
                // v: w_sub->npic+1 col vector

                solnp_float *a_diag_dx = (solnp_float *)solnp_malloc(w_sub->J->a_row * (w_sub->npic + 1) * sizeof(solnp_float));
                memcpy(a_diag_dx, w_sub->J->a, w_sub->J->a_row * (w_sub->npic + 1) * sizeof(solnp_float));

                for (i = 0; i < w_sub->npic + 1; i++)
                {
                    for (j = 0; j < w_sub->J->a_row; j++)
                    {
                        a_diag_dx[i * w_sub->J->a_row + j] *= dx[i];
                    }
                }

                solnp_float *a_diag_dx_T = (solnp_float *)solnp_malloc(w_sub->J->a_row * (w_sub->npic + 1) * sizeof(solnp_float));
                SOLNP(transpose)
                (w_sub->J->a_row, w_sub->npic + 1, a_diag_dx, a_diag_dx_T);

                solnp_float *a_diag_dx2_aT = (solnp_float *)solnp_malloc(w_sub->J->a_row * w_sub->J->a_row * sizeof(solnp_float));
                SOLNP(AB)
                (a_diag_dx2_aT, a_diag_dx, a_diag_dx_T, w_sub->J->a_row, w_sub->npic + 1, w_sub->J->a_row);

                solnp_float *dx_c = (solnp_float *)solnp_malloc((w_sub->npic + 1) * sizeof(solnp_float));
                memcpy(dx_c, c, (w_sub->npic + 1) * sizeof(solnp_float));
                for (i = 0; i < w_sub->npic + 1; i++)
                {
                    dx_c[i] *= dx[i];
                }

                SOLNP(Ax)
                (y, a_diag_dx, dx_c, w_sub->J->a_row, w_sub->npic + 1);
                SOLNP(chol)
                (w_sub->J->a_row, a_diag_dx2_aT);
                SOLNP(solve_lin_sys)
                (w_sub->J->a_row, 1, a_diag_dx2_aT, y);
                solnp_float *aT = (solnp_float *)solnp_malloc(w_sub->J->a_row * (w_sub->npic + 1) * sizeof(solnp_float));
                SOLNP(transpose)
                (w_sub->J->a_row, w_sub->npic + 1, w_sub->J->a, aT);
                SOLNP(Ax)
                (v, aT, y, w_sub->npic + 1, w_sub->J->a_row);
                SOLNP(add_scaled_array)
                (v, c, w_sub->npic + 1, -1.0);
                for (i = 0; i < w_sub->npic + 1; i++)
                {
                    v[i] *= -dx[i] * dx[i];
                }

                solnp_free(a_diag_dx);
                solnp_free(a_diag_dx_T);
                solnp_free(a_diag_dx2_aT);
                solnp_free(dx_c);
                solnp_free(aT);

                if (v[w_sub->npic] > 0)
                {
                    // calculate step size
                    z = p[w_sub->npic] / v[w_sub->npic];
                    for (i = 0; i < w->nic; i++)
                    {
                        if (v[i] < 0)
                        {
                            z = MIN(z, -(w->pb->iu[i] - p[i]) / v[i]);
                        }
                        else if (v[i] > 0)
                        {
                            z = MIN(z, (p[i] - w->pb->il[i]) / v[i]);
                        }
                    }
                    if (w->pb->Ipb[0] > 0)
                    {
                        for (i = 0; i < w->n; i++)
                        {
                            if (v[w->nic + i] < 0)
                            {
                                z = MIN(z, -(w->pb->pu[i] - p[w->nic + i]) / v[w->nic + i]);
                            }
                            else if (v[w->nic + i] > 0)
                            {
                                z = MIN(z, (p[w->nic + i] - w->pb->pl[i]) / v[w->nic + i]);
                            }
                        }
                    }

                    // update
                    if (z >= (p[w_sub->npic] / v[w_sub->npic]))
                    {
                        SOLNP(add_scaled_array)
                        (p, v, w_sub->npic + 1, -z);
                    }
                    else
                    {
                        SOLNP(add_scaled_array)
                        (p, v, w_sub->npic + 1, -0.9 * z);
                    }
                    go = p[w_sub->npic];
                    if (minit >= 10)
                    {
                        go = 0;
                    }
                }
                else
                {
                    go = 0;
                    minit = 10;
                }
            }

            if (minit >= 10)
            {
                solnp_printf("SOLNP+--> ");
                solnp_printf("The linearized problem has no feasible     \n");
                solnp_printf("         ");
                solnp_printf("solution.  The problem may not be feasible.\n");
            }

            memcpy(w->p, p, w_sub->npic * sizeof(solnp_float));
            SOLNP(Ax)
            (w_sub->b, w_sub->J->a, w->p, w_sub->J->a_row, w_sub->npic);

            // free pointers
            solnp_free(p);
            solnp_free(c);
            solnp_free(dx);
            solnp_free(y);
            solnp_free(v);
        }
    }

    // qpsolver_time += SOLNP(tocq)(&qpsolver_timer) / 1e3;
    // return 0;

    solnp_float qpsolver_time = SOLNP(tocq)(&qpsolver_timer) / 1e3;
    info->qpsolver_time += qpsolver_time;

    return qpsolver_time;
}

void solve_qp_osqp(
    SOLNPWork *w,
    SUBNPWork *w_sub,
    solnp_float *p0,
    solnp_float *y,
    solnp_float *g,
    SOLNPInfo *info)
{
    solnp_int i;

    solnp_float *u = (solnp_float *)solnp_malloc(w_sub->J->npic * sizeof(solnp_float));
    solnp_float *c = (solnp_float *)solnp_calloc(w_sub->J->npic, sizeof(solnp_float));

    solnp_float *lb = (solnp_float *)solnp_malloc(w_sub->J->npic * sizeof(solnp_float));
    solnp_float *ub = (solnp_float *)solnp_malloc(w_sub->J->npic * sizeof(solnp_float));
    for (i = 0; i < w_sub->J->npic; i++)
    {
        lb[i] = -100. / SQRTF(w_sub->J->npic);
        ub[i] = 100. / SQRTF(w_sub->J->npic);
    }
    for (i = 0; i < w_sub->J->npic; i++)
    {
        if (i < w->nic)
        {
            lb[i] = MAX(lb[i], w->pb->il[i] - w->p[i]);
            ub[i] = MIN(ub[i], w->pb->iu[i] - w->p[i]);
        }
        else
        {
            lb[i] = MAX(lb[i], w->pb->pl[i - w->nic] - w->p[i]);
            ub[i] = MIN(ub[i], w->pb->pu[i - w->nic] - w->p[i]);
        }
    }
    info->qpsolver_time += qpsolver(w_sub->J->npic, w->h, g, w_sub->J->a_row, w_sub->J->npic, w_sub->J->a, c, lb, ub, u, y);
    solnp_free(lb);
    solnp_free(ub);
    /*
    else {
        solnp_float* dx = (solnp_float*)solnp_malloc(w_sub->J->npic * sizeof(solnp_float));
        solnp_float* L = (solnp_float*)solnp_malloc(w_sub->J->npic * w_sub->J->npic * sizeof(solnp_float));
        for (i = 0; i < w_sub->J->npic; i++) {
            dx[i] = 0.01;
        }
        memcpy(L, w->h, (w->nic + w->n) * (w->nic + w->n) * sizeof(solnp_float));
        // perform cholesky to h+mu* diag(dx.^2)
        for (i = 0; i < (w->nic + w->n); i++) {
            L[i * (w->nic + w->n) + i] += w->mu * dx[i] * dx[i];
        }
        SOLNP(chol)((w->nic + w->n), L);

        // u = -c* yg;
        memcpy(u, g, w_sub->J->npic * sizeof(solnp_float));
        SOLNP(scale_array)(u, -1.0, w_sub->J->npic);
        SOLNP(solve_lin_sys)(w_sub->J->npic, 1, L, u);
        solnp_free(L);
        solnp_free(dx);
    }*/
    SOLNP(add_scaled_array)
    (u, w->p, w_sub->J->npic, 1.0);
    memcpy(p0, u, w_sub->J->npic * sizeof(solnp_float));
    proj(p0, w);
    solnp_free(c);
    solnp_free(u);

    return;
}
// void solve_qp_aff(
solnp_float solve_qp_aff(
    SOLNPWork *w,
    SUBNPWork *w_sub,
    solnp_float *p0,
    solnp_float *y,
    solnp_float *g,
    solnp_int sys,
    SOLNPInfo *info)
{

    SOLNP(timer)
    qpsolver_timer;
    SOLNP(tic)
    (&qpsolver_timer);

    solnp_float go = -1.0;
    solnp_int i, k;
    solnp_float *dx = (solnp_float *)solnp_malloc(w_sub->J->npic * sizeof(solnp_float));
    // calculate the ball
    for (i = 0; i < w_sub->J->npic; i++)
    {
        dx[i] = 0.01;
    }
    if (w->pb->Ipb[1] > 0.5)
    {
        solnp_float m = 0.01; // minimum of dx[i]
        for (i = 0; i < w_sub->mm; i++)
        {
            if (i < w->nic)
            {
                dx[i] = MAX(1 / (MIN(w->p[i] - w->pb->il[i], w->pb->iu[i] - w->p[i]) + EPS), 0.1);
            }
            else
            {
                dx[i] = MAX(1 / (MIN(w->p[i] - w->pb->pl[i - w->nic], w->pb->pu[i - w->nic] - w->p[i]) + EPS), 0.1);
            }

            if (dx[i] < m)
            {
                m = dx[i];
            }
        }
        if (w->pb->Ipb[0] <= 0)
        {
            m = MIN(m, 0.01);
            for (i = w_sub->mm; i < w_sub->J->npic; i++)
            {
                dx[i] = m;
            }
        }
    }
    w->mu /= 10; // mu is l in subnp.m, lagrange multiplier
    // solnp_float* c = (solnp_float *)solnp_malloc(w_sub->J->npic * w_sub->J->npic * sizeof(solnp_float));
    solnp_float *u = (solnp_float *)solnp_calloc((w_sub->J->npic + w_sub->J->a_row), sizeof(solnp_float));
    solnp_int affit = 0;
    while (go <= 0 && affit < 50)
    {
        // SOLNP(chol)(w_sub->J->npic,w->h,dx,c,w->mu);
        // c = inv(c);
        // yg = c' *g
        affit++;
        if (sys == 1)
        {
            solnp_float *L = (solnp_float *)solnp_malloc(w_sub->J->npic * w_sub->J->npic * sizeof(solnp_float));
            memcpy(L, w->h, (w->nic + w->n) * (w->nic + w->n) * sizeof(solnp_float));
            // perform cholesky to h+mu* diag(dx.^2)
            for (i = 0; i < (w->nic + w->n); i++)
            {
                L[i * (w->nic + w->n) + i] += w->mu * dx[i] * dx[i];
            }
            SOLNP(chol)
            ((w->nic + w->n), L);

            // u = -c* yg;
            memcpy(u, g, w_sub->J->npic * sizeof(solnp_float));
            SOLNP(scale_array)
            (u, -1.0, w_sub->J->npic);
            SOLNP(solve_lin_sys)
            (w_sub->J->npic, 1, L, u);
            if (w_sub->J->a_row > 0)
            {
                solnp_float *inv_H_aT = (solnp_float *)solnp_malloc(w_sub->J->npic * w_sub->J->a_row * sizeof(solnp_float));
                SOLNP(transpose)
                (w_sub->J->a_row, w_sub->J->npic, w_sub->J->a, inv_H_aT);
                SOLNP(solve_lin_sys)
                (w_sub->J->npic, w_sub->J->a_row, L, inv_H_aT);

                solnp_float *a_inv_H_aT = (solnp_float *)solnp_malloc((w_sub->J->a_row * w_sub->J->a_row) * sizeof(solnp_float));
                SOLNP(AB)
                (a_inv_H_aT, w_sub->J->a, inv_H_aT, w_sub->J->a_row, w_sub->J->npic, w_sub->J->a_row);

                SOLNP(chol)
                (w_sub->J->a_row, a_inv_H_aT);
                SOLNP(Ax)
                (y, w_sub->J->a, u, w_sub->J->a_row, w_sub->J->npic); // au = a  H^-1 u
                SOLNP(scale_array)
                (y, -1.0, w_sub->J->a_row);

                SOLNP(solve_lin_sys)
                (w_sub->J->a_row, 1, a_inv_H_aT, y); // a_inv_H_aT * y = a*H^-1*u

                solnp_float *inv_H_aT_y = (solnp_float *)solnp_malloc(w_sub->J->npic * sizeof(solnp_float));
                SOLNP(Ax)
                (inv_H_aT_y, inv_H_aT, y, w_sub->J->npic, w_sub->J->a_row);
                SOLNP(add_scaled_array)
                (u, inv_H_aT_y, w_sub->J->npic, 1.0);

                solnp_free(inv_H_aT);
                solnp_free(a_inv_H_aT);
                solnp_free(inv_H_aT_y);
            }
            solnp_free(L);
        }
        if (sys == 2)
        {
            // directly solve the linear system
            memcpy(u, w_sub->J->g, w_sub->npic * sizeof(solnp_float));
            SOLNP(set_as_scaled_array)
            (u, u, -1., w_sub->npic);
            for (i = 0; i < w_sub->J->a_row; i++)
            {
                u[i + w_sub->npic] = 0;
            }

            solnp_float *kkt_sys = (solnp_float *)solnp_calloc((w_sub->J->npic + w_sub->J->a_row) * (w_sub->J->a_row + w_sub->J->npic), sizeof(solnp_float));
            for (i = 0; i < w_sub->npic; i++)
            {
                memcpy(kkt_sys + i * (w_sub->J->npic + w_sub->J->a_row), w->h + i * w_sub->J->npic, w_sub->J->npic * sizeof(solnp_float));
                memcpy(kkt_sys + i * (w_sub->J->npic + w_sub->J->a_row) + w_sub->J->npic, w_sub->J->a + i * w_sub->J->a_row, w_sub->J->a_row * sizeof(solnp_float));
                kkt_sys[i + i * (w_sub->J->npic + w_sub->J->a_row)] += w->mu * dx[i] * dx[i];
            }
            for (i = 0; i < w_sub->J->a_row; i++)
            {
                for (k = 0; k < w_sub->npic; k++)
                {
                    kkt_sys[(i + w_sub->npic) * (w_sub->J->npic + w_sub->J->a_row) + k] = w_sub->J->a[i + k * w_sub->J->a_row];
                }
            }

            solnp_int sol_state = SOLNP(solve_sys_lin_sys)(w_sub->J->npic + w_sub->J->a_row, kkt_sys, u, w_sub->J->a_row);
            if (isnan(u[0]))
            {
                w->exit = 2;
                break;
            }
            if (sol_state != 0)
            {
                w->mu = 3 * w->mu;
                solnp_free(kkt_sys);
                continue;
            }

            memcpy(y, u + w_sub->npic, w_sub->J->a_row * sizeof(solnp_float));
            solnp_free(kkt_sys);
        }

        SOLNP(add_scaled_array)
        (u, w->p, w_sub->J->npic, 1.0);
        memcpy(p0, u, w_sub->J->npic * sizeof(solnp_float));
        if (isnan(p0[0]))
        {
            w->exit = 2;
            break;
        }

        if (w->pb->Ipb[1] <= 0)
        {
            go = 1;
        }
        else
        {
            go = INFINITY;
            for (i = 0; i < w_sub->mm; i++)
            {
                if (i < w->nic && go > MIN(p0[i] - w->pb->il[i], w->pb->iu[i] - p0[i]))
                {
                    go = MIN(p0[i] - w->pb->il[i], w->pb->iu[i] - p0[i]);
                }
                else if (i >= w->nic && go > MIN(p0[i] - w->pb->pl[i - w->nic], w->pb->pu[i - w->nic] - p0[i]))
                {
                    go = MIN(p0[i] - w->pb->pl[i - w->nic], w->pb->pu[i - w->nic] - p0[i]);
                }
            }
            w->mu = 3 * w->mu;
        }
    }
    solnp_free(u);
    solnp_free(dx);

    // qpsolver_time += SOLNP(tocq)(&qpsolver_timer) / 1e3;
    // return;

    solnp_float qpsolver_time = SOLNP(tocq)(&qpsolver_timer) / 1e3;
    info->qpsolver_time += qpsolver_time;
    return qpsolver_time;
}

void Gradient_descent(
    SOLNPWork *w,
    solnp_float *p0,
    solnp_float *g,
    solnp_float stepsize)
{
    memcpy(p0, w->p, (w->n + w->nic) * sizeof(solnp_float));
    SOLNP(add_scaled_array)
    (p0, g, w->n + w->nic, -stepsize);
    return;
}

solnp_float *sub_trustreg(
    solnp_float *Q,
    solnp_float *c,
    solnp_float *G,
    solnp_float delta,
    solnp_float tol)
{
    // This subroutine solve min  1/2*a'*Q*a + c'*a
    /* s.t.a'*G*a <= delta^2        (dual  <------  lam >= 0)

         the first - order condition of nonconvex case is
         (Q + lam * G)* a + c = 0
         a'*G*a <= delta^2
         lam >= 0
         lam * (a'*G*a - delta^2) = 0

         second - order condition
         (Q + lam * G) >= 0
         REMEMBER TO FREE the RETURN pointer
    */
    solnp_float *alpha_mval = (solnp_float *)solnp_malloc(3 * sizeof(solnp_float));
    solnp_int dim = 2;

    // Lower bound of Multiplier
    solnp_float a = G[0] * G[3];
    solnp_float b = G[3] * Q[0] + Q[3] * G[0];
    solnp_float c_1 = Q[0] * Q[3] - (G[1] + Q[1]) * (G[1] + Q[1]);
    solnp_float lam_lb = MAX((-b + sqrt(b * b - 4 * a * c_1)) / (2 * a), 0);
    solnp_float lam;

    // Uppder bound of Multiplier
    solnp_float lam_ub = lam_lb + 1e8;
    solnp_int tag = 1;
    // Bisection
    while ((lam_ub - lam_lb) >= tol && tag)
    {
        if (lam_lb >= 1e8)
        {
            lam = lam_lb * 1.3;
            tag = 0;
        }
        else
        {
            lam = (lam_lb + lam_ub) / 2;
        }
        solnp_float *q_lam_G = (solnp_float *)solnp_malloc(sizeof(solnp_float) * 4);
        memcpy(q_lam_G, Q, sizeof(solnp_float) * 4);
        SOLNP(add_scaled_array)
        (q_lam_G, G, 4, lam);

        memcpy(alpha_mval, c, 2 * sizeof(solnp_float));
        SOLNP(chol)
        (2, q_lam_G);
        SOLNP(solve_lin_sys)
        (2, 1, q_lam_G, alpha_mval);

        solnp_free(q_lam_G);

        solnp_float aGa = 0;
        for (solnp_int i = 0; i < 2; i++)
        {
            alpha_mval[i] = -alpha_mval[i];
            for (solnp_int j = 0; j < 2; j++)
            {
                aGa += alpha_mval[i] * alpha_mval[j] * G[i + 2 * j];
            }
        }

        if (aGa > delta * delta)
        {
            lam_lb = lam;
        }
        else if (aGa < delta * delta)
        {
            lam_ub = lam;
        }
        else
        {
            break;
        }
    }

    alpha_mval[2] = 0;
    for (solnp_int i = 0; i < 2; i++)
    {
        alpha_mval[2] += c[i] * alpha_mval[i];
        for (solnp_int j = 0; j < 2; j++)
        {
            alpha_mval[2] += alpha_mval[i] * alpha_mval[j] * Q[i + 2 * j] / 2;
        }
    }

    return alpha_mval;
}

solnp_float drsom(
    SOLNPWork *w,
    SUBNPWork *w_sub,
    SOLNPSettings *stgs,
    solnp_float *p0,
    solnp_float *g,
    solnp_float *m,
    solnp_float radius,
    solnp_float val)
{
    // This subroutine implement the drsom update.
    // g is the gradient direction and m is the momentum direction
    // radius is the interpolation radius.

    // calculate the two orthgonal direction
    solnp_float *d1 = (solnp_float *)solnp_malloc(w->n * sizeof(solnp_float));
    solnp_float *d2 = (solnp_float *)solnp_malloc(w->n * sizeof(solnp_float));
    memcpy(d1, g, w->n * sizeof(solnp_float));
    memcpy(d2, m, w->n * sizeof(solnp_float));
    solnp_set_as_scaled_array(d1, d1, 1 / solnp_norm(d1, w->n), w->n);
    solnp_add_scaled_array(d2, d1, w->n, -solnp_dot(d1, d2, w->n));
    solnp_set_as_scaled_array(d2, d2, 1 / solnp_norm(d2, w->n), w->n);

    // interpolation in two directions
    solnp_float *g_and_Q = interpolate2d(w, stgs, w_sub, radius, d1, d2, p0, val);
    solnp_float *g_2d = (solnp_float *)solnp_malloc(2 * sizeof(solnp_float));
    solnp_float *Q_2d = (solnp_float *)solnp_malloc(4 * sizeof(solnp_float));
    memcpy(g_2d, g_and_Q, 2 * sizeof(solnp_float));
    memcpy(Q_2d, g_and_Q + 2, 4 * sizeof(solnp_float));

    // calculate the trust-region subproblem
    solnp_float *Id = (solnp_float *)solnp_calloc(4, sizeof(solnp_float));
    Id[0] = 1;
    Id[2] = 1;
    solnp_float *alpha_mval = sub_trustreg(Q_2d, g_2d, Id, radius, 1e-6);

    // calculate the new point
    memcpy(p0, w->p, w->n * sizeof(solnp_float));
    solnp_add_scaled_array(p0, d1, w->n, alpha_mval[0]);
    solnp_add_scaled_array(p0, d2, w->n, alpha_mval[1]);
    solnp_float mval = alpha_mval[2];

    // Free variables
    solnp_free(g_and_Q);
    solnp_free(d1);
    solnp_free(d2);
    solnp_free(g_2d);
    solnp_free(Q_2d);
    solnp_free(Id);
    solnp_free(alpha_mval);

    return mval;
}
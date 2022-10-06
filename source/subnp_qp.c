#include "subnp.h"
#include "linalg.h"
#include "linsys.h"
#include <mkl.h>
#include "solnp_util.h"
#include "osqp.h"
#include <math.h>
solnp_float qpsolver_time = 0;

void copySOLNPCost(SOLNPCost* ob1, SOLNPCost* ob2) {
    ob1->obj = ob2->obj;
    if (ob2->nec) {
        memcpy(ob1->ec, ob2->ec, ob2->nec * sizeof(solnp_float));
    }
    else {
        ob1->ec = SOLNP_NULL;
    }
    if (ob2->nic) {
        memcpy(ob1->ic, ob2->ic, ob2->nic * sizeof(solnp_float));
    }
    else {
        ob1->ic = SOLNP_NULL;
    }
}

SOLNPCost* init_cost(solnp_int nec, solnp_int nic)
{
    SOLNPCost* c = (SOLNPCost*)solnp_calloc(1, sizeof(SOLNPCost));
    c->nec = nec;
    c->nic = nic;
    if (nec > 0) {
        c->ec = (solnp_float*)solnp_malloc(nec * sizeof(solnp_float));
    }
    else {
        c->ec = SOLNP_NULL;
    }

    if (nic > 0) {
        c->ic = (solnp_float*)solnp_malloc(nic * sizeof(solnp_float));
    }
    else {
        c->ic = SOLNP_NULL;
    }

    c->cost = SOLNP_NULL;

    return c;
}

SUBNPWork *init_work_subnp
(
    SOLNPWork *w,
    SOLNPSettings *stgs
)
{
    solnp_int i, j;
    SUBNPWork *w_sub = (SUBNPWork *)solnp_malloc(sizeof(SUBNPWork));
    
    // calloc alp
    w_sub->alp = (solnp_float *)solnp_calloc(3, sizeof(solnp_float));

    // malloc scale
    w_sub->n_scale = 1 + w->nec + w->nic + w->n;
    w_sub->scale = (solnp_float *)solnp_malloc(w_sub->n_scale * sizeof(solnp_float));
    // init ch
    w_sub->ch = 1;
    w_sub->ob_cand = init_cost(w->nec, w->nic);
    w_sub->ob_cand->obj = INFINITY;
    // init mm
    if (w->pb->Ipb[1] >= 0.5) {
        if (w->pb->Ipb[0] <= 0.5) {
            w_sub->mm = w->nic;
        } else {
            w_sub->mm = w->nic + w->n;
        }
    } else {
        w_sub->mm = 0;
    }

    // init nc and npic
    w_sub->nc = w->nc;
    w_sub->npic = w->nic + w->n;

    // malloc Jacob
    w_sub->J = (SUBNPJacob *)solnp_malloc(sizeof(SUBNPJacob));
    w_sub->J->n = w->n;
    w_sub->J->nec = w->nec;
    w_sub->J->nic = w->nic;
    w_sub->J->nc  = w->nc;
    w_sub->J->npic = w->nic + w->n;
    w_sub->J->g = (solnp_float *)solnp_calloc(w_sub->J->npic, sizeof(solnp_float));
    w_sub->J->a = (solnp_float *)solnp_calloc(w_sub->nc * (w_sub->npic + 1), sizeof(solnp_float));
    w_sub->p_cand = (solnp_float*)solnp_malloc((w->n + w->nic) * sizeof(solnp_float));

    // malloc b
    if (w->nc > 0.5) {
        w_sub->b = (solnp_float *)solnp_malloc(w->nc * sizeof(solnp_float));
    } else {
        w_sub->b = SOLNP_NULL;
    }

    return w_sub;
}

SUBNPWork *SUBNP(init)
(
    SOLNPWork *w,
    SOLNPSettings *stgs
)
{
    SUBNPWork *w_sub;
    w_sub = init_work_subnp(w, stgs);

    return w_sub;
}

solnp_float *compute_scale
(
    SOLNPWork *w,
    SOLNPSettings *stgs
)
{
    solnp_int i, n_scale = 1 + w->nec + w->nic + w->n;
    solnp_float *scale = (solnp_float *)solnp_malloc(n_scale * sizeof(solnp_float));

    // objective value and equality constraints
    if (w->nec > 0){
        scale[0] = w->ob->obj;
        solnp_float ec_norm_inf = SOLNP(norm_inf)(w->ob->ec, w->nec);
        for (i = 0; i < w->nec; i++) {
            scale[1 + i] = ec_norm_inf;
        }
    } else {
        scale[0] = 1.;
    }

    // inequality constraints and bounds of decision variables
    if (w->pb->Ipb[1] <= 0) {
        memcpy(&(scale[1 + w->nec]), w->p, (w->nic+w->n)*sizeof(solnp_float));
    } else {
        for (i = 0; i < (w->nic + w->n); i++) {
            scale[1 + w->nec + i] = 1.;
        }
    }

    // bound scale
    for (i = 0; i < n_scale; i++){
        scale[i] = MIN(MAX(ABS(scale[i]), stgs->tol), SAFEDIV_POS(1., stgs->tol));
    }

    return scale;
}

solnp_int rescale_ob
(
    SOLNPCost *ob,
    const solnp_float *scale
)
{
    solnp_int i;
    ob->obj /= scale[0];
    for (i = 0; i < ob->nec; i++) {
        ob->ec[i] /= scale[1 + i];
    }
    for (i = 0; i < ob->nic; i++) {
        ob->ic[i] /= scale[1 + ob->nec + i];
    }

    return 0;
}

solnp_int rescale_p
(
    solnp_float *p,
    solnp_float *scale,
    SOLNPWork *w
)
{
    solnp_int i;
    for (i = 0; i < (w->nic + w->n); i++) {
        p[i] /= scale[1 + w->nec + i];
    }

    return 0;
}

solnp_int rescale
(
    SUBNPWork *w_sub,
    SOLNPWork *w,
    solnp_float *scale,
    SOLNPSettings *stgs
)
{
    solnp_int i, j;

    // rescale ob
    rescale_ob(w->ob, scale);

    // rescale p
    rescale_p(w->p, scale, w);

    // rescale pb
    if (w->pb->Ipb[1] >= 0.5) {
        for (i = 0; i < w->nic; i++) {
            w->pb->il[i] /= scale[1 + w->nec + i];
            w->pb->iu[i] /= scale[1 + w->nec + i];
        }
        if (w->pb->Ipb[0] > 0.5) {
            for (i = 0; i < w->n; i++) {
                w->pb->pl[i] /= scale[1 + w->nec + w->nic + i];
                w->pb->pu[i] /= scale[1 + w->nec + w->nic + i];
            }
        }
    }

    // rescale lagrangian multipliers
    if (w->nc > 0.5) {
        for (i = 0; i < w->nc; i++) {
            w->l[i] = scale[1 + i] * w->l[i] / scale[0];
        }
    }

    // rescale Hessian
    for (i = 0; i < w_sub->npic; i++) {
        for (j = 0; j <= i; j++) {
            w->h[i + j * w_sub->npic] = scale[1 + w->nec + i] * w->h[i + j * w_sub->npic] * scale[1 + w->nec + j] / scale[0];

            // use symmetric
            if (j < i) {
                w->h[j + i * w_sub->npic] = w->h[i + j * w_sub->npic];
            }
        }
    }

    // record scale in w_sub
    for (i = 0; i < w_sub->n_scale; i++) {
        w_sub->scale[i] *= scale[i];
    }

    return 0;
}

void unscale
(
    SUBNPWork* w_sub,
    SOLNPWork* w,
    SOLNPSettings* stgs
){
    solnp_int i;
    for (i = 0; i < w_sub->n_scale; i++) {
        w_sub->scale[i] = 1 / w_sub->scale[i];
    }
    rescale(w_sub, w, w_sub->scale, stgs);
}

solnp_int update_work_subnp
(
    SUBNPWork *w_sub,
    SOLNPWork *w,
    SOLNPSettings *stgs
)
{
    solnp_int i;

    // init scale
    for (i = 0; i < w_sub->n_scale; i++) {
        w_sub->scale[i] = 1.;
    }

    // init Jacob
    for (i = 0; i < w->nic; i++) {
        // a[k+nec, k] = -1
        w_sub->J->a[i + w->nec + i * w_sub->nc] = -1.;
    }

    solnp_float *scale = compute_scale(w, stgs);
    rescale(w_sub, w, scale, stgs);

    solnp_free(scale);
    return 0;
}

solnp_float calculate_infeas_scaledob
(
    SOLNPCost* ob,
    SOLNPWork* w,
    solnp_float* x,
    solnp_float* scale
) {
    solnp_float con = 0;
    solnp_int i;
    //calculate constraints
    if (w->nic > 0.5) {
        for (i = 0; i < w->nic; i++) {
            if (ob->ic[i] < w->pb->il[i]) {
                con += (w->pb->il[i] - ob->ic[i]) * (w->pb->il[i] - ob->ic[i])*scale[1 + w->nec + i] * scale[1 + w->nec + i];
            }
            if (ob->ic[i] > w->pb->iu[i] ) {
                con += (ob->ic[i] - w->pb->iu[i] ) * (ob->ic[i] - w->pb->iu[i]) * scale[1 + w->nec + i] * scale[1 + w->nec + i];
            }
        }
    }
    if (w->nec > 0.5) {
        for (i = 0; i < w->nec; i++) {
            con += ob->ec[i] * ob->ec[i] * scale[1 + i] * scale[1 + i];
        }
    }
    for (i = 0; i < w->n; i++) {
        if (x[i] > w->pb->pu[i] ) {
            con += (x[i] - w->pb->pu[i]) * (x[i] - w->pb->pu[i] ) * scale[1 + w->nec + w->nic + i] * scale[1 + w->nec + w->nic + i];
        }
        if (x[i] < w->pb->pl[i] ) {
            con += (x[i] - w->pb->pl[i] ) * (x[i] - w->pb->pl[i] ) * scale[1 + w->nec + w->nic + i] * scale[1 + w->nec + w->nic + i];
        }
    }
    con = SQRTF(con);
    return con;
}
solnp_float calculate_infeas_unscaleob
(
    SOLNPCost* ob,
    SOLNPWork* w,
    solnp_float* x,
    solnp_float* scale
) {
    solnp_float con = 0;
    solnp_int i;
    //calculate constraints
    if (w->nic > 0.5) {
        for (i = 0; i < w->nic; i++) {
            if (ob->ic[i] < w->pb->il[i] * scale[1 + w->nec + i]) {
                con += (w->pb->il[i] * scale[1 + w->nec + i] - ob->ic[i]) * (w->pb->il[i] * scale[1 + w->nec + i] - ob->ic[i]);
            }
            if (ob->ic[i] > w->pb->iu[i] * scale[1 + w->nec + i]) {
                con += (ob->ic[i] - w->pb->iu[i] * scale[1 + w->nec + i]) * (ob->ic[i] - w->pb->iu[i] * scale[1 + w->nec + i]);
            }
        }
    }
    if (w->nec > 0.5) {
        for (i = 0; i < w->nec; i++) {
            con += ob->ec[i] * ob->ec[i];
        }
    }
    for (i = 0; i < w->n; i++) {
        if (x[i] > w->pb->pu[i] * scale[1 + w->nec + w->nic + i]) {
            con += (x[i] - w->pb->pu[i] * scale[1 + w->nec + w->nic + i]) * (x[i] - w->pb->pu[i] * scale[1 + w->nec + w->nic + i]);
        }
        if (x[i] < w->pb->pl[i] * scale[1 + w->nec + w->nic + i]) {
            con += (x[i] - w->pb->pl[i] * scale[1 + w->nec + w->nic + i]) * (x[i] - w->pb->pl[i] * scale[1 + w->nec + w->nic + i]);
        }
    }
    con = SQRTF(con);
    return con;
}
solnp_float calculate_almcrit_iq
(
    SOLNPWork* w,
    solnp_float* scale,
    solnp_float* p
) {
    solnp_float almcrit = 0;
    solnp_float temp;
    solnp_int i;
    for (i = 0; i < w->nic; i++) {
        temp = p[i] * scale[w->nec + 1 + i] - scale[0] * w->l[w->nec + i] / scale[w->nec + 1 + i];
        if (temp > w->pb->iu[i] * scale[w->nec + 1 + i]) {
            temp = w->pb->iu[i] * scale[w->nec + 1 + i];
        }
        else if (temp < w->pb->il[i] * scale[w->nec + 1 + i]) {
            temp = w->pb->il[i] * scale[w->nec + 1 + i];
        }
        almcrit += (temp - p[i] * scale[w->nec + 1 + i]) * (temp - p[i] * scale[w->nec + 1 + i]);
    }
    return almcrit;
}
void proj
(
    solnp_float* p,
    SOLNPWork* w
) {
    solnp_int i;
    for (i = 0; i < w->nic; i++) {
        if (p[i] > w->pb->iu[i]) {
            p[i] = w->pb->iu[i];
        }if (p[i] < w->pb->il[i]) {
            p[i] = w->pb->il[i];
        }
    }for (i = 0; i < w->n; i++) {
        if (p[w->nic + i] > w->pb->pu[i]) {
            p[w->nic + i] = w->pb->pu[i];
        }if (p[w->nic + i] < w->pb->pl[i]) {
            p[w->nic + i] = w->pb->pl[i];
        }
    }
}
void calculate_scaled_cost
(
    SOLNPCost *ob,
    solnp_float *p,
    const solnp_float *scale,
    SOLNPSettings* stgs,
    SOLNPWork *w
)
{
    solnp_int i;
    solnp_float *x = (solnp_float *)solnp_malloc(w->n * sizeof(solnp_float));
    
    for (i = 0; i < w->n; i++) {
        x[i] = p[w->nic + i] * scale[1 + w->nc + i];
    }
    // calculate rescaled cost
    // (*w->ob->cost)(ob, x, w->n, 0);
    w->ob->cost(ob, x, w->n);
    w->count_cost += 1;

    solnp_float con = calculate_infeas_unscaleob(ob, w, x, scale);
    
    if (con <= stgs->tol_con && w->bestobj > ob->obj) {
        w->bestcon = con;
        w->bestobj = ob->obj;
        memcpy(w->bestp, x, (w->n) * sizeof(solnp_float));
        memcpy(w->bestl, w->l, MAX(w->nc,1) * sizeof(solnp_float));
    }
    rescale_ob(ob, scale);

    if (w->count_cost >= stgs->maxfev) {
        w->exit = 1;
    }

    // free x
    solnp_free(x);
}

solnp_int calculate_Jacob
(
    SUBNPWork *w_sub,
    SOLNPWork *w,
    SOLNPSettings *stgs
)
{
    SOLNPCost *ob = init_cost(w->nec, w->nic);
    solnp_float *p = (solnp_float *)solnp_malloc((w->nic + w->n) * sizeof(solnp_float));
    memcpy(p, w->p, w_sub->npic * sizeof(solnp_float));
    solnp_float infeas_cand,temp;
    solnp_int i, j;

    for (i = 0; i < w->n; i++) {
        // perturb
        p[w->nic + i] = p[w->nic + i] + stgs->delta;
        calculate_scaled_cost(ob, p, w_sub->scale,stgs, w);

        infeas_cand = calculate_infeas_scaledob(ob, w, &p[w->nic], w_sub->scale);
        if (infeas_cand <= stgs->tol_con && MIN(p[w->nic+i]-w->pb->pl[i], w->pb->pu[i]- p[w->nic + i])>0 && w_sub->ob_cand->obj > ob->obj) {
            memcpy(w_sub->p_cand, p, w_sub->npic * sizeof(solnp_float));
            copySOLNPCost(w_sub->ob_cand, ob);
        }

        // calculate Jacobian approximately
        w_sub->J->g[w->nic + i] = (ob->obj - w->ob->obj) /  stgs->delta;
        for (j = 0; j < w->nec; j++) {
            w_sub->J->a[j + (w->nic + i) * w_sub->nc] = (ob->ec[j] - w->ob->ec[j]) / stgs->delta;
        }
        for (j = 0; j < w->nic; j++) {
            w_sub->J->a[(w->nec + j) + (w->nic + i) * w_sub->nc] = (ob->ic[j] - w->ob->ic[j]) / stgs->delta;
        }

        // restore perturb
        p[w->nic + i] = p[w->nic + i] - stgs->delta;

        temp = p[w->nic + i] * w_sub->scale[w->nc + 1 + i] - w_sub->J->g[w->nic + i] * w_sub->scale[0] / w_sub->scale[w->nc + 1 + i];
        if (w->nc > 0) {
            temp += w_sub->scale[0]* SOLNP(dot)(&w_sub->J->a[(w->nic + i) * w_sub->nc],w->l,w->nc) / w_sub->scale[w->nc + 1 + i];
        }
        if (temp > w->pb->pu[i] * w_sub->scale[w->nc + 1 + i]) {
            temp = w->pb->pu[i] * w_sub->scale[w->nc + 1 + i];
        }
        else if (temp < w->pb->pl[i] * w_sub->scale[w->nc + 1 + i]) {
            temp = w->pb->pl[i] * w_sub->scale[w->nc + 1 + i];
        }
        w->alm_crit += (temp - p[w->nic + i] * w_sub->scale[w->nc + 1 + i]) * (temp - p[w->nic + i] * w_sub->scale[w->nc + 1 + i]);
        if (w->exit == 1) {
            break;
        }
    }
    solnp_float* aT = (solnp_float*)solnp_malloc(w_sub->nc * w_sub->npic * sizeof(solnp_float));
    SOLNP(transpose)(w_sub->nc, w_sub->npic, w_sub->J->a, aT);
    solnp_float* aaT = (solnp_float*)solnp_malloc(w_sub->nc * w_sub->nc * sizeof(solnp_float));
    SOLNP(AB)(aaT, w_sub->J->a, aT, w_sub->nc, w_sub->npic, w_sub->nc);
    solnp_float *cond = (solnp_float*)solnp_malloc(sizeof(solnp_float));
    SOLNP(cond)(w_sub->nc, aaT,cond);
    solnp_free(aT);
    solnp_free(aaT);
    // calculate condition number
    // TODO: calculate condition number
     if (*cond <=  EPS) {
         solnp_printf("SOLNP+--> ");
         solnp_printf("Redundant constraints were found. Poor              \n");
         solnp_printf("         ");
         solnp_printf("intermediate results may result.  Suggest that you  \n");
         solnp_printf("         ");
         solnp_printf("remove redundant constraints and re-OPTIMIZE.       \n");
     }
    // free pointers
     solnp_free(cond);
    solnp_free(p);
    free_cost(ob);

    return 0;
}
solnp_int qpsolver
(
    solnp_int n,
    solnp_float* H,
    solnp_float* g,
    solnp_int m1,
    solnp_int m2,
    solnp_float* A,
    solnp_float* b,
    solnp_float* lb,
    solnp_float* ub,
    solnp_float* p,
    solnp_float* lambda
) {
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
    c_int  P_nnz, A_nnz;
    c_float* P_x;
    c_int* P_i;
    c_int* P_p;

    SOLNP(timer) qpsolver_timer;
    SOLNP(tic)(&qpsolver_timer);

    /*Setup two matrice*/
    if (H == SOLNP_NULL) {
        P_nnz = 0;
        P_x = OSQP_NULL;
        P_i = OSQP_NULL;
        P_p = c_malloc((n + 1) * sizeof(c_int));
        for (i = 0; i < n + 1; i++) {
            P_p[i] = 0;
        }
    }
    else {
        P_nnz = countA_sys(n, n, H);
        P_p = (c_int*)c_malloc((n + 1) * sizeof(c_int));
        P_i = (c_int*)c_malloc(P_nnz * sizeof(c_int));
        P_x = (c_float*)c_malloc(P_nnz * sizeof(c_float));
        calculate_csc_sys(n, n, H, P_x, P_i, P_p);
    }
    A_nnz = countA(m, n, A);
    c_int* A_p = (c_int*)c_malloc((n + 1) * sizeof(c_int));
    c_int* A_i = (c_int*)c_malloc(A_nnz * sizeof(c_int));
    c_float* A_x = (c_float*)c_malloc(A_nnz * sizeof(c_float));
    calculate_csc(m, n, A, A_x, A_i, A_p);

    c_float* q = (c_float*)c_malloc(n * sizeof(c_float));
    memcpy(q, g, n * sizeof(c_float));

    c_float* l = (c_float*)c_malloc((m + n) * sizeof(c_float));
    memcpy(l, b, m * sizeof(c_float));
    memcpy(&l[m], lb, n * sizeof(c_float));

    c_float* u = (c_float*)c_malloc((m + n) * sizeof(c_float));
    memcpy(u, b, m * sizeof(c_float));
    memcpy(&u[m], ub, n * sizeof(c_float));

    // Exitflag
    c_int exitflag = 0;

    // Workspace structures
    OSQPWorkspace* work;
    OSQPSettings* settings = (OSQPSettings*)c_malloc(sizeof(OSQPSettings));
    OSQPData* data = (OSQPData*)c_malloc(sizeof(OSQPData));

    // Populate data
    if (data) {
        data->n = n1;
        data->m = m + n1;
        data->P = csc_matrix(data->n, data->n, P_nnz, P_x, P_i, P_p);
        data->q = q;
        data->A = csc_matrix(data->m, data->n, A_nnz, A_x, A_i, A_p);
        data->l = l;
        data->u = u;
    }

    // Define solver settings as default
    if (settings) osqp_set_default_settings(settings);
    settings->eps_abs = 1e-6;
    settings->eps_rel = 1e-6;
    // Setup workspace
    exitflag = osqp_setup(&work, data, settings);



    // Solve Problem
    osqp_solve(work);

    memcpy(p, work->solution->x, n * sizeof(c_float));
    if (lambda) {
        memcpy(lambda, work->solution->y, m1 * sizeof(c_float));
        SOLNP(scale_array)(lambda, -1.0, m1);
    }

    // Clean workspace
    osqp_cleanup(work);
    if (data) {
        if (data->A) c_free(data->A);
        if (data->P) c_free(data->P);
        c_free(q);
        c_free(l);
        c_free(u);
        c_free(data);
    }
    c_free(A_x);
    c_free(A_i);
    c_free(A_p);

    if (P_x) {
        c_free(P_x);
    }
    if (P_i) {
        c_free(P_i);
    }
    c_free(P_p);
    if (settings)  c_free(settings);


    qpsolver_time += SOLNP(tocq)(&qpsolver_timer) / 1e3;

    return exitflag;
}
solnp_int find_int_feas_sol_osqp
(
    SUBNPWork* w_sub,
    SOLNPWork* w,
    SOLNPSettings* stgs
)
{
    w_sub->ch = -1;
    w_sub->alp[0] = stgs->tol - SOLNP(norm_inf)(&w_sub->J->a[w_sub->npic * w_sub->nc], w_sub->nc);

    if (w_sub->alp[0] <= 0) {
        w_sub->ch = 1;
        if (w->pb->Ipb[1] == 0) {
            // TODO: project onto hyperplane
            // need linear system solver
            // p0=p0-a'*((a*a')\constraint);

            // p0: w->p is a w_sub->npic vector
            // a:  w_sub->J->a is a w_sub->nc * w_sub->npic matrix
            //     in col major form
            // constraint: &w_sub->J->a[w_sub->nc * w_sub->npic] is a w_sub->nc vector


            solnp_float* aT = (solnp_float*)solnp_malloc(w_sub->nc * w_sub->npic * sizeof(solnp_float));
            SOLNP(transpose)(w_sub->nc, w_sub->npic, w_sub->J->a, aT);
            solnp_float* aaT = (solnp_float*)solnp_malloc(w_sub->nc * w_sub->nc * sizeof(solnp_float));
            SOLNP(AB)(aaT, w_sub->J->a, aT, w_sub->nc, w_sub->npic, w_sub->nc);
            SOLNP(chol)(w_sub->nc, aaT);
            solnp_float* inv_aaT_cons = (solnp_float*)solnp_malloc(w_sub->nc * sizeof(solnp_float));
            memcpy(inv_aaT_cons, &w_sub->J->a[w_sub->nc * w_sub->npic], w_sub->nc * sizeof(solnp_float));
            SOLNP(solve_lin_sys)(w_sub->nc, 1, aaT, inv_aaT_cons);
            solnp_float* aT_inv_aaT_cons = (solnp_float*)solnp_malloc(w_sub->npic * sizeof(solnp_float));
            SOLNP(Ax)(aT_inv_aaT_cons, aT, inv_aaT_cons, w_sub->npic, w_sub->nc);
            SOLNP(add_scaled_array)(w->p, aT_inv_aaT_cons, w_sub->npic, -1.0);

            solnp_free(aaT);
            solnp_free(aT);
            solnp_free(inv_aaT_cons);
            solnp_free(aT_inv_aaT_cons);


            w_sub->alp[0] = 1;
        }
        else {
            // TODO: affine scaling method to solve LP
            solnp_float* p = (solnp_float*)solnp_calloc((w_sub->npic + 1), sizeof(solnp_float));
            solnp_float* c = (solnp_float*)solnp_calloc((w_sub->npic + 1), sizeof(solnp_float));
            solnp_float* y = (solnp_float*)solnp_calloc((w_sub->nc), sizeof(solnp_float));
            solnp_float* lb = (solnp_float*)solnp_calloc((w_sub->npic + 1), sizeof(solnp_float));
            solnp_float* ub = (solnp_float*)solnp_calloc((w_sub->npic + 1), sizeof(solnp_float));
            // init p, c, dx

            c[w_sub->npic] = 1.;
            lb[w_sub->npic] = 0;
            ub[w_sub->npic] = INFINITY;
            memcpy(lb, w->pb->il, w->pb->nic * sizeof(solnp_float));
            memcpy(ub, w->pb->iu, w->pb->nic * sizeof(solnp_float));
            memcpy(&lb[w->pb->nic], w->pb->pl, w->pb->n * sizeof(solnp_float));
            memcpy(&ub[w->pb->nic], w->pb->pu, w->pb->n * sizeof(solnp_float));

            // last column of a should be **negative** constraint
            SOLNP(scale_array)(&w_sub->J->a[w->nc * w_sub->npic], -1.0, w->nc);

            qpsolver(w_sub->npic + 1, SOLNP_NULL, c, w->nc, w_sub->npic + 1, w_sub->J->a, w_sub->b, lb, ub, p, SOLNP_NULL);

            proj(p, w);
            memcpy(w->p, p, w_sub->npic * sizeof(solnp_float));
            SOLNP(Ax)(w_sub->b, w_sub->J->a, w->p, w->nc, w_sub->npic);

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


solnp_int find_int_feas_sol_aff
(
    SUBNPWork *w_sub,
    SOLNPWork *w,
    SOLNPSettings *stgs
)
{

    SOLNP(timer) qpsolver_timer;
    SOLNP(tic)(&qpsolver_timer);
    w_sub->ch = -1;
    w_sub->alp[0] = stgs->tol - SOLNP(norm_inf)(&w_sub->J->a[w_sub->npic * w_sub->nc], w_sub->nc);

    if (w_sub->alp[0] <= 0) {
        w_sub->ch = 1;
        if (w->pb->Ipb[1] == 0) {
            // TODO: project onto hyperplane
            // need linear system solver
            // p0=p0-a'*((a*a')\constraint);

            // p0: w->p is a w_sub->npic vector
            // a:  w_sub->J->a is a w_sub->nc * w_sub->npic matrix
            //     in col major form
            // constraint: &w_sub->J->a[w_sub->nc * w_sub->npic] is a w_sub->nc vector


            solnp_float *aT = (solnp_float*)solnp_malloc(w_sub->nc * w_sub->npic * sizeof(solnp_float));
            SOLNP(transpose)(w_sub->nc, w_sub->npic, w_sub->J->a, aT);
            solnp_float *aaT = (solnp_float*)solnp_malloc(w_sub->nc * w_sub->nc * sizeof(solnp_float));
            SOLNP(AB)(aaT, w_sub->J->a, aT, w_sub->nc, w_sub->npic, w_sub->nc);
            SOLNP(chol)(w_sub->nc, aaT);
            solnp_float *inv_aaT_cons = (solnp_float*)solnp_malloc(w_sub->nc * sizeof(solnp_float));
            memcpy(inv_aaT_cons, &w_sub->J->a[w_sub->nc * w_sub->npic], w_sub->nc * sizeof(solnp_float));
            SOLNP(solve_lin_sys)(w_sub->nc, 1, aaT, inv_aaT_cons);
            solnp_float *aT_inv_aaT_cons = (solnp_float*)solnp_malloc(w_sub->npic * sizeof(solnp_float));
            SOLNP(Ax)(aT_inv_aaT_cons, aT, inv_aaT_cons, w_sub->npic, w_sub->nc);
            SOLNP(add_scaled_array)(w->p, aT_inv_aaT_cons, w_sub->npic, -1.0);

            solnp_free(aaT);
            solnp_free(aT);
            solnp_free(inv_aaT_cons);
            solnp_free(aT_inv_aaT_cons);


            w_sub->alp[0] = 1;
        } else {
            // TODO: affine scaling method to solve LP
            solnp_float* p = (solnp_float*)solnp_calloc((w_sub->npic + 1), sizeof(solnp_float));
            solnp_float* c = (solnp_float*)solnp_calloc((w_sub->npic + 1), sizeof(solnp_float));
            solnp_float* dx = (solnp_float*)solnp_calloc((w_sub->npic + 1), sizeof(solnp_float));
            solnp_float* y = (solnp_float*)solnp_calloc((w_sub->nc), sizeof(solnp_float));
            solnp_float* v = (solnp_float*)solnp_calloc((w_sub->npic + 1), sizeof(solnp_float));
            solnp_int i, j, minit = 0;
            solnp_float go = 1.0, z;

            // init p, c, dx
            memcpy(p, w->p, w_sub->npic * sizeof(solnp_float));
            p[w_sub->npic] = 1.;
            c[w_sub->npic] = 1.;

            // last column of a should be **negative** constraint
            SOLNP(scale_array)(&w_sub->J->a[w->nc * w_sub->npic], -1.0, w->nc);

            for (i = 0; i <= w_sub->npic; i++) {
                dx[i] = 1.;
            }
            // POTENTIAL PROBLEM: the last column should be -constraint? -Right. This bug has been fixed.

            // affine scaling while loop
            while (go >= stgs->tol)
            {
                minit += 1;

                // form ellipse
                for (i = 0; i < w->nic; i++) {
                    dx[i] = MIN(p[i] - w->pb->il[i], w->pb->iu[i] - p[i]);
                }
                dx[w_sub->npic] = p[w_sub->npic];
                if (w->pb->Ipb[0] == 0) { // x is free
                    for (i = 0; i < w->n; i++) {
                        solnp_float dx_norm_inf = SOLNP(norm_inf)(dx, w->nic);
                        dx_norm_inf = MAX(dx_norm_inf, 100);
                        dx[w->nic + i] = dx_norm_inf;
                    }
                }
                else { // x is not free
                    for (i = 0; i < w->n; i++) {
                        dx[w->nic + i] = MIN(p[w->nic + i] - w->pb->pl[i], w->pb->pu[i] - p[w->nic + i]);
                    }
                }

                // TODO: solve for direction v
                // y=(a*diag(dx))'\(dx.*c');
                // v=dx.*(dx.*(c'-a'*y));

                // a: w_sub->J->a is a w_sub->nc * (w_sub->npic+1) matrix
                //    in col major form
                // c': w_sub->npic+1 col vector (c is row vector in matlab)
                // dx: w_sub->npic+1 col vector
                // y: w_sub->nc col vector
                // v: w_sub->npic+1 col vector


                solnp_float* a_diag_dx = (solnp_float*)solnp_malloc(w_sub->nc * (w_sub->npic + 1) * sizeof(solnp_float));
                memcpy(a_diag_dx, w_sub->J->a, w_sub->nc * (w_sub->npic + 1) * sizeof(solnp_float));


                for (i = 0; i < w_sub->npic + 1; i++) {
                    for (j = 0; j < w_sub->nc; j++) {
                        a_diag_dx[i * w_sub->nc + j] *= dx[i];
                    }
                }

                solnp_float* a_diag_dx_T = (solnp_float*)solnp_malloc(w_sub->nc * (w_sub->npic + 1) * sizeof(solnp_float));
                SOLNP(transpose)(w_sub->nc, w_sub->npic + 1, a_diag_dx, a_diag_dx_T);

                solnp_float* a_diag_dx2_aT = (solnp_float*)solnp_malloc(w_sub->nc * w_sub->nc * sizeof(solnp_float));
                SOLNP(AB)(a_diag_dx2_aT, a_diag_dx, a_diag_dx_T, w_sub->nc, w_sub->npic + 1, w_sub->nc);

                solnp_float* dx_c = (solnp_float*)solnp_malloc((w_sub->npic + 1) * sizeof(solnp_float));
                memcpy(dx_c, c, (w_sub->npic + 1) * sizeof(solnp_float));
                for (i = 0; i < w_sub->npic + 1; i++) {
                    dx_c[i] *= dx[i];
                }

                SOLNP(Ax)(y, a_diag_dx, dx_c, w_sub->nc, w_sub->npic + 1);
                SOLNP(chol)(w_sub->nc, a_diag_dx2_aT);
                SOLNP(solve_lin_sys)(w->nc, 1, a_diag_dx2_aT, y);
                solnp_float* aT = (solnp_float*)solnp_malloc(w_sub->nc * (w_sub->npic + 1) * sizeof(solnp_float));
                SOLNP(transpose)(w_sub->nc, w_sub->npic + 1, w_sub->J->a, aT);
                SOLNP(Ax)(v, aT, y, w_sub->npic + 1, w_sub->nc);
                SOLNP(add_scaled_array)(v, c, w_sub->npic + 1, -1.0);
                for (i = 0; i < w_sub->npic + 1; i++) {
                    v[i] *= -dx[i] * dx[i];
                }

                solnp_free(a_diag_dx);
                solnp_free(a_diag_dx_T);
                solnp_free(a_diag_dx2_aT);
                solnp_free(dx_c);
                solnp_free(aT);





                if (v[w_sub->npic] > 0) {
                    // calculate step size
                    z = p[w_sub->npic] / v[w_sub->npic];
                    for (i = 0; i < w->nic; i++) {
                        if (v[i] < 0) {
                            z = MIN(z, -(w->pb->iu[i] - p[i]) / v[i]);
                        }
                        else if (v[i] > 0)
                        {
                            z = MIN(z, (p[i] - w->pb->il[i]) / v[i]);
                        }
                    }
                    if (w->pb->Ipb[0] > 0) {
                        for (i = 0; i < w->n; i++) {
                            if (v[w->nic + i] < 0) {
                                z = MIN(z, -(w->pb->pu[i] - p[w->nic + i]) / v[w->nic + i]);
                            }
                            else if (v[w->nic + i] > 0)
                            {
                                z = MIN(z, (p[w->nic + i] - w->pb->pl[i]) / v[w->nic + i]);
                            }
                        }
                    }

                    // update
                    if (z >= (p[w_sub->npic] / v[w_sub->npic])) {
                        SOLNP(add_scaled_array)(p, v, w_sub->npic + 1, -z);
                    }
                    else {
                        SOLNP(add_scaled_array)(p, v, w_sub->npic + 1, -0.9 * z);
                    }
                    go = p[w_sub->npic];
                    if (minit >= 10) {
                        go = 0;
                    }
                }
                else {
                    go = 0;
                    minit = 10;
                }
            }

            if (minit >= 10) {
                solnp_printf("SOLNP+--> ");
                solnp_printf("The linearized problem has no feasible     \n");
                solnp_printf("         ");
                solnp_printf("solution.  The problem may not be feasible.\n");
            }

            memcpy(w->p, p, w_sub->npic * sizeof(solnp_float));
            SOLNP(Ax)(w_sub->b, w_sub->J->a, w->p, w->nc, w_sub->npic);

            // free pointers
            solnp_free(p);
            solnp_free(c);
            solnp_free(dx);
            solnp_free(y);
            solnp_free(v);
        }
    }
    qpsolver_time += SOLNP(tocq)(&qpsolver_timer) / 1e3;

    return 0;
}


solnp_float calculate_ALM
(
    SOLNPCost* ob,
    SOLNPSettings* stgs,
    solnp_float* p,
    const SOLNPWork* w,
    const SUBNPWork* w_sub
) { 
    
    
    solnp_float result = ob->obj;
    if (w->nc) {
        solnp_float* ap = (solnp_float*)solnp_malloc(w->nc * sizeof(solnp_float));
        SOLNP(Ax)(ap, w_sub->J->a, p, w->nc, w_sub->npic);
        if (w->nic) {
            solnp_float* ic = (solnp_float*)solnp_malloc(w->nic * sizeof(solnp_float));
            memcpy(ic, ob->ic, w->nic * sizeof(solnp_float));
            SOLNP(add_scaled_array)(ic, p, w->nic, -1.0);
            SOLNP(add_scaled_array)(ic, ap + w->nec, w->nic, -1.0);
            SOLNP(add_scaled_array)(ic, w_sub->b + w->nec, w->nic, 1.0);
            result = result + w->rho * SOLNP(norm_sq)(ic, w->nic);
            result = result - SOLNP(dot)(w->l + w->nec, ic, w->nic);
            solnp_free(ic);
        }
        if (w->nec) {
            solnp_float* ec = (solnp_float*)solnp_malloc(w->nec * sizeof(solnp_float));
            memcpy(ec, ob->ec, w->nec * sizeof(solnp_float));
            SOLNP(add_scaled_array)(ec, ap, w->nec, -1.0);
            SOLNP(add_scaled_array)(ec, w_sub->b, w->nec, 1.0);
            result = result + w->rho * SOLNP(norm_sq)(ec, w->nec);
            result = result - SOLNP(dot)(w->l, ec, w->nec);
            solnp_free(ec);
        }
        solnp_free(ap);
    }
    return result;
}
solnp_int calculate_ALMgradient
(
    SUBNPWork* w_sub,
    SOLNPWork* w,
    SOLNPSettings* stgs,
    solnp_float* g,
    solnp_float j
)
{
    SOLNPCost* obm = init_cost(w->nec, w->nic);
    solnp_float alm;
    solnp_int i;
    solnp_float infeas,temp;
    for (i = 0; i < w->n; i++) {
        solnp_float* contemp = (solnp_float*)solnp_malloc(w->nc * sizeof(solnp_float));
        // perturb
        w->p[w->nic + i] = w->p[w->nic + i] + stgs->delta;
        calculate_scaled_cost(obm, w->p, w_sub->scale,stgs, w);

        if (w->nec) {
            memcpy(contemp, obm->ec, w->nec * sizeof(solnp_float));
        }if (w->nic) {
            memcpy(&contemp[w->nec], obm->ic, w->nic * sizeof(solnp_float));
        }
       
        alm = obm->obj;
        if (w_sub->nc > 0) {
            alm = calculate_ALM(obm, stgs, w->p, w, w_sub);
        }

        infeas = calculate_infeas_scaledob(obm, w, &w->p[w->nic], w_sub->scale);
        if (infeas <= stgs->tol_con && MIN(w->p[w->nic + i] - w->pb->pl[i], w->pb->pu[i] - w->p[w->nic + i]) > 0 && w_sub->ob_cand->obj > obm->obj) {
            memcpy(w_sub->p_cand, w->p, w_sub->npic * sizeof(solnp_float));
            copySOLNPCost(w_sub->ob_cand, obm);
        }

        // restore perturb
        w->p[w->nic + i] = w->p[w->nic + i] - stgs->delta;
        // calculate gradient approximately
        g[w->nic + i] = (alm - j) / stgs->delta;
        if (w->exit == 1) {
            solnp_free(contemp);
            break;
        }


        temp = w->p[w->nic + i] * w_sub->scale[w->nc + 1 + i] - ((obm->obj-w->ob->obj)/stgs->delta) * w_sub->scale[0] / w_sub->scale[w->nc + 1 + i];
        if (w->nc > 0) {
            if (w->nec) {
                SOLNP(add_scaled_array)(contemp, w->ob->ec, w->nec, -1.0);
            }
            if (w->nic) {
                SOLNP(add_scaled_array)(&contemp[w->nec], w->ob->ic, w->nic, -1.0);
            }SOLNP(scale_array)(contemp, 1/stgs->delta, w->nc);
            temp += w_sub->scale[0] * SOLNP(dot)(contemp, w->l, w->nc) / w_sub->scale[w->nc + 1 + i];
        }
        if (temp > w->pb->pu[i] * w_sub->scale[w->nc + 1 + i]) {
            temp = w->pb->pu[i] * w_sub->scale[w->nc + 1 + i];
        }
        else if (temp < w->pb->pl[i] * w_sub->scale[w->nc + 1 + i]) {
            temp = w->pb->pl[i] * w_sub->scale[w->nc + 1 + i];
        }
        w->alm_crit += (temp - w->p[w->nic + i] * w_sub->scale[w->nc + 1 + i]) * (temp - w->p[w->nic + i] * w_sub->scale[w->nc + 1 + i]);
               
       
        solnp_free(contemp);
        
    }
    // if (w->nic > 0.5) {
    //     for (k = 0; k < w->nic; k++) {
    //         g[k] = 0;
    //     }
    // }
    free_cost(obm);
    return 0;
}

void BFGSudpate
(
    SOLNPWork* w,
    SUBNPWork* w_sub,
    SOLNPSettings* stgs,
    solnp_float* g,
    solnp_float* yg,
    solnp_float* sx
) 
{   
    SOLNP(add_scaled_array)(yg, g, w_sub->J->npic, -1.0);
    SOLNP(scale_array)(yg, -1.0, w_sub->J->npic);
    SOLNP(add_scaled_array)(sx, w->p, w_sub->J->npic, -1.0);
    SOLNP(scale_array)(sx, -1.0, w_sub->J->npic);
    solnp_float sc[2];
    solnp_float* temp = (solnp_float*)solnp_malloc(w_sub->J->npic*sizeof(solnp_float));
    SOLNP(Ax)(temp, w->h, sx, w_sub->J->npic, w_sub->J->npic);
    sc[0] = SOLNP(dot)(sx,temp, w_sub->J->npic);
    sc[1] = SOLNP(dot)(sx, yg, w_sub->J->npic);
    if (sc[0] * sc[1] > 0) {
        memcpy(sx, temp, w_sub->J->npic * sizeof(solnp_float));
        // two rank 1 updates
        SOLNP(rank1update)(w_sub->J->npic, w->h, -1 / sc[0], sx);
        SOLNP(rank1update)(w_sub->J->npic, w->h, 1 / sc[1], yg);
    }
    solnp_free(temp);
}

void solve_qp_osqp
(
    SOLNPWork* w,
    SUBNPWork* w_sub,
    solnp_float* p0,
    solnp_float* y,
    solnp_float* g
) {
    solnp_int i;

    solnp_float* u = (solnp_float*)solnp_malloc(w_sub->J->npic * sizeof(solnp_float));
    solnp_float* c = (solnp_float*)solnp_calloc(w_sub->J->npic, sizeof(solnp_float));


    solnp_float* lb = (solnp_float*)solnp_malloc(w_sub->J->npic * sizeof(solnp_float));
    solnp_float* ub = (solnp_float*)solnp_malloc(w_sub->J->npic * sizeof(solnp_float));
    for (i = 0; i < w_sub->J->npic; i++) {
        lb[i] = -100. / SQRTF(w_sub->J->npic);
        ub[i] = 100. / SQRTF(w_sub->J->npic);
    }
    for (i = 0; i < w_sub->J->npic; i++) {
        if (i < w->nic) {
            lb[i] = MAX(lb[i], w->pb->il[i] - w->p[i]);
            ub[i] = MIN(ub[i], w->pb->iu[i] - w->p[i]);
        }
        else {
            lb[i] = MAX(lb[i], w->pb->pl[i - w->nic] - w->p[i]);
            ub[i] = MIN(ub[i], w->pb->pu[i - w->nic] - w->p[i]);
        }
    }
    qpsolver(w_sub->J->npic, w->h, g, w->nc, w_sub->J->npic, w_sub->J->a, c, lb, ub, u, y);
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
    SOLNP(add_scaled_array)(u, w->p, w_sub->J->npic, 1.0);
    memcpy(p0, u, w_sub->J->npic * sizeof(solnp_float));
    proj(p0, w);
    solnp_free(c);
    solnp_free(u);

    return;
}
void solve_qp_aff
(
    SOLNPWork* w,
    SUBNPWork* w_sub,
    solnp_float* p0,
    solnp_float* y,
    solnp_float* g
) {

    SOLNP(timer) qpsolver_timer;
    SOLNP(tic)(&qpsolver_timer);


    solnp_float go = -1.0;
    solnp_int i;
    solnp_float* dx = (solnp_float*)solnp_malloc(w_sub->J->npic * sizeof(solnp_float));
    //calculate the ball
    for (i = 0; i < w_sub->J->npic; i++) {
        dx[i] = 0.01;
    }
    if (w->pb->Ipb[1] > 0.5) {
        solnp_float m = 0.01; //minimum of dx[i]
        for (i = 0; i < w_sub->mm; i++) {
            if (i < w->nic) {
                dx[i] = (1 / (MIN(w->p[i] - w->pb->il[i], w->pb->iu[i] - w->p[i]) + EPS));
            }
            else {
                dx[i] = (1 / (MIN(w->p[i] - w->pb->pl[i - w->nic], w->pb->pu[i - w->nic] - w->p[i]) + EPS));
            }

            if (dx[i] < m) {
                m = dx[i];
            }
        }
        if (w->pb->Ipb[0] <= 0) {
            m = MIN(m, 0.01);
            for (i = w_sub->mm; i < w_sub->J->npic; i++) {
                dx[i] = m;
            }
        }
    }
    w->mu /= 10;  //mu is l in subnp.m, lagrange multiplier
    //solnp_float* c = (solnp_float *)solnp_malloc(w_sub->J->npic * w_sub->J->npic * sizeof(solnp_float));
    solnp_float* u = (solnp_float*)solnp_malloc(w_sub->J->npic * sizeof(solnp_float));
    while (go <= 0 ) {
        //SOLNP(chol)(w_sub->J->npic,w->h,dx,c,w->mu);
        //c = inv(c);
        //yg = c' *g
        solnp_float* L = (solnp_float*)solnp_malloc(w_sub->J->npic * w_sub->J->npic * sizeof(solnp_float));
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
        if (w->nc > 0) {
            solnp_float* inv_H_aT = (solnp_float*)solnp_malloc(w_sub->J->npic * w->nc * sizeof(solnp_float));
            SOLNP(transpose)(w->nc, w_sub->J->npic, w_sub->J->a, inv_H_aT);
            SOLNP(solve_lin_sys)(w_sub->J->npic, w->nc, L, inv_H_aT);

            solnp_float* a_inv_H_aT = (solnp_float*)solnp_malloc((w->nc * w->nc) * sizeof(solnp_float));
            SOLNP(AB)(a_inv_H_aT, w_sub->J->a, inv_H_aT, w->nc, w_sub->J->npic, w->nc);

            SOLNP(chol)(w->nc, a_inv_H_aT);
            SOLNP(Ax)(y, w_sub->J->a, u, w->nc, w_sub->J->npic); // au = a  H^-1 u
            SOLNP(scale_array)(y, -1.0, w->nc);

            SOLNP(solve_lin_sys)(w->nc, 1, a_inv_H_aT, y); // a_inv_H_aT * y = a*H^-1*u

            solnp_float* inv_H_aT_y = (solnp_float*)solnp_malloc(w_sub->J->npic * sizeof(solnp_float));
            SOLNP(Ax)(inv_H_aT_y, inv_H_aT, y, w_sub->J->npic, w->nc);
            SOLNP(add_scaled_array)(u, inv_H_aT_y, w_sub->J->npic, 1.0);

            solnp_free(inv_H_aT);
            solnp_free(a_inv_H_aT);
            solnp_free(inv_H_aT_y);
        }
        solnp_free(L);

        SOLNP(add_scaled_array)(u, w->p, w_sub->J->npic, 1.0);
        memcpy(p0, u, w_sub->J->npic * sizeof(solnp_float));
        if isnan(p0[0]) {
            w->exit = 2;
        }

        if (w->pb->Ipb[1] <= 0) {
            go = 1;
        }
        else {
            go = INFINITY;
            for (i = 0; i < w_sub->mm; i++) {
                if (i < w->nic && go > MIN(p0[i] - w->pb->il[i], w->pb->iu[i] - p0[i])) {
                    go = MIN(p0[i] - w->pb->il[i], w->pb->iu[i] - p0[i]);
                }
                else if (i >= w->nic && go > MIN(p0[i] - w->pb->pl[i - w->nic], w->pb->pu[i - w->nic] - p0[i])) {
                    go = MIN(p0[i] - w->pb->pl[i - w->nic], w->pb->pu[i - w->nic] - p0[i]);
                }
            }
            w->mu = 3 * w->mu;
        }
    }
    solnp_free(u);
    solnp_free(dx);

    qpsolver_time += SOLNP(tocq)(&qpsolver_timer) / 1e3;

    return;
}
solnp_int linesearch
(
    SUBNPWork* w_sub,
    SOLNPWork* w,
    SOLNPSettings* stgs,
    solnp_float* j,
    solnp_float* reduce,
    solnp_float* p0,
    solnp_float* sx
)
{
    solnp_int i,tag = 0;
    solnp_float go = 1.0;
    solnp_float obm,obn,alm_cand;
    SOLNPCost* ob1 = init_cost(w->nec, w->nic);
    SOLNPCost* ob2 = init_cost(w->nec, w->nic);
    SOLNPCost* ob3 = init_cost(w->nec, w->nic);
    solnp_float** pt = (solnp_float**)solnp_malloc(3 * sizeof(solnp_float*));
    for (i = 0; i < 3; i++) {
        pt[i] = (solnp_float*)solnp_malloc(w_sub->J->npic * sizeof(solnp_float));
    }
    solnp_float sob[3];
    copySOLNPCost(ob1, w->ob);
    copySOLNPCost(ob2, w->ob);
    calculate_scaled_cost(ob3, p0, w_sub->scale,stgs, w);
    memcpy(pt[0], w->p, w_sub->J->npic * sizeof(solnp_float));
    memcpy(pt[1], w->p, w_sub->J->npic * sizeof(solnp_float));
    memcpy(pt[2], p0, w_sub->J->npic * sizeof(solnp_float));
    sob[0] = *j;
    sob[1] = *j;
    sob[2] = calculate_ALM(ob3, stgs, p0, w, w_sub);// need to modify
    w_sub->alp[0] = 0.0;
    w_sub->alp[2] = 1.0;
    solnp_int lstime = 0;
    while (go > stgs->tol && lstime < stgs->ls_time) {
        if (SOLNP(min)(sob,3) < *j || w->exit == 1) {
            break;
        }
        lstime++;
        memcpy(pt[1], p0, w_sub->J->npic * sizeof(solnp_float));
        w_sub->alp[1] = (w_sub->alp[0] + w_sub->alp[2]) / 2.0;
        SOLNP(scale_array)(pt[1], w_sub->alp[1], w_sub->npic);
        SOLNP(add_scaled_array)(pt[1], w->p, w_sub->npic, (1 - w_sub->alp[1]));
        calculate_scaled_cost(ob2, pt[1], w_sub->scale,stgs, w);
        sob[1] = calculate_ALM(ob2,stgs,pt[1],w,w_sub);
        obm = SOLNP(max)(sob,3);
        if (obm < *j) {
            obn = SOLNP(min)(sob, 3);
            go = stgs->tol * (obm - obn) / (*j - obm);
        }
        if (sob[1] >= sob[0]) {
            sob[2] = sob[1]; 
            copySOLNPCost(ob3, ob2);
            w_sub->alp[2] = w_sub->alp[1];
            memcpy(pt[2], pt[1], w_sub->J->npic * sizeof(solnp_float));
        }
        else if (sob[0] <= sob[2]) {
            sob[2] = sob[1];
            copySOLNPCost(ob3, ob2);
            w_sub->alp[2] = w_sub->alp[1];
            memcpy(pt[2], pt[1], w_sub->J->npic * sizeof(solnp_float));
        }
        else {
            sob[0] = sob[1];
            copySOLNPCost(ob1, ob2);
            w_sub->alp[0] = w_sub->alp[1];
            memcpy(pt[0], pt[1], w_sub->J->npic * sizeof(solnp_float));
        }
        if (go >= stgs->tol) {
            go = w_sub->alp[2] - w_sub->alp[0];
        }
    }

    memcpy(sx, w->p, w_sub->J->npic * sizeof(solnp_float));
    w_sub->ch = 1;
    obn = SOLNP(min)(sob, 3);
    if (w_sub->ob_cand->obj != INFINITY) {
        alm_cand = calculate_ALM(w_sub->ob_cand, stgs, w_sub->p_cand, w, w_sub);
        obn = MIN(obn, alm_cand);
    }
    if (*j <= obn) {
        tag = 1;
    }
    *reduce = (*j - obn) / MAX(1 , fabs(*j));
    if (stgs->noise) {
        if (*reduce > 30 * stgs->delta) {
            stgs->delta = stgs->delta * stgs->k_i;
        }
        if (*reduce < MAX(stgs->tol, 10 * stgs->delta)) {
            stgs->delta = MAX(stgs->delta_end, stgs->delta/stgs->k_r);
            tag = 1;
        }
    }
    else {
        if (*reduce < stgs->tol) {
            tag = 1;
        }
    }
    if (sob[0] < sob[1]) {
        *j = sob[0]; 
        memcpy(w->p, pt[0], w_sub->J->npic * sizeof(solnp_float));
        copySOLNPCost(w->ob, ob1);
    }
    else if (sob[2] < sob[1]) {
        *j = sob[2];
        memcpy(w->p, pt[2], w_sub->J->npic * sizeof(solnp_float));
        copySOLNPCost(w->ob, ob3);
    }
    else {
        *j = sob[1];
        memcpy(w->p, pt[1], w_sub->J->npic * sizeof(solnp_float));
        copySOLNPCost(w->ob, ob2);
    }
    if (w->ob->obj > w_sub->ob_cand->obj) {
        *j = alm_cand;
        w_sub->ch = 1;
        copySOLNPCost(w->ob, w_sub->ob_cand);
        memcpy(w->p, w_sub->p_cand, w_sub->npic * sizeof(solnp_float));
    }
    free_cost(ob1);
    free_cost(ob2);
    free_cost(ob3);
    for (i = 0; i < 3; i++) {
        solnp_free(pt[i]);
    }
    solnp_free(pt);
    return tag;
}

solnp_int SUBNP(solve)
(
    SUBNPWork *w_sub,
    SOLNPWork *w,
    SOLNPSettings *stgs
)
{
    solnp_int minit = 0;
    /* scale procedure */
    update_work_subnp(w_sub, w, stgs);
    w->alm_crit = calculate_almcrit_iq(w, w_sub->scale, w->p);
    /* calculate Jacobian matrices */
    if (w->nc > 0.5) {
        calculate_Jacob(w_sub, w, stgs);
        if (w->exit == 1) {
            unscale(w_sub, w, stgs);
            return 0;
        }
        // record constraint value into last column of A
        memcpy(&w_sub->J->a[w_sub->npic * w_sub->nc], w->ob->ec, w->nec * sizeof(solnp_float));
        if (w->nic > 0.5) {
            memcpy(&w_sub->J->a[w->nec + w_sub->npic * w_sub->nc], w->ob->ic, w->nic * sizeof(solnp_float));
            SOLNP(add_scaled_array)(&w_sub->J->a[w->nec + w_sub->npic * w_sub->nc], w->p, w->nic, -1.0);
        }
        SOLNP(Ax)(w_sub->b, w_sub->J->a, w->p, w_sub->nc, w_sub->npic);
        SOLNP(add_scaled_array)(w_sub->b, &w_sub->J->a[w_sub->npic * w_sub->nc], w->nc, -1.0);
    }   w->alm_crit = SQRTF(w->alm_crit);
    // corrspond to subnp line 73
    /* find interior (near-)feasible solution */
     if (w->nc > 0.5) {
         if (stgs->qpsolver == 1) {
             find_int_feas_sol_aff(w_sub, w, stgs);
         }
         else {
             find_int_feas_sol_osqp(w_sub, w, stgs);
         }
    }
    // recalculate cost if find new solution p
    // subnp.m line 139
    solnp_float j;
    if (w_sub->ch>0) {
        calculate_scaled_cost(w->ob, w->p, w_sub->scale,stgs, w);
    }
    j = calculate_ALM(w->ob, stgs, w->p, w, w_sub);
    if (w->exit == 1) {
        unscale(w_sub, w, stgs);
        return 0;
    }
    // solve su bproblem QP and BFGS update Hessian 
    // TODO: subnp.m line start at 152
    solnp_int tag;
    solnp_float* yg = (solnp_float*)solnp_malloc(w_sub->J->npic * sizeof(solnp_float));
    solnp_float* sx = (solnp_float*)solnp_malloc(w_sub->J->npic * sizeof(solnp_float)); 
    solnp_float* p0 = (solnp_float*)solnp_malloc(w_sub->J->npic * sizeof(solnp_float));
    solnp_float* y = (solnp_float*)solnp_calloc(w->nc, sizeof(solnp_float));
    solnp_float reduce;
    memcpy(p0, w->p, w_sub->J->npic* sizeof(solnp_float));
    
    while (minit < stgs->min_iter) {
        if (minit > 0) {
            w_sub->ob_cand->obj = INFINITY;
        }
        minit++;
        if (w_sub->ch > 0) {
            // calculate the gradient
            w->alm_crit = calculate_almcrit_iq(w, w_sub->scale, w->p);
           calculate_ALMgradient(w_sub, w, stgs, w_sub->J->g, j);
           if (w->exit == 1) {
               break;
           }
           w->alm_crit = SQRTF(w->alm_crit);
        }
        if (minit > 1) {
            // BFGS update
           BFGSudpate(w, w_sub, stgs, w_sub->J->g, yg, sx);
        }
        
       //solve QP subproblem
        if (stgs->qpsolver == 1) {
            solve_qp_aff(w, w_sub, p0, y, w_sub->J->g);
        }
        else {
            solve_qp_osqp(w, w_sub, p0,y, w_sub->J->g);
        }
        //Line Search 
        tag = linesearch(w_sub, w, stgs, &j, &reduce, p0, sx);
        if (w->exit >= 1) {
            break;
        }
        memcpy(yg, w_sub->J->g, w_sub->J->npic * sizeof(solnp_float));
        if (tag) {
            break;
        }
    }
    memcpy(w->l, y, w->nc * sizeof(solnp_float));
    // TODO:
    /* unscale and record information */
    unscale(w_sub, w, stgs);
    solnp_free(y);
    solnp_free(yg);
    solnp_free(sx);
    solnp_free(p0);  
    /*
    if (reduce > MAX(stgs->tol,10*stgs->delta) && w->exit == 0) {
        printf("SOLNP+--> Minor optimization routine did not converge \n         in the specified number of minor iterations.\n         You may need to increase the number of minor iterations. \n ");
    }*/
    return 0;
}

solnp_int free_Jacob(SUBNPJacob *J) {
    if (J) {
        if (J->g) {
            solnp_free(J->g);
        }
        if (J->a) {
            solnp_free(J->a);
        }
        solnp_free(J);
    }
    return 0;
}

solnp_int free_work_subnp(SUBNPWork *w_sub)
{
    if (w_sub){
        if (w_sub->alp) {
            solnp_free(w_sub->alp);
        }
        if (w_sub->scale) {
            solnp_free(w_sub->scale);
        }
        if (w_sub->J) {
            free_Jacob(w_sub->J);
        }
        if (w_sub->b) {
            solnp_free(w_sub->b);
        }
        if (w_sub->p_cand) {
            solnp_free(w_sub->p_cand);
        }
        if (w_sub->ob_cand) {
            free_cost(w_sub->ob_cand);
        }
        solnp_free(w_sub);
    }

    return 0;
}

solnp_int SUBNP(finish)
(
    SUBNPWork *w_sub,
    SOLNPWork *w,
    SOLNPSettings *stgs
)
{   
    
    if (w_sub){
        free_work_subnp(w_sub);
    }
    return 0;
}


solnp_int subnp_qp(SOLNPWork *w, SOLNPSettings *stgs, SOLNPInfo *info)
{
    SUBNPWork *w_sub = SUBNP(init)(w, stgs);

    if(w_sub) 
    {
        SUBNP(solve)(w_sub, w, stgs);
        
        SUBNP(finish)(w_sub, w, stgs);
    }

    info->qpsolver_time = qpsolver_time;
    qpsolver_time = 0;

    return 0;
}

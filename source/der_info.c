#include "der_info.h"
#include "linalg.h"
#include "linsys.h"
#include "qp_solver.h"

solnp_float* Est_noise(
    solnp_int nf,
    solnp_float* fval
) {
    // This subroutine implement the ECnoise program proposed by Argonne National Laboratory
    // Jorge More' and Stefan Wild. November 2009. 
    // The input: nf: number of points 
    // fval: array of funciton values
    // Output: fnoise_inform: the first entry is estimated noise. The second entry is infomration of the output.
    solnp_float* fnoise_inform = (solnp_float*) solnp_malloc(2*sizeof(solnp_float));
    solnp_float* level = (solnp_float*)solnp_calloc(nf - 1, sizeof(solnp_float));
    solnp_float* dsgn = (solnp_float*)solnp_calloc(nf - 1, sizeof(solnp_float));
    fnoise_inform[0] = 0;
    solnp_float gamma = 1.;
    solnp_float fmin = solnp_min(fval, nf);
    solnp_float fmax = solnp_max(fval, nf);
    if ((fmax-fmin)/(MAX(ABS(fmin),ABS(fmax))) > .1)
    {
        fnoise_inform[1] = 3;
        solnp_free(level);
        solnp_free(dsgn);
        return fnoise_inform;
    }
    for (solnp_int j = 0; j < nf-1; j++) {
        solnp_int count_zero = 0;
        for (solnp_int i = 0; i < nf - j; i++) {
            fval[i] = fval[i + 1] - fval[i];        
        }
        for (solnp_int i = 0; i < nf - 1; i++) {
            if (fval[i] == 0) {
                count_zero++;
            }
        }

        if (j == 1 && count_zero >= nf/2) {
            fnoise_inform[1] = 2;
            solnp_free(level);
            solnp_free(dsgn);
            return fnoise_inform;
        }
        gamma *= 0.5 * ((j+1.) / (2 * j + 1.));
        
        //Compute the estimates for the noise level.
        level[j] = SQRTF(gamma * SOLNP(norm_sq)(fval, nf - j) / (nf - j));

        //Determine differences in sign.
        solnp_float emin = SOLNP(min)(fval, nf - j);
        solnp_float emax = SOLNP(max)(fval, nf - j);
        if (emin * emax < 0) {
            dsgn[j] = 1;
        }
    }

    for (solnp_int k = 0; k < nf - 3; k++) {
        solnp_float emin = SOLNP(min)(&level[k], 3);
        solnp_float emax = SOLNP(max)(&level[k], 3);
        if (emax <= 4 * emin && dsgn[k]) {
            fnoise_inform[0] = level[k];
            fnoise_inform[1] = 1;
            solnp_free(level);
            solnp_free(dsgn);
            return fnoise_inform;
        }
    }

    fnoise_inform[1] = 3;
    solnp_free(level);
    solnp_free(dsgn);
    return fnoise_inform;
}

solnp_float* Est_second_diff_sub(
    solnp_float fval,
    solnp_float delta,
    solnp_float fnoise,
    solnp_float* p,
    solnp_float* d,
    SUBNPWork* w_sub,
    SOLNPWork* w,
    SOLNPSettings* stgs

) {
    solnp_float tau1 = 100;
    solnp_float tau2 = 0.1;
    solnp_float* res = (solnp_float*)solnp_malloc(2*sizeof(solnp_float));
    solnp_float delta_h;
    SOLNPCost** obm = (SOLNPCost**)solnp_malloc(1 * sizeof(SOLNPCost*));
    obm[0] = init_cost(w->nec, w->nic);
    solnp_float* ptemp = (solnp_float*)solnp_malloc(w_sub->J->npic * sizeof(solnp_float));
    memcpy(ptemp, p, (w->nic+w->n) * sizeof(solnp_float));
    solnp_add_scaled_array(ptemp, d, w->n + w->nic,delta);
    calculate_scaled_cost(obm, ptemp, w_sub->scale, stgs, w, 1);
    solnp_float alm_f = obm[0]->obj;

    if (w_sub->nc > 0)
    {
        alm_f = calculate_ALM(obm[0], stgs, ptemp, w, w_sub);
    }

    solnp_add_scaled_array(ptemp, d, w->n + w->nic, -2*delta);
    calculate_scaled_cost(obm, ptemp, w_sub->scale, stgs, w, 1);
    solnp_float alm_b = obm[0]->obj;

    if (w_sub->nc > 0)
    {
        alm_b = calculate_ALM(obm[0], stgs, ptemp, w, w_sub);
    }

    delta_h = ABS(alm_f + alm_b - fval * 2);
    res[0] = delta_h / (delta * delta);
    if (delta_h >= fnoise * tau1 && ABS(alm_f - fval) <= tau2 * MAX(ABS(fval), ABS(alm_f)) \
        && ABS(alm_b - fval) <= tau2 * MAX(ABS(fval), ABS(alm_b))) {
        res[1] = 1;
    }
    else {
        res[1] = 0;
    }

    solnp_free(ptemp);
    free_cost(obm[0]);
    solnp_free(obm);
    return res;
}

solnp_float Est_second_diff(
    solnp_float fnoise,
    solnp_float fval,
    solnp_float* p,
    solnp_float* d,
    SUBNPWork* w_sub,
    SOLNPWork* w,
    SOLNPSettings* stgs
) {
    // This subroutine estimate the second-order derivative of a function  

    solnp_float ha = pow(fnoise, 0.25);
    solnp_float mu;
    solnp_float* mua = Est_second_diff_sub(fval, ha, fnoise, p, d, w_sub, w, stgs);
    if (mua[1] == 1){
        mu = mua[0];
        solnp_free(mua);
        return mu;
    }
    solnp_float hb = pow(fnoise/mua[0], 0.25);
    solnp_float* mub = Est_second_diff_sub(fval, hb, fnoise, p, d, w_sub, w, stgs);
    if (mub[1] == 1 || ABS(mua[0]-mub[0]) <= 0.5*mub[0] ) {
        mu = mub[0];
        solnp_free(mua);
        solnp_free(mub);
        return mu;
    }

    //The noise is too large.
    solnp_free(mua);
    solnp_free(mub);
    return 1;
}

solnp_float calculate_delta(
    solnp_int nf,
    solnp_float fval,
    solnp_float* p,
    SUBNPWork* w_sub,
    SOLNPWork* w,
    SOLNPSettings* stgs
) {
    // This subroutine is used to calculate the step size of finite-difference method
    solnp_float* fvals = (solnp_float*)solnp_malloc(nf * sizeof(solnp_float));
    SOLNPCost** obm = (SOLNPCost**)solnp_malloc(1 * sizeof(SOLNPCost*));
    obm[0] = init_cost(w->nec, w->nic);
    solnp_float h_default = stgs->h;
    solnp_float* d = (solnp_float*)solnp_malloc((w->n + w->nic) * sizeof(solnp_float));
    Uniform_sphere(d, w->n, 1);

    for (solnp_int i = 0; i < (nf - 1) / 2; i++) {
        solnp_float* ptemp = (solnp_float*)solnp_malloc(w_sub->J->npic * sizeof(solnp_float));
        memcpy(ptemp, p, (w->nic + w->n) * sizeof(solnp_float));
        solnp_add_scaled_array(ptemp, d, w->n, h_default);
        calculate_scaled_cost(obm, ptemp, w_sub->scale, stgs, w, 1);
        fvals[i] = obm[0]->obj;

        if (w_sub->nc > 0)
        {
            fvals[i] = calculate_ALM(obm[0], stgs, ptemp, w, w_sub);
        }

        solnp_add_scaled_array(ptemp, d, w->n, -2*h_default);
        calculate_scaled_cost(obm, ptemp, w_sub->scale, stgs, w, 1);
        fvals[nf-i-1] = obm[0]->obj;

        if (w_sub->nc > 0)
        {
            fvals[nf - i - 1] = calculate_ALM(obm[0], stgs, ptemp, w, w_sub);
        }
        solnp_free(ptemp);
    }
    fvals[(nf - 1) / 2] = fval;
    // Estimate noise
    solnp_float* res_ecnoise = Est_noise(nf, fvals);
    if (res_ecnoise[1] != 1) {
        if (res_ecnoise[1] == 2) {
            stgs->h *= 10;
        }
        else {
            stgs->h /= 10;
        }

        solnp_free(res_ecnoise);
        solnp_free(d);
        solnp_free(fvals);
        free_cost(obm[0]);
        return -1.;
    }
    // Estimate second-order derivative
    solnp_float v2 = Est_second_diff(res_ecnoise[0], fval, p, d, w_sub, w, stgs);
    solnp_float delta = pow(8, 0.25) * pow(res_ecnoise[0] / v2, 0.5);

    solnp_free(res_ecnoise);
    solnp_free(d);
    solnp_free(fvals);
    free_cost(obm[0]);
    solnp_free(obm);
    return delta;
}

solnp_float calculate_infeas_scaledob(
    SOLNPCost *ob,
    SOLNPWork *w,
    solnp_float *x,
    solnp_float *scale)
{
    solnp_float con = 0;
    solnp_int i;
    // calculate constraints
    if (w->nic > 0.5)
    {
        for (i = 0; i < w->nic; i++)
        {
            if (ob->ic[i] < w->pb->il[i])
            {
                con += (w->pb->il[i] - ob->ic[i]) * (w->pb->il[i] - ob->ic[i]) * scale[1 + w->nec + i] * scale[1 + w->nec + i];
            }
            if (ob->ic[i] > w->pb->iu[i])
            {
                con += (ob->ic[i] - w->pb->iu[i]) * (ob->ic[i] - w->pb->iu[i]) * scale[1 + w->nec + i] * scale[1 + w->nec + i];
            }
        }
    }
    if (w->nec > 0.5)
    {
        for (i = 0; i < w->nec; i++)
        {
            con += ob->ec[i] * ob->ec[i] * scale[1 + i] * scale[1 + i];
        }
    }
    for (i = 0; i < w->pb->n; i++)
    {
        if (x[i] > w->pb->pu[i])
        {
            con += (x[i] - w->pb->pu[i]) * (x[i] - w->pb->pu[i]) * scale[1 + w->nec + w->nic + i] * scale[1 + w->nec + w->nic + i];
        }
        if (x[i] < w->pb->pl[i])
        {
            con += (x[i] - w->pb->pl[i]) * (x[i] - w->pb->pl[i]) * scale[1 + w->nec + w->nic + i] * scale[1 + w->nec + w->nic + i];
        }
    }
    con = SQRTF(con);
    return con;
}

solnp_float calculate_infeas_scaledob_l1(
    SOLNPCost *ob,
    SOLNPWork *w,
    solnp_float *x)
{
    solnp_float con = 0;
    solnp_int i;
    // calculate constraints
    if (w->nic > 0.5)
    {
        for (i = 0; i < w->nic; i++)
        {
            if (ob->ic[i] < w->pb->il[i])
            {
                con += (w->pb->il[i] - ob->ic[i]) * (w->pb->il[i] - ob->ic[i]);
            }
            if (ob->ic[i] > w->pb->iu[i])
            {
                con += (ob->ic[i] - w->pb->iu[i]) * (ob->ic[i] - w->pb->iu[i]);
            }
        }
    }
    if (w->nec > 0.5)
    {
        for (i = 0; i < w->nec; i++)
        {
            con += (ob->ec[i] + x[i + w->n - 2 * w->nec] - x[i + w->n - w->nec]) * (ob->ec[i] + x[i + w->n - 2 * w->nec] - x[i + w->n - w->nec]);
        }
    }
    for (i = 0; i < w->n - 2 * w->nec; i++)
    {
        if (x[i] > w->pb->pu[i])
        {
            con += (x[i] - w->pb->pu[i]) * (x[i] - w->pb->pu[i]);
        }
        if (x[i] < w->pb->pl[i])
        {
            con += (x[i] - w->pb->pl[i]) * (x[i] - w->pb->pl[i]);
        }
    }
    con = SQRTF(con);
    return con;
}

solnp_float calculate_infeas_unscaleob(
    SOLNPCost *ob,
    SOLNPWork *w,
    solnp_float *x,
    solnp_float *scale)
{
    solnp_float con = 0;
    solnp_int i;
    // calculate constraints
    if (w->nic > 0.5)
    {
        for (i = 0; i < w->nic; i++)
        {
            if (ob->ic[i] < w->pb->il[i] * scale[1 + w->nec + i])
            {
                con += (w->pb->il[i] * scale[1 + w->nec + i] - ob->ic[i]) * (w->pb->il[i] * scale[1 + w->nec + i] - ob->ic[i]);
            }
            if (ob->ic[i] > w->pb->iu[i] * scale[1 + w->nec + i])
            {
                con += (ob->ic[i] - w->pb->iu[i] * scale[1 + w->nec + i]) * (ob->ic[i] - w->pb->iu[i] * scale[1 + w->nec + i]);
            }
        }
    }
    if (w->nec > 0.5)
    {
        for (i = 0; i < w->nec; i++)
        {
            con += ob->ec[i] * ob->ec[i];
        }
    }
    for (i = 0; i < w->pb->n; i++)
    {
        if (x[i] > w->pb->pu[i] * scale[1 + w->nec + w->nic + i])
        {
            con += (x[i] - w->pb->pu[i] * scale[1 + w->nec + w->nic + i]) * (x[i] - w->pb->pu[i] * scale[1 + w->nec + w->nic + i]);
        }
        if (x[i] < w->pb->pl[i] * scale[1 + w->nec + w->nic + i])
        {
            con += (x[i] - w->pb->pl[i] * scale[1 + w->nec + w->nic + i]) * (x[i] - w->pb->pl[i] * scale[1 + w->nec + w->nic + i]);
        }
    }
    con = SQRTF(con);
    return con;
}
solnp_float calculate_almcrit_iq(
    SOLNPWork *w,
    solnp_float *scale,
    solnp_float *p)
{
    solnp_float almcrit = 0;
    solnp_float temp;
    solnp_int i;
    for (i = 0; i < w->nic; i++)
    {
        temp = p[i] * scale[w->nec + 1 + i] - scale[0] * w->l[w->nec + i] / scale[w->nec + 1 + i];
        if (temp > w->pb->iu[i] * scale[w->nec + 1 + i])
        {
            temp = w->pb->iu[i] * scale[w->nec + 1 + i];
        }
        else if (temp < w->pb->il[i] * scale[w->nec + 1 + i])
        {
            temp = w->pb->il[i] * scale[w->nec + 1 + i];
        }
        almcrit += (temp - p[i] * scale[w->nec + 1 + i]) * (temp - p[i] * scale[w->nec + 1 + i]);
    }
    return almcrit;
}
void calculate_scaled_cost(
    SOLNPCost **ob,
    solnp_float *p,
    const solnp_float *scale,
    SOLNPSettings *stgs,
    SOLNPWork *w,
    solnp_int nfeval)
{
    solnp_int i, j;
    solnp_int len = stgs->rescue ? w->n - 2 * w->nec : w->n;
    w->pb->n = len;
    solnp_float *x = (solnp_float *)solnp_malloc(len * nfeval * sizeof(solnp_float));
    for (j = 0; j < nfeval; j++)
    {
        for (i = 0; i < len; i++)
        {
            x[i + j * len] = p[i + j * w->n] * scale[1 + w->nc + i];
        }
    }
    // calculate rescaled cost
    // (*w->ob->cost)(ob, x, w->n, 0);
    w->ob->cost(ob, x, len, nfeval, 0);
    w->count_cost += nfeval;

    for (j = 0; j < nfeval; j++)
    {
        solnp_float con = calculate_infeas_unscaleob(ob[j], w, &x[j * len], scale);
        if (con <= stgs->tol_con && w->bestobj > ob[j]->obj)
        {
            w->bestcon = con;
            w->bestobj = ob[j]->obj;
            memcpy(w->bestp, &x[j * len], len * sizeof(solnp_float));
            memcpy(w->bestl, w->l, MAX(w->nc, 1) * sizeof(solnp_float));
            copySOLNPCost(w->bestob, ob[j]);
        }

        if (con < w->best_fea_con)
        {
            w->best_fea_con = con;
            memcpy(w->best_fea_p, &x[j * len], len * sizeof(solnp_float));
            memcpy(w->best_fea_l, w->l, MAX(w->nc, 1) * sizeof(solnp_float));
            copySOLNPCost(w->best_fea_ob, ob[j]);
        }
        if (stgs->rescue && w->nec)
        {
            // ob[j]->obj = SOLNP(norm_sq)(w->p + w->nic + w->n - w->nec, w->nec);
            //  Modify equality constraints and objective
            for (i = 0; i < w->nec; i++)
            {
                ob[j]->obj += w->pen_l1 * (p[i + w->n - 2 * w->nec + j * w->n] * scale[1 + w->nc + i + w->n - 2 * w->nec] + p[i + w->n - w->nec + j * w->n] * scale[1 + w->nc + i + w->n - w->nec]);
            }

            for (i = 0; i < w->nec; i++)
            {
                ob[j]->ec[i] += -p[i + w->n - 2 * w->nec + j * w->n] * scale[1 + w->nc + i + w->n - 2 * w->nec] + p[i + w->n - w->nec + j * w->n] * scale[1 + w->nc + i + w->n - w->nec];
            }
        }
        rescale_ob(ob[j], scale);
    }
    if (w->count_cost >= stgs->maxfev)
    {
        w->exit = 1;
    }
    w->pb->n = w->n;
    // free x
    solnp_free(x);
}
/*
void calculate_scaled_cost_rescue
(
    SOLNPCost** ob,
    solnp_float* p,
    solnp_float* slack,
    const solnp_float* scale,
    SOLNPSettings* stgs,
    SOLNPWork* w,
    solnp_int nfeval
)
{
    solnp_int i, j;
    solnp_int len = w->n - w->nec;
    solnp_float* x = (solnp_float*)solnp_malloc(len * nfeval * sizeof(solnp_float));
    for (j = 0; j < nfeval; j++) {
        for (i = 0; i < len; i++) {
            x[i + j * len] = p[i + j * len] * scale[1 + w->nc + i];
        }
    }
    // calculate rescaled cost
    // (*w->ob->cost)(ob, x, len, 0);
    w->ob->cost(ob, x, len, nfeval, 0);
    w->count_cost += nfeval;

    for (j = 0; j < nfeval; j++) {
        if (stgs->rescue) {
            w->n -= w->nec;
        }

        solnp_float con = calculate_infeas_unscaleob(ob[j], w, &x[j * len], scale);
        if (stgs->rescue) {
            w->n += w->nec;
        }

        if (con <= stgs->tol_con && w->bestobj > ob[j]->obj) {
            w->bestcon = con;
            w->bestobj = ob[j]->obj;
            memcpy(w->bestp, &x[j * len], (len) * sizeof(solnp_float));
            for (i = 0; i < w->nc; i++) {
                w->bestl[i] = w->l[i] * scale[0] / scale[1 + i];
            }
        }
        if (con < w->best_fea_con) {
            w->best_fea_con = con;
            memcpy(w->best_fea_p, &x[j * len], (len) * sizeof(solnp_float));
            for (i = 0; i < w->nc; i++) {
                w->best_fea_l[i] = w->l[i] * scale[0] / scale[1 + i];
            }
        }
        if (stgs->rescue && w->nec) {
            //ob[j]->obj = SOLNP(norm_sq)(w->p + w->nic + w->n - w->nec, w->nec);
            ob[j]->obj = 0;
            for (i = 0; i < w->nec; i++) {
                ob[j]->obj += slack[i] * scale[1 + w->nc + w->n - w->nec + i]* slack[i] * scale[1 + w->nc + w->n - w->nec + i];
            }

            for (i = 0; i < w->nec; i++) {
                ob[j]->ec[i] -= slack[i] * scale[1 + w->nc + w->n - w->nec + i];
            }
        }
        rescale_ob(ob[j], scale);
    }
    if (w->count_cost >= stgs->maxfev) {
        w->exit = 1;
    }

    // free x
    solnp_free(x);
}*/

solnp_float line_search_merit(
    SOLNPCost *ob,
    SOLNPWork *w,
    SUBNPWork *w_sub,
    solnp_float *p)
{
    // Input: ob after modification
    // Reminder: The function is for rescue case! I will write the  general function later.
    // Output: L1 exact penlaty function value
    solnp_float result;
    solnp_int i;

    result = ob->obj * w_sub->scale[0];
    for (i = 0; i < w->nec; i++)
    {
        // recover objective
        result -= w->pen_l1 * (p[w->nic + w->n + i - 2 * w->nec] * w_sub->scale[1 + w->nc + i + w->n - 2 * w->nec] + p[w->nic + w->n + i - 1 * w->nec] * w_sub->scale[1 + w->nc + i + w->n - 1 * w->nec]);
        // equality penalty
        result += w->pen_l1 * ABS(ob->ec[i] + p[w->nic + w->n + i - 2 * w->nec] * w_sub->scale[1 + w->nc + i + w->n - 2 * w->nec] - p[w->nic + w->n + i - 1 * w->nec] * w_sub->scale[1 + w->nc + i + w->n - 1 * w->nec]) * w_sub->scale[1 + i];
    }
    for (i = 0; i < w->nic; i++)
    {
        // ineqality penalty
        if (ob->ic[i] < w->pb->il[i])
        {
            result += w->pen_l1 * (w->pb->il[i] - ob->ic[i]) * w_sub->scale[w->nec + 1 + i];
        }
        if (ob->ic[i] > w->pb->iu[i])
        {
            result += w->pen_l1 * (ob->ic[i] - w->pb->iu[i]) * w_sub->scale[w->nec + 1 + i];
        }
    }
    return result;
}

solnp_int calculate_scaled_grad(
    solnp_float *g,
    solnp_float *p,
    const solnp_float *scale,
    SOLNPSettings *stgs,
    SOLNPWork *w)
{
    solnp_float *x = (solnp_float *)solnp_malloc(w->n * sizeof(solnp_float));
    solnp_int i, j;

    for (i = 0; i < w->n; i++)
    {
        x[i] = p[i] * scale[1 + w->nc + i];
    }
    w->ob->grad(g, x, w->n, w->nc + 1, 0);
    for (i = 0; i < w->n; i++)
    {
        g[i] = g[i] / scale[0] * scale[1 + w->nc + i];
    }
    for (j = 0; j < w->nc; j++)
    {
        for (i = 0; i < w->n; i++)
        {
            g[w->n * (j + 1) + i] = g[w->n * (j + 1) + i] * scale[1 + w->nc + i] / scale[1 + j];
        }
    }
    w->count_grad++;
    solnp_free(x);
}

solnp_int calculate_scaled_grad_random(
    solnp_float *g,
    solnp_float *p,
    SOLNPCost *ob_p,
    const solnp_float *scale,
    SOLNPSettings *stgs,
    SOLNPWork *w)
{
    solnp_float *x = (solnp_float *)solnp_malloc(w->n * sizeof(solnp_float));
    solnp_float *x_temp = (solnp_float *)solnp_malloc(w->n * sizeof(solnp_float));
    solnp_float *diff = (solnp_float *)solnp_malloc(w->n * sizeof(solnp_float));
    SOLNPCost *ob = init_cost(w->nec, w->nic);
    SOLNPCost* ob_backward = init_cost(w->nec, w->nic);
    solnp_int *index = SOLNP_NULL;
    solnp_int i, j, k;

    // unscale x
    for (i = 0; i < w->n; i++)
    {
        x[i] = p[i] * scale[1 + w->nc + i];
        g[i] = 0;
    }
    // unscale ob
    unscale_ob(ob_p, scale);
    // Estimate the gradient by random sampling, Batchsize = stgs->batchsize

    solnp_int batch = MIN(stgs->batchsize, stgs->maxfev - w->count_cost);

    if (batch == 0)
    {
        w->exit = 1;
        return 0;
    }
    if (stgs->rs == 2)
    {
        index = (solnp_int *)solnp_calloc(sizeof(solnp_int), w->n);
    }
    for (i = 0; i < batch; i++)
    {
        if (stgs->rs == 1)
        {
            // Generate random Guassian Vector
            Uniform_sphere(diff, w->n, 1.);
            memcpy(x_temp, x, w->n * sizeof(solnp_float));
            SOLNP(add_scaled_array)
            (x_temp, diff, w->n, stgs->delta);
            // Estimate Gradient
            w->ob->cost(&ob, x_temp, w->n, 1, 0);
            w->count_cost += 1;

            if (stgs->cen_diff) {
                SOLNP(add_scaled_array)
                    (x_temp, diff, w->n, -2*stgs->delta);
                // Estimate Gradient
                w->ob->cost(&ob_backward, x_temp, w->n, 1, 0);
                w->count_cost += 1;
            }

            for (k = 0; k < w->n; k++)
            {
                if (stgs->cen_diff) {
                    // Use CENTRE difference to calculate gradient
                    g[k] += (ob->obj - ob_backward->obj) / (2*stgs->delta) * diff[k];
                    for (j = 1; j < w->nec + 1; j++)
                    {
                        g[k + j * w->n] += (ob->ec[j - 1] - ob_backward->ec[j - 1]) / (2 * stgs->delta) * diff[k];
                    }
                    for (j = w->nec + 1; j < w->nec + 1 + w->nic; j++)
                    {
                        g[k + j * w->n] = (ob->ic[j - w->nec - 1] - ob_backward->ic[j - w->nec - 1]) / (2 * stgs->delta) * diff[k];
                    }
                }
                else {
                    // use forward difference to calcualate gradient
                    g[k] += (ob->obj - ob_p->obj) / stgs->delta * diff[k];
                    for (j = 1; j < w->nec + 1; j++)
                    {
                        g[k + j * w->n] += (ob->ec[j - 1] - ob_p->ec[j - 1]) / stgs->delta * diff[k];
                    }
                    /*  if (isnan(g[k])) {
                          g[k] = g[k];
                      }*/
                    for (j = w->nec + 1; j < w->nec + 1 + w->nic; j++)
                    {
                        g[k + j * w->n] = (ob->ic[j - w->nec - 1] - ob_p->ic[j - w->nec - 1]) / stgs->delta * diff[k];
                    }
                }
            }
        }
        else if (stgs->rs == 2)
        {

            solnp_int r_ind = rand() % w->n;
            while (index[r_ind])
            {
                r_ind = rand() % w->n;
            }
            index[r_ind] = 1;
            memcpy(x_temp, x, w->n * sizeof(solnp_float));
            x_temp[r_ind] += stgs->delta;
            w->ob->cost(&ob, x_temp, w->n, 1, 0);
            w->count_cost += 1;
            if (stgs->cen_diff) {
                // Use CENTRE differnce to calculate the gradient
                x_temp[r_ind] -= 2 * stgs->delta;
                w->ob->cost(&ob_backward, x_temp, w->n, 1, 0);
                w->count_cost += 1;
                g[r_ind] = (ob->obj - ob_backward->obj) / (2*stgs->delta);
            }
            else {
                // Use forward difference to calculate the gradient
                g[r_ind] = (ob->obj - ob_p->obj) / stgs->delta;
            }
        }
    }

    // Average
    for (k = 0; k < w->n * (1 + w->nc); k++)
    {
        g[k] = g[k] / MAX(i, 1);
    }

    solnp_float g_norm = SOLNP(norm)(g, w->n);
    // SOLNP(set_as_scaled_array)(g, g, 1 / g_norm, w->n);

    // scale g
    for (i = 0; i < w->n; i++)
    {
        g[i] = g[i] / scale[0] * scale[1 + w->nc + i];
    }
    for (j = 0; j < w->nc; j++)
    {
        for (i = 0; i < w->n; i++)
        {
            g[w->n * (j + 1) + i] = g[w->n * (j + 1) + i] * scale[1 + w->nc + i] / scale[1 + j];
        }
    }

    // sacle ob
    rescale_ob(ob_p, scale);

    w->count_grad++;
    solnp_free(x);
    solnp_free(x_temp);
    solnp_free(diff);
    solnp_free(index);
    free_cost(ob);
    free_cost(ob_backward);
}

solnp_int calculate_scaled_hess(
    solnp_float *h,
    solnp_float *p,
    const solnp_float *scale,
    SOLNPSettings *stgs,
    SOLNPWork *w)
{
    solnp_float *x = (solnp_float *)solnp_malloc(w->n * sizeof(solnp_float));
    solnp_int i, j, k;

    for (i = 0; i < w->n; i++)
    {
        x[i] = p[i] * scale[1 + w->nc + i];
    }
    w->ob->hess(h, x, w->n, w->nc + 1, 0);

    // rescaled
    for (i = 0; i < w->n; i++)
    {
        for (j = 0; j < w->n; j++)
        {
            h[i + j * w->n] = h[i + j * w->n] / scale[0] * scale[1 + w->nc + i] * scale[1 + w->nc + j];
        }
    }
    for (k = 1; k < w->nc + 1; k++)
    {
        for (i = 0; i < w->n; i++)
        {
            for (j = 0; j < w->n; j++)
            {
                h[i + j * w->n + k * w->n * w->n] *= scale[1 + w->nc + i] * scale[1 + w->nc + j] / scale[k];
            }
        }
    }
    w->count_hess++;
    solnp_free(x);
}
void calculate_alm_criterion(
    SOLNPWork *w,
    SUBNPWork *w_sub,
    solnp_float *grad)
{
    // This function is used for calculate the gradient of ALM when provided gradient information
    solnp_int i;
    solnp_float *temp_alm = (solnp_float *)solnp_malloc(w_sub->npic * sizeof(solnp_float));
    solnp_float *gg_y = (solnp_float *)solnp_malloc(w->n * sizeof(solnp_float));
    solnp_float *temp_l = (solnp_float *)solnp_malloc(w->nic * sizeof(solnp_float));
    memcpy(temp_l, w->l + w->nec, sizeof(solnp_float) * w->nic);
    memcpy(temp_alm, w->p, w_sub->npic * sizeof(solnp_float));
    for (i = 0; i < w->nic; i++)
    {
        temp_l[i] *= w_sub->scale[0] / (w_sub->scale[1 + w->nec + i] * w_sub->scale[1 + w->nec + i]);
    }
    // calculate ALM stop criterion
    SOLNP(add_scaled_array)
    (temp_alm + w->nic, grad, w->n, -1.);
    if (w->nc)
    {
        SOLNP(Ax)
        (gg_y, grad + w->n, w->l, w->n, w->nc);
        SOLNP(add_scaled_array)
        (temp_alm + w->nic, gg_y, w->n, 1.);
        SOLNP(add_scaled_array)
        (temp_alm, temp_l, w->nic, -1.);
    }
    proj(temp_alm, w);
    SOLNP(add_scaled_array)
    (temp_alm, w->p, w_sub->npic, -1.);
    for (i = 0; i < w->n; i++)
    {
        temp_alm[i + w->nic] *= w_sub->scale[0] / w_sub->scale[1 + w->nc + i];
    }
    for (i = 0; i < w->nic; i++)
    {
        temp_alm[i] *= w_sub->scale[1 + w->nec + i];
    }
    w->alm_crit = SOLNP(norm)(temp_alm, w_sub->npic);
    solnp_free(temp_alm);
    solnp_free(gg_y);
    solnp_free(temp_l);
}

solnp_int calculate_Jacob_zero(
    SUBNPWork *w_sub,
    SOLNPWork *w,
    SOLNPSettings *stgs)
{
    solnp_int len = stgs->rescue ? w->n - 2 * w->nec : w->n;
    w->pb->n = len;
    SOLNPCost **ob = (SOLNPCost **)solnp_malloc(1 * sizeof(SOLNPCost *));
    solnp_float *p = (solnp_float *)solnp_malloc(1 * w->n * sizeof(solnp_float));
    memcpy(p, &w->p[w->nic], w->n * sizeof(solnp_float));

    solnp_float infeas_cand, temp;
    solnp_float *atemp = (solnp_float *)solnp_malloc(w->nc * len * sizeof(solnp_float));

    solnp_int i, j, tag;

    tag = 1;
    //// Initiation
    // for (j = 0; j < len; j++) {
    //     ob[j] = init_cost(w->nec, w->nic);
    //     memcpy(&p[j * w->n], &w->p[w->nic], w->n * sizeof(solnp_float));
    //     p[j * w->n + j] += stgs->delta;
    // }
    //// parallel calculation
    // calculate_scaled_cost(ob, p, w_sub->scale, stgs, w, len);
    ob[0] = init_cost(w->nec, w->nic);

    for (i = 0; i < len; i++)
    {

        p[i] += stgs->delta;

        calculate_scaled_cost(ob, p, w_sub->scale, stgs, w, 1);

        infeas_cand = calculate_infeas_scaledob(ob[0], w, p, w_sub->scale);
        if (infeas_cand <= stgs->tol_con && MIN(p[i] - w->pb->pl[i], w->pb->pu[i] - p[i]) > 0 && w_sub->ob_cand->obj > ob[0]->obj)
        {
            memcpy(w_sub->p_cand, w->p, w->nic * sizeof(solnp_float));
            memcpy(&w_sub->p_cand[w->nic], p, w->n * sizeof(solnp_float));
            copySOLNPCost(w_sub->ob_cand, ob[0]);
        }
        p[i] -= stgs->delta;
        // calculate Jacobian approximately
        w_sub->J->g[w->nic + i] = (ob[0]->obj - w->ob->obj) / stgs->delta;
        for (j = 0; j < w->nec; j++)
        {
            atemp[j + i * w_sub->nc] = (ob[0]->ec[j] - w->ob->ec[j]) / stgs->delta;
        }
        for (j = 0; j < w->nic; j++)
        {
            atemp[(w->nec + j) + (i)*w_sub->nc] = (ob[0]->ic[j] - w->ob->ic[j]) / stgs->delta;
        }

        if (w_sub->J->g[w->nic + i] != 0)
        {
            tag = 0;
        }
        // restore perturb

        temp = w->p[w->nic + i] * w_sub->scale[w->nc + 1 + i] - w_sub->J->g[w->nic + i] * w_sub->scale[0] / w_sub->scale[w->nc + 1 + i];
        if (w->nc > 0)
        {
            temp += w_sub->scale[0] * SOLNP(dot)(&atemp[i * w_sub->nc], w->l, w->nc) / w_sub->scale[w->nc + 1 + i];
        }
        if (temp > w->pb->pu[i] * w_sub->scale[w->nc + 1 + i])
        {
            temp = w->pb->pu[i] * w_sub->scale[w->nc + 1 + i];
        }
        else if (temp < w->pb->pl[i] * w_sub->scale[w->nc + 1 + i])
        {
            temp = w->pb->pl[i] * w_sub->scale[w->nc + 1 + i];
        }
        w->alm_crit += (temp - w->p[w->nic + i] * w_sub->scale[w->nc + 1 + i]) * (temp - w->p[w->nic + i] * w_sub->scale[w->nc + 1 + i]);
    }
    // copy atemp to J->a
    if (w->n > w->nec || stgs->rescue)
    {
        memcpy(&w_sub->J->a[w->nic * w->nc], atemp, w->nc * len * sizeof(solnp_float));
    }
    else
    {
        for (i = 0; i < w->n; i++)
        {
            for (j = 0; j < w->n - 1; j++)
            {
                w_sub->J->a[(w->nic + i) * w_sub->J->a_row + j] = atemp[i * w->nc + w_sub->constr_index[j]];
            }
        }
    }
    solnp_free(atemp);

    if (tag)
    {
        w->const_time++;
    }

    if (stgs->rescue == 0)
    {
        // Determin the condtion number of the Jocobi
        solnp_float *aT = (solnp_float *)solnp_malloc(w_sub->J->a_row * w_sub->npic * sizeof(solnp_float));
        SOLNP(transpose)
        (w_sub->J->a_row, w_sub->npic, w_sub->J->a, aT);
        solnp_float *aaT = (solnp_float *)solnp_malloc(w_sub->J->a_row * w_sub->J->a_row * sizeof(solnp_float));
        SOLNP(AB)
        (aaT, w_sub->J->a, aT, w_sub->J->a_row, w_sub->npic, w_sub->J->a_row);
        solnp_float *cond = (solnp_float *)solnp_malloc(sizeof(solnp_float));
        SOLNP(cond)
        (w_sub->J->a_row, aaT, cond);
        solnp_free(aT);
        solnp_free(aaT);
        // calculate condition number
        // TODO: calculate condition number
        if (*cond <= EPS)
        {
            solnp_printf("SOLNP+--> ");
            solnp_printf("Redundant constraints were found. Poor              \n");
            solnp_printf("         ");
            solnp_printf("intermediate results may result.  Suggest that you  \n");
            solnp_printf("         ");
            solnp_printf("remove redundant constraints and re-OPTIMIZE.       \n");
        }
        // free pointers
        solnp_free(cond);
    }
    solnp_free(p);
    free_cost(ob[0]);
    solnp_free(ob);
    w->pb->n = w->n;
    return 0;
}

solnp_int calculate_Jacob_zero_rescue(
    SUBNPWork *w_sub,
    SOLNPWork *w,
    SOLNPSettings *stgs)
{
    solnp_int i;
    calculate_Jacob_zero(w_sub, w, stgs);
    /*
    for (i = 0; i < w->n - 2*w->nec; i++) {
        w_sub->J->g[i] = 0;
    }*/
    for (i = 0; i < 2 * w->nec; i++)
    {
        w_sub->J->g[w->nic + w->n - 2 * w->nec + i] = w_sub->scale[1 + w->nc + w->n - 2 * w->nec + i]; // w_sub->scale[0];
    }
    for (i = 0; i < w->nec; i++)
    {
        w_sub->J->a[w_sub->J->a_row * (w->nic + w->n - 2 * w->nec) + i + i * w_sub->J->a_row] = -w_sub->scale[1 + w->nc + w->n - 2 * w->nec + i] / w_sub->scale[1 + i];
    }
    for (i = w->nec; i < 2 * w->nec; i++)
    {
        w_sub->J->a[w_sub->J->a_row * (w->nic + w->n - 2 * w->nec + i) + i - w->nec] = w_sub->scale[1 + w->nc + w->n - 2 * w->nec + i] / w_sub->scale[1 + i - w->nec];
    }
    return 0;
}

solnp_int calculate_Jacob_first(
    SUBNPWork *w_sub,
    SOLNPWork *w,
    SOLNPSettings *stgs)
{
    solnp_int i;
    solnp_float *grad = (solnp_float *)solnp_malloc((w->nc + 1) * w->n * sizeof(solnp_float));

    // calculate_scaled_grad(grad, w->p + w->nic, w_sub->scale, stgs, w);
    if (stgs->rs)
    {
        // use random sampling to calcualte gradient
        calculate_scaled_grad_random(grad, w->p + w->nic, w->ob, w_sub->scale, stgs, w);
    }
    else
    {
        // use coordinate-wise finite difference to calculate gradient
        calculate_scaled_grad(grad, w->p + w->nic, w_sub->scale, stgs, w);
    }
    if (w->exit == 1)
    {
        return 0;
    }
    for (i = 0; i < w->nic; i++)
    {
        w_sub->J->anew[i + (i + w->nec) * w_sub->npic] = -1.;
        w_sub->J->a[i + w->nec + i * w->nc] = -1;
        w_sub->J->g[i] = 0;
    }
    for (i = 0; i < w->nc; i++)
    {
        memcpy(w_sub->J->anew + w->nic + i * w_sub->npic, grad + (i + 1) * w->n, w->n * sizeof(solnp_float));
    }

    SOLNP(transpose)
    (w->n, w->nc, grad + w->n, w_sub->J->a + w->nc * w->nic);
    memcpy(w_sub->J->g + w->nic, grad, w->n * sizeof(solnp_float));

    calculate_alm_criterion(w, w_sub, grad);
    solnp_free(grad);
}

solnp_int calculate_Jacob_first_rescue(
    SUBNPWork *w_sub,
    SOLNPWork *w,
    SOLNPSettings *stgs)
{
    solnp_int i;
    if (stgs->rescue)
    {
        w->n -= 2 * w->nec;
    }
    calculate_Jacob_first(w_sub, w, stgs);
    if (stgs->rescue)
    {
        w->n += 2 * w->nec;
    }
    /*
    for (i = 0; i < w->n - 2*w->nec; i++) {
        w_sub->J->g[i] = 0;
    }*/
    for (i = 0; i < 2 * w->nec; i++)
    {
        w_sub->J->g[w->nic + w->n - 2 * w->nec + i] = w_sub->scale[1 + w->nc + w->n - 2 * w->nec + i]; // w_sub->scale[0];
    }
    for (i = 0; i < w->nec; i++)
    {
        w_sub->J->a[w_sub->J->a_row * (w->nic + w->n - 2 * w->nec) + i + i * w_sub->J->a_row] = -w_sub->scale[1 + w->nc + w->n - 2 * w->nec + i] / w_sub->scale[1 + i];
    }
    for (i = w->nec; i < 2 * w->nec; i++)
    {
        w_sub->J->a[w_sub->J->a_row * (w->nic + w->n - 2 * w->nec + i) + i - w->nec] = w_sub->scale[1 + w->nc + w->n - 2 * w->nec + i] / w_sub->scale[1 + i - w->nec];
    }
    return 0;
}

solnp_float calculate_ALM(
    SOLNPCost *ob,
    SOLNPSettings *stgs,
    solnp_float *p,
    const SOLNPWork *w,
    const SUBNPWork *w_sub)
{

    solnp_float result = ob->obj;
    if (w->nc)
    {
        solnp_float *ap = (solnp_float *)solnp_malloc(w_sub->J->a_row * sizeof(solnp_float));
        SOLNP(Ax)
        (ap, w_sub->J->a, p, w_sub->J->a_row, w_sub->npic);
        if (w->nic)
        {
            solnp_float *ic = (solnp_float *)solnp_malloc(w->nic * sizeof(solnp_float));
            memcpy(ic, ob->ic, w->nic * sizeof(solnp_float));
            SOLNP(add_scaled_array)
            (ic, p, w->nic, -1.0);
            if (w->n > w->nec || stgs->rescue)
            {
                SOLNP(add_scaled_array)
                (ic, ap + w->nec, w->nic, -1.0);
                SOLNP(add_scaled_array)
                (ic, w_sub->b + w->nec, w->nic, 1.0);
            }
            result = result + w->rho * SOLNP(norm_sq)(ic, w->nic);
            result = result - SOLNP(dot)(w->l + w->nec, ic, w->nic);
            solnp_free(ic);
        }
        if (w->nec)
        {
            solnp_float *ec = (solnp_float *)solnp_malloc(w->nec * sizeof(solnp_float));
            memcpy(ec, ob->ec, w->nec * sizeof(solnp_float));
            if (w->n > w->nec || stgs->rescue)
            {
                SOLNP(add_scaled_array)
                (ec, ap, w->nec, -1.0);
                SOLNP(add_scaled_array)
                (ec, w_sub->b, w->nec, 1.0);
            }
            result = result + w->rho * SOLNP(norm_sq)(ec, w->nec);
            result = result - SOLNP(dot)(w->l, ec, w->nec);
            solnp_free(ec);
        }
        solnp_free(ap);
    }
    return result;
}
/*
solnp_float calculate_ALM_rescue
(
    SOLNPCost* ob,
    SOLNPSettings* stgs,
    solnp_float* p,
    const SOLNPWork* w,
    const SUBNPWork* w_sub
) {

    solnp_int i;
    solnp_float result = 0;
    for (i = 0; i < w->nec; i++) {
        result += (w_sub->scale[1 + w->nc + w->n - w->nec + i] * w->p[w->n - w->nec + i]) * (w_sub->scale[1 + w->nc + w->n - w->nec + i] * w->p[w->n - w->nec + i]);
    }
    result /= w_sub->scale[0];
    if (w->nc) {
        solnp_float* ap = (solnp_float*)solnp_malloc(w_sub->J->a_row * sizeof(solnp_float));
        SOLNP(Ax)(ap, w_sub->J->a, p, w_sub->J->a_row, w_sub->npic);
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
}*/
solnp_int calculate_ALMgradient_zero(
    SUBNPWork *w_sub,
    SOLNPWork *w,
    SOLNPSettings *stgs,
    solnp_float *g,
    solnp_float j)
{   
    // Calculate the gradient of Augmented Lagrangian Function

    SOLNPCost **obm = (SOLNPCost **)solnp_malloc(1 * sizeof(SOLNPCost *));
    SOLNPCost** obm_backward = (SOLNPCost**)solnp_malloc(1 * sizeof(SOLNPCost*));
    solnp_int i, tag;
    tag = 1;
    solnp_int len = stgs->rescue ? w->n - 2 * w->nec : w->n;
    solnp_float alm,alm_backward;
    solnp_float infeas, temp;
    solnp_float *contemp = (solnp_float *)solnp_malloc(w->nc * sizeof(solnp_float));
    solnp_float *p = (solnp_float *)solnp_malloc(1 * w->n * sizeof(solnp_float));
    memcpy(p, &w->p[w->nic], w->n * sizeof(solnp_float));

    obm[0] = init_cost(w->nec, w->nic);
    obm_backward[0] = init_cost(w->nec, w->nic);


    for (i = 0; i < len; i++)
    {

        p[i] += stgs->delta;
        calculate_scaled_cost(obm, p, w_sub->scale, stgs, w, 1);

        // contemp and temp variable is used for calculate the projected alm gradient(as stop criterion)
       
        if (w->nec)
        {
            memcpy(contemp, obm[0]->ec, w->nec * sizeof(solnp_float));
        }
        if (w->nic)
        {
            memcpy(&contemp[w->nec], obm[0]->ic, w->nic * sizeof(solnp_float));
        }

        if (stgs->rescue == 0 && obm[0]->obj != w->ob->obj)
        {
            tag = 0;
        }

        // calculate the ALM function value
        alm = obm[0]->obj;
        solnp_float *ptemp = (solnp_float *)solnp_malloc(w_sub->J->npic * sizeof(solnp_float));
        memcpy(ptemp, w->p, w->nic * sizeof(solnp_float));
        memcpy(&ptemp[w->nic], p, w->n * sizeof(solnp_float));
        if (w_sub->nc > 0)
        {
            alm = calculate_ALM(obm[0], stgs, ptemp, w, w_sub);
        }
        solnp_free(ptemp);

        // calculate the infeasibility
        infeas = calculate_infeas_scaledob(obm[0], w, p, w_sub->scale);

        {
            // Record the best point
            if (infeas <= stgs->tol_con && MIN(p[i] - w->pb->pl[i], w->pb->pu[i] - p[i]) > 0 && w_sub->ob_cand->obj > obm[0]->obj)
            {
                memcpy(w_sub->p_cand, w->p, w->nic * sizeof(solnp_float));
                memcpy(&w_sub->p_cand[w->nic], p, w->n * sizeof(solnp_float));
                copySOLNPCost(w_sub->ob_cand, obm[0]);
            }

            // calculate gradient approximately
            if (stgs->cen_diff) {
                //use centre difference to calculate the gradient
                p[i] -= 2 * stgs->delta;
                calculate_scaled_cost(obm_backward, p, w_sub->scale, stgs, w, 1);
                alm_backward = obm_backward[0]->obj;
                solnp_float* ptemp = (solnp_float*)solnp_malloc(w_sub->J->npic * sizeof(solnp_float));
                memcpy(ptemp, w->p, w->nic * sizeof(solnp_float));
                memcpy(&ptemp[w->nic], p, w->n * sizeof(solnp_float));
                if (w_sub->nc > 0)
                {
                    alm_backward = calculate_ALM(obm[0], stgs, ptemp, w, w_sub);
                }
                solnp_free(ptemp);

                g[w->nic + i] = (alm - alm_backward) / (2*stgs->delta);

                p[i] += 2 * stgs->delta;
            }
            else {
                g[w->nic + i] = (alm - j) / stgs->delta;
            }
        }
        /*
        if (w->exit == 1) {
            solnp_free(contemp);
            break;
        }*/
        p[i] -= stgs->delta;

        temp = w->p[w->nic + i] * w_sub->scale[w->nc + 1 + i] - ((obm[0]->obj - w->ob->obj) / stgs->delta) * w_sub->scale[0] / w_sub->scale[w->nc + 1 + i];
        if (w->nc > 0)
        {
            if (w->nec)
            {
                SOLNP(add_scaled_array)
                (contemp, w->ob->ec, w->nec, -1.0);
            }
            if (w->nic)
            {
                SOLNP(add_scaled_array)
                (&contemp[w->nec], w->ob->ic, w->nic, -1.0);
            }
            SOLNP(scale_array)
            (contemp, 1 / stgs->delta, w->nc);
            temp += w_sub->scale[0] * SOLNP(dot)(contemp, w->l, w->nc) / w_sub->scale[w->nc + 1 + i];
        }
        if (temp > w->pb->pu[i] * w_sub->scale[w->nc + 1 + i])
        {
            temp = w->pb->pu[i] * w_sub->scale[w->nc + 1 + i];
        }
        else if (temp < w->pb->pl[i] * w_sub->scale[w->nc + 1 + i])
        {
            temp = w->pb->pl[i] * w_sub->scale[w->nc + 1 + i];
        }

        {
            w->alm_crit += (temp - w->p[w->nic + i] * w_sub->scale[w->nc + 1 + i]) * (temp - w->p[w->nic + i] * w_sub->scale[w->nc + 1 + i]);
        }
    }
    // if (w->nic > 0.5) {
    //     for (k = 0; k < w->nic; k++) {
    //         g[k] = 0;
    //     }
    // }

    if (tag)
    {
        w->const_time++;
    }

    solnp_free(p);
    solnp_free(contemp);
    free_cost(obm[0]);
    free_cost(obm_backward[0]);
    solnp_free(obm);
    solnp_free(obm_backward);
    return 0;
}

solnp_int calculate_ALMgradient_zero_rescue(
    SUBNPWork *w_sub,
    SOLNPWork *w,
    SOLNPSettings *stgs,
    solnp_float *g,
    solnp_float j)
{
    calculate_ALMgradient_zero(w_sub, w, stgs, g, j);
    for (solnp_int i = w->n - 2 * w->nec; i < w->n - 1 * w->nec; i++)
    {
        w_sub->J->g[w->nic + i] = w->pen_l1 * w_sub->scale[1 + w->nc + i] + w->l[i - (w->n - 2 * w->nec)] * w_sub->scale[1 + w->nc + i] / w_sub->scale[1 + i - (w->n - 2 * w->nec)] - 2 * w->rho * w->ob->ec[i - (w->n - 2 * w->nec)] * w_sub->scale[1 + w->nc + i] / w_sub->scale[1 + i - (w->n - 2 * w->nec)];
    }
    for (solnp_int i = w->n - 1 * w->nec; i < w->n; i++)
    {
        w_sub->J->g[w->nic + i] = w->pen_l1 * w_sub->scale[1 + w->nc + i] - w->l[i - (w->n - w->nec)] * w_sub->scale[1 + w->nc + i] / w_sub->scale[1 + i - (w->n - w->nec)] + 2 * w->rho * w->ob->ec[i - (w->n - w->nec)] * w_sub->scale[1 + w->nc + i] / w_sub->scale[1 + i - (w->n - w->nec)];
    }
    return 0;
}

solnp_int calculate_ALMgradient_first(
    SUBNPWork *w_sub,
    SOLNPWork *w,
    SOLNPSettings *stgs,
    solnp_float *g)
{
    solnp_int i;
    solnp_float *grad = (solnp_float *)solnp_malloc((w->nc + 1) * w->n * sizeof(solnp_float));
    solnp_float *temp_r = (solnp_float *)solnp_malloc(w->nc * sizeof(solnp_float));
    solnp_float *temp_l = (solnp_float *)solnp_malloc(w->nc * w->n * sizeof(solnp_float));

    if (stgs->rs)
    {
        calculate_scaled_grad_random(grad, w->p + w->nic, w->ob, w_sub->scale, stgs, w);
    }
    else
    {
        calculate_scaled_grad(grad, w->p + w->nic, w_sub->scale, stgs, w);
    }
    //
    if (w->exit == 1)
    {
        return 0;
    }

    for (i = 0; i < w->nic; i++)
    {
        g[i] = 0;
        w_sub->J->anew[i + (i + w->nec) * w_sub->npic] = -1.;
    }
    for (i = 0; i < w->nc; i++)
    {
        memcpy(w_sub->J->anew + w->nic + i * w_sub->npic, grad + (i + 1) * w->n, w->n * sizeof(solnp_float));
    }

    SOLNP(transpose)
    (w->n, w->nc, grad + w->n, temp_l);
    memcpy(g + w->nic, grad, w->n * sizeof(solnp_float));

    // calculate ALM stop criterion
    calculate_alm_criterion(w, w_sub, grad);
    solnp_free(grad);

    if (w->nc > 0)
    {
        SOLNP(add_scaled_array)
        (temp_l, w_sub->J->a + w->nc * w->nic, w->nc * w->n, -1.);
        SOLNP(set_as_scaled_array)
        (temp_r, w->l, -1., w->nc);
        SOLNP(add_scaled_array)
        (temp_r, w->ob->ec, w->nec, w->rho);
        SOLNP(add_scaled_array)
        (temp_r + w->nec, w->ob->ic, w->nic, w->rho);
        SOLNP(add_scaled_array)
        (temp_r + w->nec, w->p, w->nic, -w->rho);
        solnp_float *temp = (solnp_float *)solnp_malloc(w_sub->npic * sizeof(solnp_float));
        SOLNP(AB)
        (temp, temp_r, temp_l, 1, w->nc, w->n);
        SOLNP(add_scaled_array)
        (g + w->nic, temp, w_sub->npic, 1.);
        solnp_free(temp);
    }

    solnp_free(temp_r);
    solnp_free(temp_l);
}

solnp_int calculate_ALMgradient_first_rescue(
    SUBNPWork *w_sub,
    SOLNPWork *w,
    SOLNPSettings *stgs,
    solnp_float *g,
    solnp_float j)
{
    if (stgs->rescue)
    {
        w->n -= 2 * w->nec;
    }
    calculate_ALMgradient_first(w_sub, w, stgs, g);
    if (stgs->rescue)
    {
        w->n += 2 * w->nec;
    }
    for (solnp_int i = w->n - 2 * w->nec; i < w->n - 1 * w->nec; i++)
    {
        w_sub->J->g[w->nic + i] = w->pen_l1 * w_sub->scale[1 + w->nc + i] + w->l[i - (w->n - 2 * w->nec)] * w_sub->scale[1 + w->nc + i] / w_sub->scale[1 + i - (w->n - 2 * w->nec)] - 2 * w->rho * w->ob->ec[i - (w->n - 2 * w->nec)] * w_sub->scale[1 + w->nc + i] / w_sub->scale[1 + i - (w->n - 2 * w->nec)];
    }
    for (solnp_int i = w->n - 1 * w->nec; i < w->n; i++)
    {
        w_sub->J->g[w->nic + i] = w->pen_l1 * w_sub->scale[1 + w->nc + i] - w->l[i - (w->n - w->nec)] * w_sub->scale[1 + w->nc + i] / w_sub->scale[1 + i - (w->n - w->nec)] + 2 * w->rho * w->ob->ec[i - (w->n - w->nec)] * w_sub->scale[1 + w->nc + i] / w_sub->scale[1 + i - (w->n - w->nec)];
    }
    return 0;
}

solnp_int calculate_ALM_hess(
    SUBNPWork *w_sub,
    SOLNPWork *w,
    SOLNPSettings *stgs,
    solnp_float *h)
{
    solnp_int i, j, k;
    solnp_float *hess = (solnp_float *)solnp_malloc(w->n * w->n * (w->nc + 1) * sizeof(solnp_float));
    calculate_scaled_hess(hess, w->p + w->nic, w_sub->scale, stgs, w);

    for (k = 1; k < w->nc + 1; k++)
    {
        solnp_float con_value;
        if (k - 1 < w->nec)
        {
            con_value = w->ob->ec[k - 1];
        }
        else
        {
            con_value = w->ob->ic[k - 1 - w->nec] - w->p[k - 1 - w->nec];
        }
        SOLNP(add_scaled_array)
        (hess, hess + k * w->n * w->n, w->n * w->n, w->rho * con_value - w->l[k - 1]);
    }

    for (i = 0; i < w_sub->npic; i++)
    {
        for (j = 0; j < w_sub->npic; j++)
        {
            if (i < w->nic || j < w->nic)
            {
                h[i + j * w_sub->npic] = 0;
            }
            else
            {
                h[i + j * w_sub->npic] = hess[(i - w->nic) + (j - w->nic) * w->n];
            }
        }
    }
    for (k = 0; k < w->nc; k++)
    {
        SOLNP(rank1update)
        (w_sub->npic, h, w->rho, w_sub->J->anew + k * w_sub->npic);
    }
    solnp_free(hess);
}

void BFGSudpate(
    SOLNPWork *w,
    SUBNPWork *w_sub,
    SOLNPSettings *stgs,
    solnp_float *g,
    solnp_float *yg,
    solnp_float *sx)
{
    SOLNP(add_scaled_array)
    (yg, g, w_sub->J->npic, -1.0);
    SOLNP(scale_array)
    (yg, -1.0, w_sub->J->npic);
    SOLNP(add_scaled_array)
    (sx, w->p, w_sub->J->npic, -1.0);
    SOLNP(scale_array)
    (sx, -1.0, w_sub->J->npic);
    solnp_float sc[2];
    solnp_float *temp = (solnp_float *)solnp_malloc(w_sub->J->npic * sizeof(solnp_float));
    SOLNP(Ax)
    (temp, w->h, sx, w_sub->J->npic, w_sub->J->npic);
    sc[0] = SOLNP(dot)(sx, temp, w_sub->J->npic);
    sc[1] = SOLNP(dot)(sx, yg, w_sub->J->npic);
    if (sc[0] * sc[1] > 0)
    {
        memcpy(sx, temp, w_sub->J->npic * sizeof(solnp_float));
        // two rank 1 updates
        SOLNP(rank1update)
        (w_sub->J->npic, w->h, -1 / sc[0], sx);
        SOLNP(rank1update)
        (w_sub->J->npic, w->h, 1 / sc[1], yg);
    }
    solnp_free(temp);
}

// void LBFGS_H_inv_v
//(
//     solnp_int m,
//     solnp_int n,
//     solnp_int start_index,
//     solnp_float* V,
//     solnp_float* Y,
//     solnp_float* S
//) {
//     // This subrontine calculate the H^{-1} vector product given the past data
// }

solnp_float fun_along_2d(
    SOLNPWork *w,
    SOLNPSettings *stgs,
    SUBNPWork *w_sub,
    solnp_float *d1,
    solnp_float *d2,
    solnp_float *x,
    solnp_float *coeff)
{
    // This subroutine calculate the function value along two direction
    solnp_float *p = (solnp_float *)solnp_malloc(w->n * sizeof(solnp_float));
    solnp_float val;

    // calculate the target point
    memcpy(p, x, w->n * sizeof(solnp_float));
    SOLNP(add_scaled_array)
    (p, d1, w->n, coeff[0]);
    SOLNP(add_scaled_array)
    (p, d2, w->n, coeff[1]);

    // calculate the value at p
    SOLNPCost *Ob = init_cost(w->nec, w->nic);
    calculate_scaled_cost(&Ob, p, w_sub->scale, stgs, w, 1);
    val = Ob->obj;

    free_cost(Ob);
    solnp_free(p);
    return val;
}

solnp_float *interpolate2d(
    SOLNPWork *w,
    SOLNPSettings *stgs,
    SUBNPWork *w_sub,
    solnp_float radius,
    solnp_float *d1,
    solnp_float *d2,
    solnp_float *x,
    solnp_float val)
{
    // This subroutine do quadratic interpolation to the function.
    // Input:   two directions d1 and d2
    //          current point x
    //          current value val
    //          interpolation radius radius
    // Output:  linear coeffience g and quadratic coeffience Q, stored in the vector g_and_Q
    // REMEMBER TO FREE g_and_Q after using it!
    solnp_float *g_and_Q = (solnp_float *)solnp_malloc(6 * sizeof(solnp_float));
    solnp_float *init_pt = (solnp_float *)solnp_malloc(10 * sizeof(solnp_float));
    solnp_int i, j, k;
    // Set up init_pt

    for (i = 0; i < 5; i++)
    {
        init_pt[2 * i] = cos(2 * i * PI / 5) * radius;
        init_pt[2 * i + 1] = sin(2 * i * PI / 5) * radius;
    }

    solnp_float *coeff_matrix = (solnp_float *)solnp_malloc(25 * sizeof(solnp_float));
    solnp_float *fval = (solnp_float *)solnp_malloc(5 * sizeof(solnp_float));

    // Set up the coeff matrix and favl
    for (i = 0; i < 5; i++)
    {
        coeff_matrix[i] = init_pt[2 * i];
        coeff_matrix[i + 5] = init_pt[2 * i + 1];
        fval[i] = fun_along_2d(w, stgs, w_sub, d1, d2, x, &init_pt[2 * i]) - val;
        for (j = 0; j < 2; j++)
        {
            for (k = 0; k <= j; k++)
            {
                if (k == j)
                {
                    coeff_matrix[i + 5 * (j * (j + 1) / 2 + k + 2)] = init_pt[2 * i + j] * init_pt[2 * i + k];
                }
                else
                {
                    coeff_matrix[i + 5 * (j * (j + 1) / 2 + k + 2)] = 2 * init_pt[2 * i + j] * init_pt[2 * i + k];
                }
            }
        }
    }

    // Solve linear equations coeff*x = fval
    SOLNP(solve_general_lin_sys)
    (5, coeff_matrix, fval);
    // Record the result
    memcpy(g_and_Q, fval, 2 * sizeof(solnp_float));
    g_and_Q[2] = fval[2];
    g_and_Q[3] = fval[3];
    g_and_Q[4] = fval[3];
    g_and_Q[5] = fval[4];

    // Free variable and return
    solnp_free(fval);
    solnp_free(coeff_matrix);
    solnp_free(init_pt);
    return g_and_Q;
}

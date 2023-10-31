#include "solnp.h"
#include "matrix.h"
#include "mex.h"
#include "linalg.h"
#include "solnp_util.h"

solnp_float cost_time = 0; // global cost timer

// #if !(DLONG > 0)
// #ifndef DLONG
// // this memory must be freed
// static solnp_int *cast_to_solnp_int_arr(mwIndex *arr, solnp_int len) {
//     solnp_int i;
//     solnp_int *arr_out = (solnp_int *)solnp_malloc(sizeof(solnp_int) * len);
//     for (i = 0; i < len; i++) {
//         arr_out[i] = (solnp_int)arr[i];
//     }
//     return arr_out;
// }
// #endif

// #if SFLOAT > 0
// /* this memory must be freed */
// static solnp_float *cast_to_solnp_float_arr(double *arr, solnp_int len) {
//     solnp_int i;
//     solnp_float *arr_out = (solnp_float *)solnp_malloc(sizeof(solnp_float) * len);
//     for (i = 0; i < len; i++) {
//         arr_out[i] = (solnp_float)arr[i];
//     }
//     return arr_out;
// }

// static double *cast_to_double_arr(solnp_float *arr, solnp_int len) {
//     solnp_int i;
//     double *arr_out = (double *)solnp_malloc(sizeof(double) * len);
//     for (i = 0; i < len; i++) {
//         arr_out[i] = (double)arr[i];
//     }
//     return arr_out;
// }
// #endif

static void set_output_field(mxArray **pout, solnp_float *out, solnp_int m, solnp_int n)
{

    *pout = mxCreateDoubleMatrix(m, n, mxREAL);
    for (int i = 0; i < m * n; i++)
    {
        mxGetPr(*pout)[i] = (double)out[i];
    }
}

void mexCallMatlabCost(SOLNPCost **c, solnp_float *p, solnp_int np, solnp_int nfeval, solnp_int action, mxArray *fun)
{
    // initiate cost function
    // action == 0 : do nothing
    // action == 1 : initiate cost
    // action == -1 : close cost handle
    static mxArray *cost_fun = NULL;
    if (action == 1)
    {
        cost_fun = fun;
        return;
    }
    if (action == -1)
    {
        cost_fun = NULL;
        return;
    }

    SOLNP(timer)
    cost_timer;
    SOLNP(tic)
    (&cost_timer);
    solnp_int i, j;
    // mxArray *lhs, *rhs[3];
    mxArray *lhs, *rhs[2];

    nfeval = 1;

    rhs[0] = cost_fun;
    rhs[1] = mxCreateDoubleMatrix(np, nfeval, mxREAL);
    // rhs[2] = mxCreateDoubleMatrix(1, 1, mxREAL);

    for (j = 0; j < nfeval; j++)
    {
        for (i = 0; i < np; i++)
        {
            mxGetPr(rhs[1])[i + j * np] = (double)p[i + j * np];
        }
    }
    // *mxGetPr(rhs[2]) = (double)nfeval;
    lhs = SOLNP_NULL;
    // mexCallMATLAB(1, &lhs, 3, rhs, "feval");
    mexCallMATLAB(1, &lhs, 2, rhs, "feval");

    for (j = 0; j < nfeval; j++)
    {
        // double *result = mxGetPr(lhs);

        c[j]->obj = (solnp_float)mxGetPr(lhs)[j * (1 + c[j]->nic + c[j]->nec)];

        for (i = 0; i < c[j]->nec; i++)
        {
            c[j]->ec[i] = (solnp_float)mxGetPr(lhs)[i + 1 + j * (1 + c[j]->nic + c[j]->nec)];
        }
        for (i = 0; i < c[j]->nic; i++)
        {
            c[j]->ic[i] = (solnp_float)mxGetPr(lhs)[i + 1 + c[j]->nec + j * (1 + c[j]->nic + c[j]->nec)];
        }
    }
    mxDestroyArray(rhs[1]);
    mxDestroyArray(rhs[2]);

    if (lhs != SOLNP_NULL)
    {
        mxDestroyArray(lhs);
    }
    cost_time += SOLNP(tocq)(&cost_timer) / 1e3;
}

void mexCallMatlabGradient(solnp_float *g, solnp_float *p, solnp_int np, solnp_int ngeval, solnp_int action, mxArray *fun)
{
    // initiate cost function
    // action == 0 : do nothing
    // action == 1 : initiate cost
    // action == -1 : close cost handle
    static mxArray *grad = NULL;
    if (action == 1)
    {
        grad = fun;
        return;
    }
    if (action == -1)
    {
        grad = NULL;
        return;
    }

    SOLNP(timer)
    cost_timer;
    SOLNP(tic)
    (&cost_timer);
    solnp_int i;
    mxArray *lhs, *rhs[2];

    rhs[0] = grad;
    rhs[1] = mxCreateDoubleMatrix(np, 1, mxREAL);
    for (i = 0; i < np; i++)
    {
        mxGetPr(rhs[1])[i] = (double)p[i];
    }

    lhs = SOLNP_NULL;
    mexCallMATLAB(1, &lhs, 2, rhs, "feval");
    for (i = 0; i < np * ngeval; i++)
    {
        // double *result = mxGetPr(lhs);
        g[i] = (solnp_float)mxGetPr(lhs)[i];
    }
    mxDestroyArray(rhs[1]);

    if (lhs != SOLNP_NULL)
    {
        mxDestroyArray(lhs);
    }
    cost_time += SOLNP(tocq)(&cost_timer) / 1e3;
}

void mexCallMatlabHessian(solnp_float *h, solnp_float *p, solnp_int np, solnp_int nheval, solnp_int action, mxArray *fun)
{
    // initiate cost function
    // action == 0 : do nothing
    // action == 1 : initiate cost
    // action == -1 : close cost handle

    static mxArray *hess = NULL;
    if (action == 1)
    {
        hess = fun;
        return;
    }
    if (action == -1)
    {
        hess = NULL;
        return;
    }

    SOLNP(timer)
    cost_timer;
    SOLNP(tic)
    (&cost_timer);
    solnp_int i;
    mxArray *lhs, *rhs[2];
    rhs[0] = hess;

    rhs[1] = mxCreateDoubleMatrix(np, 1, mxREAL);
    for (i = 0; i < np; i++)
    {
        mxGetPr(rhs[1])[i] = (double)p[i];
    }

    lhs = SOLNP_NULL;
    mexCallMATLAB(1, &lhs, 2, rhs, "feval");
    for (i = 0; i < np * np * nheval; i++)
    {
        // double *result = mxGetPr(lhs);
        h[i] = (solnp_float)mxGetPr(lhs)[i];
    }
    mxDestroyArray(rhs[1]);

    if (lhs != SOLNP_NULL)
    {
        mxDestroyArray(lhs);
    }
    cost_time += SOLNP(tocq)(&cost_timer) / 1e3;
}

// mex file for matlab function: [p,jh,l,h,ic]=solnp(cnstr,op,l,h)
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    SOLNP(timer)
    total_timer;
    SOLNP(tic)
    (&total_timer);

    if (nrhs < 1)
    {
        mexErrMsgTxt("Syntax error");
    }

    const mwSize one[1] = {1};
    const int num_sol_fields = 17;
    const char *sol_fields[] = {"p", "best_fea_p", "jh", "ch", "l", "h", "ic", "iter", "count_cost", "count_grad", "count_hess", "constraint", "obj", "status", "solve_time", "restart_time", "cost_his"};

    const mxArray *cnstr;
    const mxArray *op;
    const mxArray *fun;

    const mxArray *pblmx;
    solnp_float *pbl;

    const mxArray *pbumx;
    solnp_float *pbu;

    const mxArray *p0mx;
    solnp_float *p;

    const mxArray *iblmx;
    solnp_float *ibl;

    const mxArray *ibumx;
    solnp_float *ibu;

    const mxArray *ib0mx;
    solnp_float *ib0;

    solnp_int ilc = 0;
    solnp_int iuc = 0;
    solnp_int ic;

    mxArray *tmp;

    const size_t *dims;
    solnp_int i;

    const mxArray *lmex = SOLNP_NULL;
    const mxArray *hmex = SOLNP_NULL;
    solnp_float *l;
    solnp_float *h;

    solnp_int ns;
    solnp_int np;
    solnp_int nic = 0;
    solnp_int nec = 0;
    solnp_int nc = 0;

    solnp_int *Ipc = (solnp_int *)solnp_malloc(2 * sizeof(solnp_int));
    memset(Ipc, 0, 2 * sizeof(solnp_int));
    solnp_int *Ipb = (solnp_int *)solnp_malloc(2 * sizeof(solnp_int));
    memset(Ipb, 0, 2 * sizeof(solnp_int));

    char *om = (char *)solnp_malloc(sizeof(char) * 512); // error message
    sprintf(om, "SOLNP+--> ");

    cnstr = prhs[0];

    pblmx = (mxArray *)mxGetField(cnstr, 0, "pbl");
    if (pblmx && !mxIsEmpty(pblmx))
    {
        Ipc[0] = 1;
        Ipb[0] = 1;
        ns = (solnp_int)mxGetNumberOfDimensions(pblmx);
        dims = mxGetDimensions(pblmx);
        np = (solnp_int)dims[0];

        if (ns > 1 && dims[0] == 1)
        {
            np = (solnp_int)dims[1];
        }

        pbl = (solnp_float *)solnp_malloc(sizeof(solnp_float) * np);
        for (i = 0; i < np; i++)
        {
            pbl[i] = (solnp_float)mxGetPr(pblmx)[i];
        }
    }

    pbumx = (mxArray *)mxGetField(cnstr, 0, "pbu");
    if (pbumx && !mxIsEmpty(pbumx))
    {
        Ipc[1] = 1;
        Ipb[0] = 1;
        ns = (solnp_int)mxGetNumberOfDimensions(pbumx);
        dims = mxGetDimensions(pbumx);
        np = (solnp_int)dims[0];

        if (ns > 1 && dims[0] == 1)
        {
            np = (solnp_int)dims[1];
        }

        pbu = (solnp_float *)solnp_malloc(sizeof(solnp_float) * np);
        for (i = 0; i < np; i++)
        {
            pbu[i] = (solnp_float)mxGetPr(pbumx)[i];
        }
    }

    p0mx = (mxArray *)mxGetField(cnstr, 0, "p0");
    if (p0mx && !mxIsEmpty(p0mx))
    {

        ns = (solnp_int)mxGetNumberOfDimensions(p0mx);
        dims = mxGetDimensions(p0mx);
        np = (solnp_int)dims[0];

        if (ns > 1 && dims[0] == 1)
        {
            np = (solnp_int)dims[1];
        }

        p = (solnp_float *)solnp_malloc(sizeof(solnp_float) * np);
        for (i = 0; i < np; i++)
        {
            p[i] = (solnp_float)mxGetPr(p0mx)[i];
        }
    }

    if (Ipc[0] && Ipc[1])
    {

        if (!(p0mx && !mxIsEmpty(p0mx)))
        {
            p = (solnp_float *)solnp_malloc(sizeof(solnp_float) * np);
            memcpy(p, pbl, sizeof(solnp_float) * np);
            SOLNP(add_scaled_array)
            (p, pbu, np, 1);
            SOLNP(scale_array)
            (p, 0.5, np);
        }

        for (i = 0; i < np; i++)
        {
            if (mxIsInf(p[i]))
            {
                sprintf(om + strlen(om), "The user does not provide initiate point! \n");
                mexErrMsgTxt(om);
            }
        }
    }
    else
    {
        if (!(p0mx && !mxIsEmpty(p0mx)))
        {
            sprintf(om + strlen(om), "The user does not provide initiate point! \n");
            mexErrMsgTxt(om);
        }
        if (Ipc[0] == 0)
        {

            pbl = (solnp_float *)solnp_malloc(sizeof(solnp_float) * np);
            for (i = 0; i < np; i++)
            {
                pbl[i] = -INFINITY;
            }
        }
        if (Ipc[1] == 0)
        {
            pbu = (solnp_float *)solnp_malloc(sizeof(solnp_float) * np);
            for (i = 0; i < np; i++)
            {
                pbu[i] = INFINITY;
            }
        }
    }

    if (Ipb[0] > 0.5)
    {

        for (i = 0; i < np; i++)
        {
            if (pbu[i] <= pbl[i])
            {
                sprintf(om + strlen(om), "The lower bounds of the parameter constraints \n");
                sprintf(om + strlen(om), "          must be strictly less than the upper bounds. \n");
                mexErrMsgTxt(om);
            }
            else if (p[i] <= pbl[i] || p[i] >= pbu[i])
            {
                sprintf(om + strlen(om), "Initial parameter values must be within the bounds \n");
                mexErrMsgTxt(om);
            }
        }
    }

    iblmx = (mxArray *)mxGetField(cnstr, 0, "ibl");
    if (iblmx && !mxIsEmpty(iblmx))
    {
        ilc = 1;
        ns = (solnp_int)mxGetNumberOfDimensions(iblmx);
        dims = mxGetDimensions(iblmx);
        nic = (solnp_int)dims[0];

        if (ns > 1 && dims[0] == 1)
        {
            nic = (solnp_int)dims[1];
        }

        ibl = (solnp_float *)solnp_malloc(sizeof(solnp_float) * nic);
        for (i = 0; i < nic; i++)
        {
            ibl[i] = (solnp_float)mxGetPr(iblmx)[i];
        }
    }

    ibumx = (mxArray *)mxGetField(cnstr, 0, "ibu");
    if (ibumx && !mxIsEmpty(ibumx))
    {
        iuc = 1;
        ns = (solnp_int)mxGetNumberOfDimensions(ibumx);
        dims = mxGetDimensions(ibumx);
        nic = (solnp_int)dims[0];

        if (ns > 1 && dims[0] == 1)
        {
            nic = (solnp_int)dims[1];
        }

        ibu = (solnp_float *)solnp_malloc(sizeof(solnp_float) * nic);
        for (i = 0; i < nic; i++)
        {
            ibu[i] = (solnp_float)mxGetPr(ibumx)[i];
        }
    }

    ib0mx = (mxArray *)mxGetField(cnstr, 0, "ib0");
    if (ib0mx && !mxIsEmpty(ib0mx))
    {
        ns = (solnp_int)mxGetNumberOfDimensions(ib0mx);
        dims = mxGetDimensions(ib0mx);
        nic = (solnp_int)dims[0];

        if (ns > 1 && dims[0] == 1)
        {
            nic = (solnp_int)dims[1];
        }

        ib0 = (solnp_float *)solnp_malloc(sizeof(solnp_float) * nic);
        for (i = 0; i < nic; i++)
        {
            ib0[i] = (solnp_float)mxGetPr(ib0mx)[i];
        }
    }

    ic = ilc || iuc;

    if (ilc && iuc)
    {

        if (!(ib0mx && !mxIsEmpty(ib0mx)))
        {

            ib0 = (solnp_float *)solnp_malloc(sizeof(solnp_float) * nic);
            memcpy(ib0, ibl, sizeof(solnp_float) * nic);
            SOLNP(add_scaled_array)
            (ib0, ibu, nic, 1);
            SOLNP(scale_array)
            (ib0, 0.5, nic);

            for (i = 0; i < nic; i++)
            {
                if (mxIsInf(ib0[i]))
                {
                    sprintf(om + strlen(om), "The user does not provided initiate value of inequality constrains! \n");
                    mexErrMsgTxt(om);
                }
            }
        }

        for (i = 0; i < nic; i++)
        {
            if (ibu[i] <= ibl[i])
            {
                sprintf(om + strlen(om), "The lower bounds of the inequality constraints \n");
                sprintf(om + strlen(om), "          must be strictly less than the upper bounds. \n");
                mexErrMsgTxt(om);
            }
        }
    }

    else if (ic == 1)
    {

        if (!(ib0mx && !mxIsEmpty(ib0mx)))
        {
            sprintf(om + strlen(om), "The user does not provided initiate value of inequality constrains! \n");
            mexErrMsgTxt(om);
        }

        if (ilc == 0 && nic > 0)
        {
            ibl = (solnp_float *)solnp_malloc(sizeof(solnp_float) * nic);
            for (i = 0; i < nic; i++)
            {
                ibl[i] = -INFINITY;
            }
        }

        if (iuc == 0 && nic > 0)
        {
            ibu = (solnp_float *)solnp_malloc(sizeof(solnp_float) * nic);
            for (i = 0; i < nic; i++)
            {
                ibu[i] = INFINITY;
            }
        }
    }

    if (ic)
    {
        for (i = 0; i < nic; i++)
        {
            if (ib0[i] <= ibl[i] || ib0[i] >= ibu[i])
            {
                sprintf(om + strlen(om), "Initial inequalities must be within the bounds \n");
                mexErrMsgTxt(om);
            }
        }
    }

    if (Ipb[0] + nic >= 0.5)
    {
        Ipb[1] = 1;
    }

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

    solnp_free(pbl);
    solnp_free(pbu);
    if (nic > 0)
    {
        solnp_free(ibl);
        solnp_free(ibu);
    }
    solnp_free(Ipc);
    solnp_free(Ipb);

    SOLNPSettings *stgs = (SOLNPSettings *)solnp_malloc(sizeof(SOLNPSettings));
    SOLNP(set_default_settings)
    (stgs, np);
    stgs->maxfev = 500 * np;
    op = prhs[1];

    tmp = mxGetField(op, 0, "rho");
    if (tmp != SOLNP_NULL)
    {
        stgs->rho = (solnp_float)*mxGetPr(tmp);
    }

    tmp = mxGetField(op, 0, "rescue");
    if (tmp != SOLNP_NULL)
    {
        stgs->rescue = (solnp_int)*mxGetPr(tmp);
    }
    else
    {
        if (nec >= np)
        {
            stgs->rescue = 1;
        }
    }

    tmp = mxGetField(op, 0, "pen_l1");
    if (tmp != SOLNP_NULL)
    {
        stgs->pen_l1 = (solnp_float)*mxGetPr(tmp);
    }

    tmp = mxGetField(op, 0, "max_iter");
    if (tmp != SOLNP_NULL)
    {
        stgs->max_iter = (solnp_int)*mxGetPr(tmp);
    }

    tmp = mxGetField(op, 0, "min_iter");
    if (tmp != SOLNP_NULL)
    {
        stgs->min_iter = (solnp_int)*mxGetPr(tmp);
    }

    tmp = mxGetField(op, 0, "max_iter_rescue");
    if (tmp != SOLNP_NULL)
    {
        stgs->max_iter_rescue = (solnp_int)*mxGetPr(tmp);
    }
    else
    {
        stgs->max_iter_rescue = stgs->max_iter;
    }

    tmp = mxGetField(op, 0, "min_iter_rescue");
    if (tmp != SOLNP_NULL)
    {
        stgs->min_iter_rescue = (solnp_int)*mxGetPr(tmp);
    }
    else
    {
        stgs->min_iter_rescue = stgs->min_iter;
    }

    tmp = mxGetField(op, 0, "delta");
    if (tmp != SOLNP_NULL)
    {
        stgs->delta = (solnp_float)*mxGetPr(tmp);
    }

    tmp = mxGetField(op, 0, "step_ratio");
    if (tmp != SOLNP_NULL)
    {
        stgs->step_ratio = (solnp_float)*mxGetPr(tmp);
    }

    tmp = mxGetField(op, 0, "max_fev");
    if (tmp != SOLNP_NULL)
    {
        stgs->maxfev = (solnp_float)*mxGetPr(tmp);
    }
    tmp = mxGetField(op, 0, "tol");
    if (tmp != SOLNP_NULL)
    {
        stgs->tol = (solnp_float)*mxGetPr(tmp);
    }

    tmp = mxGetField(op, 0, "gd_step");
    if (tmp != SOLNP_NULL)
    {
        stgs->gd_step = (solnp_float)*mxGetPr(tmp);
    }

    tmp = mxGetField(op, 0, "tol_con");
    if (tmp != SOLNP_NULL)
    {
        stgs->tol_con = (solnp_float)*mxGetPr(tmp);
    }

    tmp = mxGetField(op, 0, "ls_time");
    if (tmp != SOLNP_NULL)
    {
        stgs->ls_time = (solnp_int)*mxGetPr(tmp);
    }

    tmp = mxGetField(op, 0, "batchsize");
    if (tmp != SOLNP_NULL)
    {
        stgs->batchsize = (solnp_int)*mxGetPr(tmp);
    }
    else
    {
        stgs->batchsize = MAX(MIN(50, np / 4), 1);
    }

    tmp = mxGetField(op, 0, "tol_restart");
    if (tmp != SOLNP_NULL)
    {
        stgs->tol_restart = (solnp_float)*mxGetPr(tmp);
    }

    tmp = mxGetField(op, 0, "re_time");
    if (tmp != SOLNP_NULL)
    {
        stgs->re_time = (solnp_int)*mxGetPr(tmp);
    }

    tmp = mxGetField(op, 0, "delta_end");
    if (tmp != SOLNP_NULL)
    {
        stgs->delta_end = (solnp_float)*mxGetPr(tmp);
    }

    tmp = mxGetField(op, 0, "maxfev");
    if (tmp != SOLNP_NULL)
    {
        stgs->maxfev = (solnp_int)*mxGetPr(tmp);
    }
    tmp = mxGetField(op, 0, "noise");
    if (tmp != SOLNP_NULL)
    {
        stgs->noise = (solnp_int)*mxGetPr(tmp);
    }
    if ((!stgs->noise) && mxGetField(op, 0, "delta") == SOLNP_NULL)
    {
        stgs->delta = 1e-5;
    }
    tmp = mxGetField(op, 0, "qpsolver");
    if (tmp != SOLNP_NULL)
    {
        stgs->qpsolver = (solnp_int)*mxGetPr(tmp);
    }

    tmp = mxGetField(op, 0, "scale");
    if (tmp != SOLNP_NULL)
    {
        stgs->scale = (solnp_int)*mxGetPr(tmp);
    }

    tmp = mxGetField(op, 0, "bfgs");
    if (tmp != SOLNP_NULL)
    {
        stgs->bfgs = (solnp_int)*mxGetPr(tmp);
    }

    tmp = mxGetField(op, 0, "rs");
    if (tmp != SOLNP_NULL)
    {
        stgs->rs = (solnp_int)*mxGetPr(tmp);
        stgs->grad = 1;
    }

    tmp = mxGetField(op, 0, "drsom");
    if (tmp != SOLNP_NULL)
    {
        stgs->drsom = (solnp_int)*mxGetPr(tmp);
    }

    tmp = mxGetField(op, 0, "cen_diff");
    if (tmp != SOLNP_NULL)
    {
        stgs->cen_diff = (solnp_int)*mxGetPr(tmp);
    }

    tmp = mxGetField(op, 0, "k_i");
    if (tmp != SOLNP_NULL)
    {
        stgs->k_i = (solnp_float)*mxGetPr(tmp);
    }
    tmp = mxGetField(op, 0, "k_r");
    if (tmp != SOLNP_NULL)
    {
        stgs->k_r = (solnp_float)*mxGetPr(tmp);
    }
    tmp = mxGetField(op, 0, "c_r");
    if (tmp != SOLNP_NULL)
    {
        stgs->c_r = (solnp_float)*mxGetPr(tmp);
    }
    tmp = mxGetField(op, 0, "c_i");
    if (tmp != SOLNP_NULL)
    {
        stgs->c_i = (solnp_float)*mxGetPr(tmp);
    }
    tmp = mxGetField(op, 0, "ls_way");
    if (tmp != SOLNP_NULL)
    {
        stgs->ls_way = (solnp_int)*mxGetPr(tmp);
    }

    if (stgs->rescue)
    {
        stgs->min_iter = 1;
    }

    fun = prhs[2];
    tmp = mxGetField(fun, 0, "cost");
    if (tmp != SOLNP_NULL)
    {
        mexCallMatlabCost(SOLNP_NULL, SOLNP_NULL, 0, 0, 1, tmp);
    }
    else
    {
        sprintf(om + strlen(om), "The user does not provided cost function in the fun structure \n");
        mexErrMsgTxt(om);
    }

    if (stgs->noise == 3 && stgs->step_ratio <= 0)
    {
        sprintf(om + strlen(om), "The step_ratio must be a positive number! \n");
        mexErrMsgTxt(om);
    }

    // first call of cost function to get nec
    mxArray *lhs, *rhs[3];
    solnp_int m;
    solnp_float *ob;

    rhs[0] = tmp;
    rhs[1] = mxCreateDoubleMatrix(np, 1, mxREAL);
    rhs[2] = mxCreateDoubleMatrix(1, 1, mxREAL);

    for (i = 0; i < np; i++)
    {
        mxGetPr(rhs[1])[i] = p[i];
    }
    *mxGetPr(rhs[2]) = 1;

    lhs = SOLNP_NULL;

    SOLNP(timer)
    cost_timer;
    SOLNP(tic)
    (&cost_timer);
    mexCallMATLAB(1, &lhs, 3, rhs, "feval");
    cost_time += SOLNP(tocq)(&cost_timer) / 1e3;

    m = mxGetM(lhs);
    if (m == 1)
    {
        m = mxGetN(lhs);
    }

    tmp = mxGetField(fun, 0, "grad");
    if (tmp != SOLNP_NULL)
    {
        mexCallMatlabGradient(SOLNP_NULL, SOLNP_NULL, 0, 0, 1, tmp);
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

    tmp = mxGetField(fun, 0, "hess");
    if (tmp != SOLNP_NULL)
    {
        mexCallMatlabHessian(SOLNP_NULL, SOLNP_NULL, 0, 0, 1, tmp);
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

    ob = (solnp_float *)solnp_malloc(sizeof(solnp_float) * m);
    for (i = 0; i < m; i++)
    {
        ob[i] = (solnp_float)mxGetPr(lhs)[i];
    }
    nec = m - 1 - nic;
    nc = m - 1;

    mxDestroyArray(rhs[1]);
    mxDestroyArray(rhs[2]);
    if (lhs != SOLNP_NULL)
    {
        mxDestroyArray(lhs);
    }

    if (nrhs > 3)
    {
        lmex = prhs[3];
    }
    if (lmex && !mxIsEmpty(lmex))
    {
        m = mxGetM(lmex);
        if (m == 1)
        {
            m = mxGetN(lmex);
        }
        l = (solnp_float *)solnp_malloc(sizeof(solnp_float) * m);

        for (i = 0; i < m; i++)
        {
            l[i] = mxGetPr(lmex)[i];
        }
    }
    else
    {
        if (nc > 0)
        {
            l = (solnp_float *)solnp_malloc(sizeof(solnp_float) * nc);
            for (i = 0; i < nc; i++)
            {
                l[i] = 0;
            }
        }
        else
        {
            l = (solnp_float *)solnp_malloc(sizeof(solnp_float) * 1);
            l[0] = 0;
        }
    }

    if (nrhs > 4)
    {
        hmex = prhs[4];
    }
    if (hmex && !mxIsEmpty(hmex))
    {

        h = (solnp_float *)solnp_malloc(sizeof(solnp_float) * mxGetM(hmex) * mxGetN(hmex));

        for (i = 0; i < mxGetM(hmex) * mxGetN(hmex); i++)
        {
            h[i] = mxGetPr(hmex)[i];
        }
    }
    else
    {
        if (stgs->bfgs)
        {
            h = (solnp_float *)solnp_malloc(sizeof(solnp_float) * (np + nic) * (np + nic));
            memset(h, 0, sizeof(solnp_float) * (np + nic) * (np + nic));

            for (i = 0; i < np + nic; i++)
            {
                h[i * (np + nic) + i] = 1;
            }
        }
        else
        {
            h = SOLNP_NULL;
        }
    }

    SOLNPCost *cost = malloc_cost(nec, nic, &mexCallMatlabCost, &mexCallMatlabGradient, &mexCallMatlabHessian);
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

    SOLNPIput *input = (SOLNPIput *)solnp_malloc(sizeof(SOLNPIput));
    input->cnstr = constraint;
    input->stgs = stgs;
    input->l = l;
    input->h = h;
    input->n = np;

    SOLNPSol *sol = (SOLNPSol *)solnp_malloc(sizeof(SOLNPSol));
    SOLNPInfo *info = (SOLNPInfo *)solnp_malloc(sizeof(SOLNPInfo));

    solnp_float *ib0_p = (solnp_float *)solnp_malloc(sizeof(solnp_float) * (np + nic));

    if (nic > 0)
    {
        memcpy(ib0_p, ib0, nic * sizeof(solnp_float));
        solnp_free(ib0);
    }

    memcpy(&ib0_p[nic], p, np * sizeof(solnp_float));

    solnp_free(p);

    info->total_time = 0;

    // solnp_printf("input:\n");
    // solnp_printf("input->cnstr:\n");
    // if (nic > 0)
    // {
    //     solnp_printf("input->cnstr->il:\n");

    //     for (
    //         int i = 0;
    //         i < nic;
    //         i++)
    //     {
    //         solnp_printf("%f ", input->cnstr->il[i]);
    //     }
    //     solnp_printf("\n");

    //     solnp_printf("input->cnstr->iu:\n");

    //     for (
    //         int i = 0;
    //         i < nic;
    //         i++)
    //     {
    //         solnp_printf("%f ", input->cnstr->iu[i]);
    //     }
    //     solnp_printf("\n");

    //     // memcpy(constraint->il, ibl, nic * sizeof(solnp_float));
    //     // memcpy(constraint->iu, ibu, nic * sizeof(solnp_float));
    // }
    // else
    // {
    //     solnp_printf("nic = 0\n");
    // }

    // if (np > 0)
    // {
    //     solnp_printf("input->cnstr->pl:\n");

    //     for (
    //         int i = 0;
    //         i < np;
    //         i++)
    //     {
    //         solnp_printf("%f ", input->cnstr->pl[i]);
    //     }
    //     solnp_printf("\n");

    //     solnp_printf("input->cnstr->pu:\n");

    //     for (
    //         int i = 0;
    //         i < np;
    //         i++)
    //     {
    //         solnp_printf("%f ", input->cnstr->pu[i]);
    //     }
    //     solnp_printf("\n");

    //     // memcpy(constraint->pl, pbl, np*sizeof(solnp_float));
    //     // memcpy(constraint->pu, pbu, np*sizeof(solnp_float));
    // }
    // else
    // {
    //     solnp_printf("np = 0\n");
    // }

    // solnp_printf("input->cnstr->Ipc:\n");
    // for (int i = 0; i < 2; i++)
    // {
    //     solnp_printf("%f ", input->cnstr->Ipc[i]);
    // }
    // solnp_printf("\n");

    // solnp_printf("input->cnstr->Ipb:\n");
    // for (int i = 0; i < 2; i++)
    // {
    //     solnp_printf("%f ", input->cnstr->Ipb[i]);
    // }
    // solnp_printf("\n");

    // // memcpy(constraint->Ipc, Ipc, 2*sizeof(solnp_float));
    // // memcpy(constraint->Ipb, Ipb, 2*sizeof(solnp_float));

    // solnp_printf("ib0_p:\n");
    // for (int i = 0; i < np + nic; i++)
    // {
    //     if (i == 358 || i == np + nic - 163)
    //     {
    //         solnp_printf("\n");
    //     }
    //     solnp_printf("%f ", ib0_p[i]);
    // }
    // solnp_printf("\n");

    // solnp_printf("input->stgs:\n");
    // solnp_printf("pen_l1: %f\n", input->stgs->pen_l1);
    // solnp_printf("rho: %f\n", input->stgs->rho);
    // solnp_printf("max_iter: %d\n", input->stgs->max_iter);
    // solnp_printf("min_iter: %d\n", input->stgs->min_iter);
    // solnp_printf("max_iter_rescue: %d\n", input->stgs->max_iter_rescue);
    // solnp_printf("min_iter_rescue: %d\n", input->stgs->min_iter_rescue);
    // solnp_printf("delta: %f\n", input->stgs->delta);
    // solnp_printf("tol: %f\n", input->stgs->tol);
    // solnp_printf("tol_con: %f\n", input->stgs->tol_con);
    // solnp_printf("ls_time: %d\n", input->stgs->ls_time);
    // solnp_printf("tol_restart: %f\n", input->stgs->tol_restart);
    // solnp_printf("re_time: %d\n", input->stgs->re_time);
    // solnp_printf("delta_end: %f\n", input->stgs->delta_end);
    // solnp_printf("maxfev: %d\n", input->stgs->maxfev);
    // solnp_printf("noise: %d\n", input->stgs->noise);
    // solnp_printf("qpsolver: %d\n", input->stgs->qpsolver);
    // solnp_printf("k_r: %f\n", input->stgs->k_r);
    // solnp_printf("k_i: %f\n", input->stgs->k_i);
    // solnp_printf("c_r: %f\n", input->stgs->c_r);
    // solnp_printf("c_i: %f\n", input->stgs->c_i);
    // solnp_printf("batchsize: %d\n", input->stgs->batchsize);
    // solnp_printf("hess: %d\n", input->stgs->hess);
    // solnp_printf("grad: %d\n", input->stgs->grad);
    // solnp_printf("rescue: %d\n", input->stgs->rescue);
    // solnp_printf("ls_way: %d\n", input->stgs->ls_way);
    // solnp_printf("bfgs: %d\n", input->stgs->bfgs);
    // solnp_printf("rs: %d\n", input->stgs->rs);
    // solnp_printf("scale: %d\n", input->stgs->scale);
    // solnp_printf("drsom: %d\n", input->stgs->drsom);
    // solnp_printf("\n");

    // solnp_printf("\nnec: %d", nec);
    // solnp_printf("nic: %d\n", nic);
    // solnp_printf("Ipb[0]: %f\n", constraint->Ipb[0]);
    // solnp_printf("Ipb[1]: %f\n", constraint->Ipb[1]);
    // solnp_printf("Ipc[0]: %f\n", constraint->Ipc[0]);
    // solnp_printf("Ipc[1]: %f\n", constraint->Ipc[1]);

    solnp_int status = SOLNP(main)(input, cost, ib0_p, sol, info);

    info->total_time += SOLNP(tocq)(&total_timer) / 1e3;
    info->cost_time = cost_time;
    cost_time = 0;

    printf("total time is:%e\n", info->total_time);
    // printf("\ntime for calling q  p solver is:%f(%.2f%%)", info->qpsolver_time, info->qpsolver_time/info->total_time * 100);
    // printf("\ntime for calculating cost function is:%f(%.2f%%)\n", info->cost_time, info->cost_time / info->total_time * 100);

    solnp_free(ib0_p);

    plhs[0] = mxCreateStructArray(1, one, num_sol_fields, sol_fields);

    tmp = mxCreateDoubleMatrix(1, 1, mxREAL);
    mxSetField(plhs[0], 0, "iter", tmp);
    *mxGetPr(tmp) = (solnp_float)sol->iter;

    tmp = mxCreateDoubleMatrix(1, 1, mxREAL);
    mxSetField(plhs[0], 0, "count_cost", tmp);
    *mxGetPr(tmp) = (solnp_float)(sol->count_cost);

    tmp = mxCreateDoubleMatrix(1, 1, mxREAL);
    mxSetField(plhs[0], 0, "count_grad", tmp);
    *mxGetPr(tmp) = (solnp_float)(sol->count_grad);

    tmp = mxCreateDoubleMatrix(1, 1, mxREAL);
    mxSetField(plhs[0], 0, "count_hess", tmp);
    *mxGetPr(tmp) = (solnp_float)(sol->count_hess);

    tmp = mxCreateDoubleMatrix(1, 1, mxREAL);
    mxSetField(plhs[0], 0, "constraint", tmp);
    *mxGetPr(tmp) = sol->constraint;

    tmp = mxCreateDoubleMatrix(1, 1, mxREAL);
    mxSetField(plhs[0], 0, "restart_time", tmp);
    *mxGetPr(tmp) = sol->restart_time;

    tmp = mxCreateDoubleMatrix(1, 1, mxREAL);
    mxSetField(plhs[0], 0, "obj", tmp);
    *mxGetPr(tmp) = sol->obj;

    tmp = mxCreateDoubleMatrix(1, 1, mxREAL);
    mxSetField(plhs[0], 0, "status", tmp);
    *mxGetPr(tmp) = (solnp_float)sol->status;

    tmp = mxCreateDoubleMatrix(1, 1, mxREAL);
    mxSetField(plhs[0], 0, "solve_time", tmp);
    *mxGetPr(tmp) = (solnp_float)info->total_time;

    set_output_field(&tmp, sol->p, np, 1);
    mxSetField(plhs[0], 0, "p", tmp);

    set_output_field(&tmp, sol->best_fea_p, np, 1);
    mxSetField(plhs[0], 0, "best_fea_p", tmp);

    set_output_field(&tmp, sol->ic, MAX(nic, 1), 1);
    mxSetField(plhs[0], 0, "ic", tmp);

    set_output_field(&tmp, sol->count_h, sol->iter + 1, 1);
    mxSetField(plhs[0], 0, "cost_his", tmp);

    set_output_field(&tmp, sol->jh, sol->iter + 1, 1);
    mxSetField(plhs[0], 0, "jh", tmp);

    set_output_field(&tmp, sol->ch, sol->iter + 1, 1);
    mxSetField(plhs[0], 0, "ch", tmp);

    set_output_field(&tmp, sol->l, MAX(nc, 1), 1);
    mxSetField(plhs[0], 0, "l", tmp);

    if (stgs->bfgs)
    {
        set_output_field(&tmp, sol->h, np + nic, np + nic);
    }
    else
    {
        set_output_field(&tmp, sol->h, 1, 1);
    }
    mxSetField(plhs[0], 0, "h", tmp);

    solnp_free(info);
    SOLNP(free_sol)
    (sol);

    // close function handle
    mexCallMatlabCost(SOLNP_NULL, SOLNP_NULL, 0, 0, -1, SOLNP_NULL);
    mexCallMatlabGradient(SOLNP_NULL, SOLNP_NULL, 0, 0, -1, SOLNP_NULL);
    mexCallMatlabHessian(SOLNP_NULL, SOLNP_NULL, 0, 0, -1, SOLNP_NULL);
}

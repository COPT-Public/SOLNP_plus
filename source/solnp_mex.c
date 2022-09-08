#include "solnp.h"
#include "matrix.h"
#include "mex.h"
#include "linalg.h"
#include "solnp_util.h"

solnp_float cost_time = 0; // global cost timer

#if def !(DLONG > 0)
/* this memory must be freed */
static solnp_int *cast_to_solnp_int_arr(mwIndex *arr, solnp_int len) {
    solnp_int i;
    solnp_int *arr_out = (solnp_int *)solnp_malloc(sizeof(solnp_int) * len);
    for (i = 0; i < len; i++) {
        arr_out[i] = (solnp_int)arr[i];
    }
    return arr_out;
}
#endif

#if SFLOAT > 0
/* this memory must be freed */
static solnp_float *cast_to_solnp_float_arr(double *arr, solnp_int len) {
    solnp_int i;
    solnp_float *arr_out = (solnp_float *)solnp_malloc(sizeof(solnp_float) * len);
    for (i = 0; i < len; i++) {
        arr_out[i] = (solnp_float)arr[i];
    }
    return arr_out;
}

static double *cast_to_double_arr(solnp_float *arr, solnp_int len) {
    solnp_int i;
    double *arr_out = (double *)solnp_malloc(sizeof(double) * len);
    for (i = 0; i < len; i++) {
        arr_out[i] = (double)arr[i];
    }
    return arr_out;
}
#endif


static void set_output_field(mxArray **pout, solnp_float *out, solnp_int len) {
    *pout = mxCreateDoubleMatrix(0, 0, mxREAL);
#if SFLOAT > 0
    mxSetPr(*pout, cast_to_double_arr(out, len));
    solnp_free(out);
#else
    mxSetPr(*pout, out);
#endif
    mxSetM(*pout, len);
    mxSetN(*pout, 1);
}


void mexCallMatlabCost(SOLNPCost *c, solnp_float *p, solnp_int np, solnp_int it){


    SOLNP(timer) cost_timer;
    SOLNP(tic)(&cost_timer);

    mxArray *lhs, *rhs[2];
    solnp_int i;

    rhs[0] =  mxCreateDoubleMatrix(np, 1, mxREAL);

    for(i=0; i<np; i++){
        mxGetPr(rhs[0])[i] = (double)p[i];
    }

    rhs[1] =  mxCreateDoubleMatrix(1, 1, mxREAL);
    *mxGetPr(rhs[1]) = (int)it;

    lhs = SOLNP_NULL;

    mexCallMATLAB(1, &lhs, 1, rhs, "cost");

    // double *result = mxGetPr(lhs);

    c->obj = (solnp_float)mxGetPr(lhs)[0];

    for(i = 0; i<c->nec; i++){
        c->ec[i] = (solnp_float)mxGetPr(lhs)[i+1];
    }
    for(i = 0; i<c->nic; i++){
        c->ic[i] = (solnp_float)mxGetPr(lhs)[i+1+c->nec];
    }

    mxDestroyArray(rhs[0]);
    mxDestroyArray(rhs[1]);

    if(lhs != SOLNP_NULL){
        mxDestroyArray(lhs);
    }

    cost_time += SOLNP(tocq)(&cost_timer) / 1e3;

}

// mex file for matlab function: [p,jh,l,h,ic]=solnp(cnstr,op,l,h)
void  mexFunction(int  nlhs, mxArray* plhs[], int  nrhs, const  mxArray* prhs[])
{
    SOLNP(timer) total_timer;
    SOLNP(tic)(&total_timer);

    if(nrhs < 1){
        mexErrMsgTxt("Syntax error");
    }

    const mwSize one[1] = {1};
    const int num_sol_fields = 11;
    const char *sol_fields[] = {"p", "jh","ch", "l", "h", "ic", "iter", "count_cost","constraint","obj","status"};

    const mxArray *cnstr;
    const mxArray *op;

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


    solnp_float *Ipc = (solnp_float*)solnp_malloc(2 * sizeof(solnp_float));
    memset(Ipc, 0, 2 * sizeof(solnp_float));
    solnp_float *Ipb = (solnp_float*)solnp_malloc(2 * sizeof(solnp_float));
    memset(Ipb, 0, 2 * sizeof(solnp_float));

    char* om = (char*)solnp_malloc(sizeof(char) * 512); // error message
    sprintf(om, "SOLNP--> ");


    cnstr = prhs[0];

    pblmx = (mxArray *)mxGetField(cnstr, 0, "pbl");
    if (pblmx && !mxIsEmpty(pblmx)){
        Ipc[0] = 1;
        Ipb[0] = 1;
		ns = (solnp_int)mxGetNumberOfDimensions(pblmx);
		dims = mxGetDimensions(pblmx);
        np = (solnp_int)dims[0];

        if (ns > 1 && dims[0] == 1) {
		    np = (solnp_int)dims[1];
		}

        pbl = (solnp_float *)solnp_malloc(sizeof(solnp_float) * np);
		for (i = 0; i < np; i++) {
		    pbl[i] = (solnp_float)mxGetPr(pblmx)[i];
		}
    }


    pbumx = (mxArray *)mxGetField(cnstr, 0, "pbu");
    if (pbumx && !mxIsEmpty(pbumx)){
        Ipc[1] = 1;
        Ipb[0] = 1;
		ns = (solnp_int)mxGetNumberOfDimensions(pbumx);
		dims = mxGetDimensions(pbumx);
        np = (solnp_int)dims[0];

        if (ns > 1 && dims[0] == 1) {
		    np = (solnp_int)dims[1];
		}

        pbu = (solnp_float *)solnp_malloc(sizeof(solnp_float) * np);
		for (i = 0; i < np; i++) {
		    pbu[i] = (solnp_float)mxGetPr(pbumx)[i];
		}
    }

    p0mx = (mxArray *)mxGetField(cnstr, 0, "p0");
    if (p0mx && !mxIsEmpty(p0mx)){

		ns = (solnp_int)mxGetNumberOfDimensions(p0mx);
		dims = mxGetDimensions(p0mx);
        np = (solnp_int)dims[0];

        if (ns > 1 && dims[0] == 1) {
		    np = (solnp_int)dims[1];
		}

        p = (solnp_float *)solnp_malloc(sizeof(solnp_float) * np);
		for (i = 0; i < np; i++) {
		    p[i] = (solnp_float)mxGetPr(p0mx)[i];
		}
    }

    if(Ipc[0] && Ipc[1]){

        if(!(p0mx && !mxIsEmpty(p0mx))){
            p = (solnp_float *)solnp_malloc(sizeof(solnp_float) * np);
            memcpy(p, pbl, sizeof(solnp_float) * np);
            SOLNP(add_scaled_array)(p, pbu, np, 1);
            SOLNP(scale_array)(p, 0.5, np);
        }
        
        for(i = 0; i < np; i++){
            if(mxIsInf(p[i])){
                sprintf(om + strlen(om), "The user does not provide initiate point! \n");
                mexErrMsgTxt(om);
            }
        }
    }
    else{
        if(!(p0mx && !mxIsEmpty(p0mx))){
            sprintf(om + strlen(om), "The user does not provide initiate point! \n");
            mexErrMsgTxt(om);
        }
        if(Ipc[0] == 0){

            pbl = (solnp_float *)solnp_malloc(sizeof(solnp_float) * np);
            for(i = 0; i < np; i++){
                pbl[i] = -INFINITY;
            }
        }
        if(Ipc[1] == 0){
            pbu = (solnp_float *)solnp_malloc(sizeof(solnp_float) * np);
            for(i = 0; i < np; i++){
                pbu[i] = INFINITY;
            }
        }
    }

    if(Ipb[0] > 0.5){

        for(i = 0; i < np; i++){
            if(pbu[i] <= pbl[i]){
                sprintf(om + strlen(om), "The lower bounds of the parameter constraints \n");
                sprintf(om + strlen(om), "          must be strictly less than the upper bounds. \n");
                mexErrMsgTxt(om);
            }
            else if(p[i] <= pbl[i] || p[i] >= pbu[i]){
                sprintf(om + strlen(om), "Initial parameter values must be within the bounds \n");
                mexErrMsgTxt(om);
            }
        }
    }


    iblmx = (mxArray *)mxGetField(cnstr, 0, "ibl");
    if (iblmx && !mxIsEmpty(iblmx)){
        ilc = 1;
		ns = (solnp_int)mxGetNumberOfDimensions(iblmx);
		dims = mxGetDimensions(iblmx);
        nic = (solnp_int)dims[0];

        if (ns > 1 && dims[0] == 1) {
		    nic = (solnp_int)dims[1];
		}

        ibl = (solnp_float *)solnp_malloc(sizeof(solnp_float) * nic);
		for (i = 0; i < nic; i++) {
		    ibl[i] = (solnp_float)mxGetPr(iblmx)[i];
		}
    }
    
    ibumx = (mxArray *)mxGetField(cnstr, 0, "ibu");
    if (ibumx && !mxIsEmpty(ibumx)){
        iuc = 1;
		ns = (solnp_int)mxGetNumberOfDimensions(ibumx);
		dims = mxGetDimensions(ibumx);
        nic = (solnp_int)dims[0];

        if (ns > 1 && dims[0] == 1) {
		    nic = (solnp_int)dims[1];
		}

        ibu = (solnp_float *)solnp_malloc(sizeof(solnp_float) * nic);
		for (i = 0; i < nic; i++) {
		    ibu[i] = (solnp_float)mxGetPr(ibumx)[i];
		}
    }

    ib0mx = (mxArray *)mxGetField(cnstr, 0, "ib0");
    if (ib0mx && !mxIsEmpty(ib0mx)){
		ns = (solnp_int)mxGetNumberOfDimensions(ib0mx);
		dims = mxGetDimensions(ib0mx);
        nic = (solnp_int)dims[0];

        if (ns > 1 && dims[0] == 1) {
		    nic = (solnp_int)dims[1];
		}

        ib0 = (solnp_float *)solnp_malloc(sizeof(solnp_float) * nic);
		for (i = 0; i < nic; i++) {
		    ib0[i] = (solnp_float)mxGetPr(ib0mx)[i];
		}
    }

    ic = ilc || iuc;

    if(ilc && iuc){

        if(!(ib0mx && !mxIsEmpty(ib0mx))){

            ib0 = (solnp_float *)solnp_malloc(sizeof(solnp_float) * nic);
            memcpy(ib0, ibl, sizeof(solnp_float) * nic);
            SOLNP(add_scaled_array)(ib0, ibu, nic, 1);
            SOLNP(scale_array)(ib0, 0.5, nic);

            for(i = 0; i < nic; i++){
                if(mxIsInf(ib0[i])){
                    sprintf(om + strlen(om), "The user does not provided initiate value of inequality constrains! \n");
                    mexErrMsgTxt(om);
                }
            }
        }

        for(i=0; i<nic; i++){
            if(ibu[i] <= ibl[i]){
                sprintf(om + strlen(om), "The lower bounds of the inequality constraints \n");
                sprintf(om + strlen(om), "          must be strictly less than the upper bounds. \n");
                mexErrMsgTxt(om);
            }
        }
    }

    else if(ic == 1){

        if(!(ib0mx && !mxIsEmpty(ib0mx))){
            sprintf(om + strlen(om), "The user does not provided initiate value of inequality constrains! \n");
            mexErrMsgTxt(om);
        }

        if(ilc == 0 && nic > 0){
            ibl = (solnp_float *)solnp_malloc(sizeof(solnp_float) * nic);
            for(i=0; i<nic; i++){
                ibl[i] = -INFINITY;
            }
        }

        if(iuc == 0 && nic > 0){
            ibu = (solnp_float *)solnp_malloc(sizeof(solnp_float) * nic);
            for(i=0; i<nic; i++){
                ibu[i] = INFINITY;
            }
        }
    }

    if(ic){
        for(i = 0; i<nic; i++){
            if(ib0[i] <= ibl[i] || ib0[i] >= ibu[i]){
                sprintf(om + strlen(om), "Initial inequalities must be within the bounds \n");
                mexErrMsgTxt(om);
            }
        }
    }

    if(Ipb[0] + nic >= 0.5){
        Ipb[1] = 1;
    }

    SOLNPConstraint *constraint = malloc_constriant(np, nic);
    if (nic > 0) {
        memcpy(constraint->il, ibl, nic * sizeof(solnp_float));
        memcpy(constraint->iu, ibu, nic * sizeof(solnp_float));
    }

    if(np > 0){
        memcpy(constraint->pl, pbl, np*sizeof(solnp_float));
        memcpy(constraint->pu, pbu, np*sizeof(solnp_float));
    } 


    memcpy(constraint->Ipc, Ipc, 2*sizeof(solnp_float));
    memcpy(constraint->Ipb, Ipb, 2*sizeof(solnp_float));


    solnp_free(pbl);
    solnp_free(pbu);
    if(nic > 0){
        solnp_free(ibl);
        solnp_free(ibu);
    }
    solnp_free(Ipc);
    solnp_free(Ipb);


    SOLNPSettings *stgs = (SOLNPSettings *)solnp_malloc(sizeof(SOLNPSettings));
    set_default_settings(stgs);
    stgs->maxfev = 500 * np;
    op = prhs[1];

    tmp = mxGetField(op, 0, "rho");
    if (tmp != SOLNP_NULL) {
        stgs->rho = (solnp_float)*mxGetPr(tmp);
    }

    tmp = mxGetField(op, 0, "max_iter");
    if (tmp != SOLNP_NULL) {
        stgs->max_iter = (solnp_int)*mxGetPr(tmp);
    }

    tmp = mxGetField(op, 0, "min_iter");
    if (tmp != SOLNP_NULL) {
        stgs->min_iter = (solnp_int)*mxGetPr(tmp);
    }

    tmp = mxGetField(op, 0, "delta");
    if (tmp != SOLNP_NULL) {
        stgs->delta = (solnp_float)*mxGetPr(tmp);
    }

    tmp = mxGetField(op, 0, "tol");
    if (tmp != SOLNP_NULL) {
        stgs->tol = (solnp_float)*mxGetPr(tmp);
    }

    tmp = mxGetField(op, 0, "tol_con");
    if (tmp != SOLNP_NULL) {
        stgs->tol_con = (solnp_float)*mxGetPr(tmp);
    }

    tmp = mxGetField(op, 0, "ls_time");
    if (tmp != SOLNP_NULL) {
        stgs->ls_time = (solnp_int)*mxGetPr(tmp);
    }

    tmp = mxGetField(op, 0, "tol_restart");
    if (tmp != SOLNP_NULL) {
        stgs->tol_restart = (solnp_float)*mxGetPr(tmp);
    }

    tmp = mxGetField(op, 0, "re_time");
    if (tmp != SOLNP_NULL) {
        stgs->re_time = (solnp_int)*mxGetPr(tmp);
    }

    tmp = mxGetField(op, 0, "delta_end");
    if (tmp != SOLNP_NULL) {
        stgs->delta_end = (solnp_float)*mxGetPr(tmp);
    }

    tmp = mxGetField(op, 0, "maxfev");
    if (tmp != SOLNP_NULL) {
        stgs->maxfev = (solnp_int)*mxGetPr(tmp);
    }
    tmp = mxGetField(op, 0, "noise");
    if (tmp != SOLNP_NULL) {
        stgs->noise = (solnp_int)*mxGetPr(tmp);
    }
    tmp = mxGetField(op, 0, "qpsolver");
    if (tmp != SOLNP_NULL) {
        stgs->qpsolver = (solnp_int)*mxGetPr(tmp);
    }
    if ((!stgs->noise) && stgs->delta == 1.) {
        stgs->delta = 1e-5;
    }

    // first call of cost function to get nec
    mxArray *lhs, *rhs;
    solnp_int m;
    solnp_float *ob;

    rhs =  mxCreateDoubleMatrix(np, 1, mxREAL);

    for(i=0; i<np; i++){
        mxGetPr(rhs)[i] = p[i];
    }


    lhs = SOLNP_NULL;


    SOLNP(timer) cost_timer;
    SOLNP(tic)(&cost_timer);
    mexCallMATLAB(1, &lhs, 1, &rhs, "cost");
    cost_time += SOLNP(tocq)(&cost_timer) / 1e3;

    m = mxGetM(lhs);
    if(m == 1){
        m = mxGetN(lhs);
    }

    ob = (solnp_float*)solnp_malloc(sizeof(solnp_float) * m);
    for(i = 0; i < m; i++){
        ob[i] = (solnp_float)mxGetPr(lhs)[i];
    }
    nec = m - 1 - nic;
    nc = m - 1;

    mxDestroyArray(rhs);
    if(lhs != SOLNP_NULL){
        mxDestroyArray(lhs);
    }


    if (nrhs > 2) {
        lmex = prhs[2];
    }
    if (lmex && !mxIsEmpty(lmex)){
        m = mxGetM(lmex);
        if(m == 1){
            m = mxGetN(lmex);
        }
        l = (solnp_float*)solnp_malloc(sizeof(solnp_float) * m);

        for(i=0; i<m; i++){
            l[i] = mxGetPr(lmex)[i];
        }
    }
    else{
        if(nc > 0){
            l = (solnp_float*)solnp_malloc(sizeof(solnp_float) * nc);
            for(i=0; i<nc; i++){
                l[i] = 0;
            }
        }
        else{
            l = (solnp_float*)solnp_malloc(sizeof(solnp_float) * 1);
            l[0] = 0;
        }
    }

    if (nrhs > 3) {
        hmex = prhs[3];
    }
    if (hmex && !mxIsEmpty(hmex)){

        h = (solnp_float*)solnp_malloc(sizeof(solnp_float) * mxGetM(hmex) * mxGetN(hmex));

        for(i=0; i<mxGetM(hmex) * mxGetN(hmex); i++){
            h[i] = mxGetPr(hmex)[i];
        }
    }
    else{

        h = (solnp_float*)solnp_malloc(sizeof(solnp_float) * (np + nic) * (np + nic));
        memset(h, 0 , sizeof(solnp_float) * (np + nic) * (np + nic));

        for(i = 0; i < np+nic; i++){
            h[i * (np+nic) + i] = 1;
        }

    }



    SOLNPCost *cost = malloc_cost(nec, nic, &mexCallMatlabCost);
    cost->obj = ob[0];
    for(i=0; i<nec; i++){
        cost->ec[i] = ob[1 + i];
    }
    for(i=0; i<nic; i++){
        cost->ic[i] = ob[1 + nec + i];
    }

    solnp_free(ob);
    

    SOLNPIput *input = (SOLNPIput*)solnp_malloc(sizeof(SOLNPIput));
    input->cnstr = constraint;
    input->stgs = stgs;
    input->l = l;
    input->h = h;
    input->n = np;

    SOLNPSol *sol = (SOLNPSol*)solnp_malloc(sizeof(SOLNPSol));
    SOLNPInfo *info = (SOLNPInfo*)solnp_malloc(sizeof(SOLNPInfo));

    solnp_float *ib0_p = (solnp_float*)solnp_malloc(sizeof(solnp_float) * (np + nic));
    
    if (nic > 0) {
        memcpy(ib0_p, ib0, nic * sizeof(solnp_float));
        solnp_free(ib0);
    }

    memcpy(&ib0_p[nic], p, np * sizeof(solnp_float));

    solnp_free(p);

    info->total_time = 0;

    solnp_int status = SOLNP(main)(input, cost, ib0_p, sol, info);

    info->total_time += SOLNP(tocq)(&total_timer) / 1e3;
    info->cost_time = cost_time;
    cost_time = 0;

    printf("total time is:%e\n", info->total_time);
   // printf("\ntime for calling q  p solver is:%f(%.2f%%)", info->qpsolver_time, info->qpsolver_time/info->total_time * 100);
   // printf("\ntime for calculating cost function is:%f(%.2f%%)\n", info->cost_time, info->cost_time / info->total_time * 100);
 
    solnp_free(info);
    solnp_free(ib0_p);

    /* output sol */
    plhs[0] = mxCreateStructArray(1, one, num_sol_fields, sol_fields);

    tmp = mxCreateDoubleMatrix(1, 1, mxREAL);
    mxSetField(plhs[0], 0, "iter", tmp);
    *mxGetPr(tmp) = (solnp_float)sol->iter;

    set_output_field(&tmp, sol->p, np);
    mxSetField(plhs[0], 0, "p", tmp);

    set_output_field(&tmp, sol->ic, MAX(nic,1));
    mxSetField(plhs[0], 0, "ic", tmp);

    set_output_field(&tmp, sol->jh, sol->iter + 1);
    mxSetField(plhs[0], 0, "jh", tmp);

    set_output_field(&tmp, sol->ch, sol->iter + 1);
    mxSetField(plhs[0], 0, "ch", tmp);

    set_output_field(&tmp, sol->l, MAX(nc,1));
    mxSetField(plhs[0], 0, "l", tmp);

    tmp = mxCreateDoubleMatrix(1, 1, mxREAL);
    mxSetField(plhs[0], 0, "count_cost", tmp);
    *mxGetPr(tmp) = (solnp_float)(sol->count_cost);
   
    tmp = mxCreateDoubleMatrix(1, 1, mxREAL);
    mxSetField(plhs[0], 0, "constraint", tmp);
    *mxGetPr(tmp) = sol->constraint;

    tmp = mxCreateDoubleMatrix(1, 1, mxREAL);
    mxSetField(plhs[0], 0, "obj", tmp);
    *mxGetPr(tmp) = sol->obj;

    tmp = mxCreateDoubleMatrix(1, 1, mxREAL);
    mxSetField(plhs[0], 0, "status", tmp);
    *mxGetPr(tmp) = (solnp_float)sol->status;


    /*
        set output h
    */
    tmp = mxCreateDoubleMatrix(0, 0, mxREAL);
    #if SFLOAT > 0
        mxSetPr(tmp, cast_to_double_arr(tmp, (np+nic)*(np+nic)));
        solnp_free(sol->h);
    #else
        mxSetPr(tmp, sol->h);
    #endif
    mxSetM(tmp, np+nic);
    mxSetN(tmp, np+nic);
    mxSetField(plhs[0], 0, "h", tmp);
}
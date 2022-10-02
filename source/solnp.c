#include "solnp.h"
#include "linalg.h"
#include "subnp.h"
#include "solnp_util.h"

SOLNPWork *init_work
(
    SOLNPIput *input,
    SOLNPCost *cost
)
{
    SOLNPWork *w = (SOLNPWork *)solnp_malloc(sizeof(SOLNPWork));
    w->n = input->n;
    w->nec = cost->nec;
    w->nic = cost->nic;
    w->nc = w->nec + w->nic;
    w->ob = cost;
    w->pb = input->cnstr;
    w->rho = input->stgs->rho;

    w->constraint = (solnp_float*)solnp_malloc(w->nc * sizeof(solnp_float));
    w->p = (solnp_float*)solnp_malloc((w->n + w->nic) * sizeof(solnp_float));
    w->l = (solnp_float*)solnp_malloc(MAX(1, w->nc) * sizeof(solnp_float));
    w->bestl = (solnp_float*)solnp_malloc(MAX(1, w->nc) * sizeof(solnp_float));
    w->h = (solnp_float*)solnp_malloc((w->nic+w->n)* (w->nic + w->n) * sizeof(solnp_float));
    w->jh = (solnp_float*)solnp_malloc((input->stgs->max_iter+1) * sizeof(solnp_float));
    w->ch = (solnp_float*)solnp_malloc((input->stgs->max_iter + 1) * sizeof(solnp_float));
    w->bestp = (solnp_float*)solnp_malloc((w->n) * sizeof(solnp_float));
    w->count_cost = 1;
    w->exit = 0;
    return w;
}


solnp_int update_work
(
    SOLNPWork *w,
    SOLNPIput *input,
    solnp_float *ib0_p
)
{   
    solnp_int nec = w->nec;
    solnp_int nic = w->nic;
    solnp_int nc = w->nc;
    solnp_int i;

    memcpy(w->p, ib0_p, (w->n + w->nic)*sizeof(solnp_float));
    w->bestobj = INFINITY;
    if(nc > 0.5){
        memcpy(w->l, input->l, nc*sizeof(solnp_float));
        memcpy(w->constraint, w->ob->ec, nec*sizeof(solnp_float));
        memcpy(&w->constraint[nec], w->ob->ic, nic*sizeof(solnp_float));

        if(nic > 0.5){
            for(i=0; i<nic; i++){
                if(w->constraint[nec+i] <= w->pb->il[i] || w->constraint[nec+i] >= w->pb->iu[i]) break;
            }
            if(i==nic) memcpy(w->p, &w->constraint[nec], nic*sizeof(solnp_float));

            SOLNP(add_scaled_array)(&w->constraint[nec], w->p, nic, -1.0);
        }

        w->cons_nm1 = SOLNP(norm)(w->constraint, w->nc);
        w->bestcon = w->cons_nm1;
        if(MAX(w->cons_nm1 - 10*input->stgs->tol, nic) <= 0){
            w->rho = 0;
        }

    }
    else{
        w->l[0] = 0;
        w->constraint = SOLNP_NULL;
    }

    w->mu = w->n;
    w->j = w->ob->obj;
    w->jh[0] = w->j;
    w->ch[0] = w->cons_nm1;
    memcpy(w->h, input->h, (w->nic+w->n)* (w->nic + w->n) * sizeof(solnp_float));
    
}



SOLNPWork *SOLNP(init)(
    SOLNPIput *input,
    SOLNPCost *cost,
    solnp_float *ib0_p
)
{
    SOLNPWork *w;
    w = init_work(input, cost);


    return w;
}


solnp_int SOLNP(solve)
(
    SOLNPWork *w,
    SOLNPIput *input,
    solnp_float *ib0_p,
    SOLNPSol *sol,
    SOLNPInfo *info
)
{
    solnp_int i = 0;
    solnp_int j;

    solnp_int n = w->n;
    solnp_int nc = w->nc;
    solnp_int nec = w->nec;
    solnp_int nic = w->nic;
    solnp_int restart = 0;
    SOLNPSettings* stgs = input->stgs;
    solnp_float delta0 = stgs->delta;
    update_work(w, input, ib0_p);

   
   
    for(i=0; i<stgs->max_iter; i++)
    {
        subnp_qp(w, stgs, info);

        //w->ob->cost(w->ob, &w->p[w->nic], n, i);
        w->obj_dif = (w->j - w->ob->obj) / MAX(ABS(w->ob->obj), 1);
        w->j = w->ob->obj;

        if(nc > 0.5){

            memcpy(w->constraint, w->ob->ec, nec*sizeof(solnp_float));
            memcpy(&w->constraint[nec], w->ob->ic, nic*sizeof(solnp_float));

            if(nic > 0.5){

                for(j=0; j<nic; j++){
                    if(w->constraint[nec+j] <= w->pb->il[j] || w->constraint[nec+j] >= w->pb->iu[j]) break;
                }
                if(j == nic) memcpy(w->p, &w->constraint[nec], nic*sizeof(solnp_float));

                SOLNP(add_scaled_array)(&w->constraint[nec], w->p, nic, -1.0);
            }

            w->cons_nm2 = SOLNP(norm)(w->constraint, nc);

            if(w->cons_nm2 < 10*stgs->tol){
                w->rho = 0;
                w->mu = MIN(w->mu, stgs->tol);
            }

            if(w->cons_nm2 < 5*w->cons_nm1){
                w->rho /= 5;
            }
            else if(w->cons_nm2 > 10*w->cons_nm1){
                w->rho = 5*MAX(w->rho, SQRTF(stgs->tol));
            }
            if(w->exit == 0 && MAX(w->obj_dif, w->cons_nm1-w->cons_nm2) <= 0 && restart < stgs->re_time && stgs->delta <= MAX(3 * stgs->tol,stgs->delta_end)){
                SOLNP(scale_array)(w->l, 0, nc);
                restart++;
                if (stgs->noise) {
                    stgs->delta = delta0;
                }
                //h=diag(diag(h));
                if (stgs->noise) {
                    stgs->delta = delta0;
                }
                for(int col=0; col<w->nic + w->n; col++){
                    for(int row=0; row<w->nic + w->n; row++){
                        if(col != row){
                            w->h[col * (w->nic + w->n) + row] = 0;
                        }
                    }
                }
            }

            w->cons_nm1 = w->cons_nm2;
        }

        w->jh[1 + i] = w->j;    
        w->ch[1 + i] = w->cons_nm1;

        if ((w->cons_nm1<= stgs->tol_con &&  w->obj_dif <= stgs->tol && stgs->delta <= MAX(stgs->tol,stgs->delta_end)) || w->exit ) {            
            if (w->exit == 0 && restart < stgs->re_time && w->alm_crit > stgs->tol_restart) {
                restart++;
                if (stgs->noise) {
                    stgs->delta = delta0;
                }
                //h=diag(diag(h));
                for (int col = 0; col < w->nic + w->n; col++) {
                    for (int row = 0; row < w->nic + w->n; row++) {
                        if (col != row) {
                            w->h[col * (w->nic + w->n) + row] = 0;
                        }
                    }
                }
            }
            else {
                i++;
                break;
            }
        }
    }

    sol->p = (solnp_float*)solnp_malloc(n * sizeof(solnp_float));
    sol->l = (solnp_float*)solnp_malloc(MAX(1, nc) * sizeof(solnp_float));
    if (w->bestobj != INFINITY) {
        memcpy(sol->p, w->bestp, n * sizeof(solnp_float));
        memcpy(sol->l, w->bestl, MAX(nc, 1) * sizeof(solnp_float));
    }
    else {
        memcpy(sol->p, &w->p[w->nic], n * sizeof(solnp_float));
        memcpy(sol->l, w->l, MAX(nc, 1) * sizeof(solnp_float));
        w->bestcon = w->cons_nm2;
        w->bestobj = w ->ob->obj;
    }


    if(w->cons_nm1 <= stgs->tol_con && w->obj_dif <= stgs->tol && stgs->delta <= MAX(stgs->tol,stgs->delta_end)){
        sol->status = 1;//Success
        printf("SOLNP+--> Success! Completed in %d iterations\n", i);
        printf("         The infeasibility is %e.\n", w->bestcon);
    }
    else {
        if (w->exit == 1) {
            sol->status = 0;
            printf("SOLNP+--> Exiting after maximum number of function evaluation. Tolerance not achieved.\n");
            printf("         The infeasibility is %e.\n", w->bestcon);
            printf("         SOLNP has restarted %d times.\n", restart);
        }else if (w->exit == 2) {
            sol->status = -3;
            printf("SOLNP+--> Exiting because of unknown error. Tolerance not achieved.\n");
            printf("         The infeasibility is %e.\n", w->bestcon);
        }else if (w->cons_nm1 > stgs->tol_con) {
            sol->status = -1;//Fail to find a feasible point. 
            printf("SOLNP+--> Exiting after maximum number of iterations. Tolerance not achieved.\n");
            printf("         The infeasibility is %e. SOLNP has restarted %d times.\n", w->bestcon,restart);
        }else if (w->obj_dif > stgs->tol) {
            sol->status = -2; // Fail to converge
            printf("SOLNP+--> Exiting after maximum number of iterations. Tolerance of infeasibility achieved.\n");
            printf("         The infeasibility is %e. SOLNP has restarted %d times.\n", w->bestcon, restart);
            printf("         SOLNP fails to converge.\n");
        }
    }

    sol->iter = i;

    sol->ic = (solnp_float*)solnp_malloc(MAX(nic,1)*sizeof(solnp_float));
    memcpy(sol->ic, w->p, MAX(nic,1)*sizeof(solnp_float));

    sol->jh = (solnp_float*)solnp_malloc((i+1)*sizeof(solnp_float));
    memcpy(sol->jh, w->jh, (i+1)*sizeof(solnp_float));

    sol->ch = (solnp_float*)solnp_malloc((i + 1) * sizeof(solnp_float));
    memcpy(sol->ch, w->ch, (i + 1) * sizeof(solnp_float));
    
    sol->h = (solnp_float*)solnp_malloc((w->nic + w->n) * (w->nic + w->n) * sizeof(solnp_float));
    memcpy(sol->h, w->h, (w->nic + w->n) * (w->nic + w->n)*sizeof(solnp_float));

    sol->obj = w->bestobj;

    sol->count_cost = w->count_cost;

    sol->constraint = w->bestcon;

    return i;

}



solnp_int free_cost(SOLNPCost *cost)
{
    if(cost){
        if(cost->ec){
            solnp_free(cost->ec);
        }
        if(cost->ic){
            solnp_free(cost->ic);
        }
        cost->cost = SOLNP_NULL;
        solnp_free(cost);
    }

    return 0;
}

solnp_int free_constraint(SOLNPConstraint *cnstr)
{
    if(cnstr){
        if(cnstr->il){
            solnp_free(cnstr->il);
        }
        if(cnstr->iu){
            solnp_free(cnstr->iu);
        }
        if(cnstr->pl){
            solnp_free(cnstr->pl);
        }
        if(cnstr->pu){
            solnp_free(cnstr->pu);
        }
        if(cnstr->Ipc){
            solnp_free(cnstr->Ipc);
        }
        if(cnstr->Ipb){
            solnp_free(cnstr->Ipb);
        }
        solnp_free(cnstr);
    }

    return 0;
}


solnp_int free_work(SOLNPWork *w)
{
    if(w){
        if(w->ob){
            free_cost(w->ob);
        }
        if(w->pb){
            free_constraint(w->pb);
        }
        if(w->h){
            solnp_free(w->h);
        }
        if(w->constraint){
            solnp_free(w->constraint);
        }
        if(w->jh){
            solnp_free(w->jh);
        }
        if (w->ch) {
            solnp_free(w->ch);
        }
        if(w->p){
            solnp_free(w->p);
        }
        if(w->l){
            solnp_free(w->l);
        }
        if (w->bestp) {
            solnp_free(w->bestp);
        }
        if (w->bestl) {
            solnp_free(w->bestl);
        }
        solnp_free(w);
    }

    return 0;
}


solnp_int free_input(SOLNPIput *input)
{
    if(input){
        // already freed in SOLNPWork
        //if(input->cnstr){
        //    free_constraint(input->cnstr);
        //} 
        if(input->stgs){
            solnp_free(input->stgs);
        }
        if(input->h){
            //free_hessian(input->h);
            solnp_free(input->h);
        }
        if(input->l){
            solnp_free(input->l);
        }
        solnp_free(input);
    }

    return 0;
}


solnp_int SOLNP(finish)(SOLNPWork *w, SOLNPIput *input)
{
    if(w){
        free_work(w);
    }
    if(input){
        free_input(input);
    }

    return 0;
}




solnp_int SOLNP(main)
(
    SOLNPIput *input,
    SOLNPCost *cost,
    solnp_float *ib0_p,
    SOLNPSol *sol,
    SOLNPInfo *info
)
{
    solnp_int status;
    info->qpsolver_time = 0;



    SOLNPWork *w = SOLNP(init)(input, cost, ib0_p);

    if(w)
    {
        SOLNP(solve)(w, input, ib0_p, sol, info);
        SOLNP(finish)(w, input);
        cost = SOLNP_NULL;
    }

    return 0;
}

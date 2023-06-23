#include "solnp.h"
#include "linalg.h"
#include "subnp.h"
#include "solnp_util.h"
#include "rescue.h"
#include "der_info.h"
SOLNPWork *init_work(
	SOLNPIput *input,
	SOLNPCost *cost)
{
	SOLNPWork *w = (SOLNPWork *)solnp_malloc(sizeof(SOLNPWork));
	/*if (cost->nec >= input->n) {
		input->stgs->rescue = 1;
	}
	else {
		input->stgs->rescue = 0;
	}*/
	if (input->stgs->rs) {
		input->stgs->bfgs = 0;
		input->stgs->scale = 0;
		input->stgs->noise = 0;
		input->stgs->min_iter = 1;
		input->stgs->re_time = input->stgs->max_iter;
	}
	w->n = input->stgs->rescue ? input->n + 2 * cost->nec : input->n;
	w->nec = cost->nec;
	w->nic = cost->nic;
	w->nc = w->nec + w->nic;
	w->ob = cost;
	w->pb = input->cnstr;
	w->rho = input->stgs->rho;
	w->pen_l1 = input->stgs->pen_l1;
	w->restart = 0;

	w->constraint = (solnp_float *)solnp_malloc(w->nc * sizeof(solnp_float));
	w->p = (solnp_float *)solnp_malloc((w->n + w->nic) * sizeof(solnp_float));
	w->p_old = SOLNP_NULL;

	w->l = (solnp_float *)solnp_malloc(MAX(1, w->nc) * sizeof(solnp_float));
	w->best_fea_l = (solnp_float *)solnp_malloc(MAX(1, w->nc) * sizeof(solnp_float));
	w->bestl = (solnp_float *)solnp_calloc(MAX(1, w->nc), sizeof(solnp_float));
	// w->h = (solnp_float*)solnp_calloc((w->nic + w->n) * (w->nic + w->n), sizeof(solnp_float));
	w->jh = (solnp_float *)solnp_malloc((input->stgs->max_iter + 2) * sizeof(solnp_float));
	w->ch = (solnp_float *)solnp_malloc((input->stgs->max_iter + 2) * sizeof(solnp_float));
	w->count_h = (solnp_float *)solnp_malloc((input->stgs->max_iter + 2) * sizeof(solnp_float));
	w->bestp = (solnp_float *)solnp_malloc((input->n) * sizeof(solnp_float));
	w->best_fea_p = (solnp_float *)solnp_malloc((input->n) * sizeof(solnp_float));

	w->best_fea_ob = init_cost(w->nec, w->nic);
	w->bestob = init_cost(w->nec, w->nic);
	w->radius = 1;

	copySOLNPCost(w->best_fea_ob, w->ob);
	w->best_fea_ob->cost = w->ob->cost;

	w->count_cost = 1;
	w->const_time = 0;
	w->count_grad = 0;
	w->count_hess = 0;
	w->exit = 0;

	// Set the Random Seed
	srand((unsigned int)time(NULL));

	return w;
}

solnp_int update_work(
	SOLNPWork *w,
	SOLNPIput *input,
	solnp_float *ib0_p)
{
	solnp_int nec = w->nec;
	solnp_int nic = w->nic;
	solnp_int nc = w->nc;
	solnp_int i, j;

	memcpy(w->p, ib0_p, (input->n + w->nic) * sizeof(solnp_float));
	memcpy(w->best_fea_p, w->p + w->nic, (input->n) * sizeof(solnp_float));
	w->bestobj = INFINITY;
	if (nc > 0.5)
	{
		memcpy(w->l, input->l, nc * sizeof(solnp_float));
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
		w->cons_nm1_orgin = w->cons_nm1;
		if (MAX(w->cons_nm1 - 10 * input->stgs->tol, nic) <= 0)
		{
			w->rho = 0;
		}
	}
	else
	{
		w->l[0] = 0;
		w->constraint = SOLNP_NULL;
	}

	w->mu = w->n;
	w->j = w->ob->obj;
	w->count_h[0] = 1;
	w->jh[0] = w->j;
	w->ch[0] = w->cons_nm1;
	w->best_fea_con = w->cons_nm1;
	// memcpy(w->h, input->h, (w->nic+w->n - 2 * w->nec)* (w->nic + w->n - 2 * w->nec) * sizeof(solnp_float));
	if (input->stgs->bfgs)
	{
		w->h = (solnp_float *)solnp_calloc((w->nic + w->n) * (w->nic + w->n), sizeof(solnp_float));
		for (i = 0; i < input->n + w->nic; i++)
		{
			for (j = 0; j < input->n + w->nic; j++)
			{
				w->h[i + j * w->n] = input->h[i + j * (input->n)];
			}
		}
	}
	else
	{
		w->h = SOLNP_NULL;
	}
	if (input->stgs->rescue == 1)
	{
		// initiate slack variable for equality constraints
		for (i = w->nic + w->n - 2 * w->nec; i < w->nic + w->n; i++)
		{
			w->p[i] = 1;
			if (input->stgs->bfgs)
			{
				w->h[i * (w->nic + w->n) + i] = 1;
			}
		}
		// Modify objective function
		for (i = 0; i < w->nec; i++)
		{
			w->ob->obj += w->pen_l1 * (w->p[i + w->nic + w->n - 2 * w->nec] + w->p[i + w->nic + w->n - w->nec]);
		}

		// Modify bound
		if (w->nec > 0)
		{
			w->pb->n = w->n;
			solnp_float *pl_temp, *pu_temp;
			pl_temp = w->pb->pl;
			pu_temp = w->pb->pu;
			w->pb->pl = (solnp_float *)solnp_malloc(w->n * sizeof(solnp_float));
			w->pb->pu = (solnp_float *)solnp_malloc(w->n * sizeof(solnp_float));
			memcpy(w->pb->pl, pl_temp, (w->n - 2 * w->nec) * sizeof(solnp_float));
			memcpy(w->pb->pu, pu_temp, (w->n - 2 * w->nec) * sizeof(solnp_float));
			w->pb->Ipb[0] = 1;
			w->pb->Ipb[1] = 1;
			for (i = 0; i < 2 * w->nec; i++)
			{
				w->pb->pl[i + w->nic + w->n - 2 * w->nec] = 0;
				w->pb->pu[i + w->nic + w->n - 2 * w->nec] = INFINITY;
			}
			solnp_free(pl_temp);
			solnp_free(pu_temp);
		}
	}
}
/*
SOLNPWork* rescue
(
	SOLNPWork* w,
	SOLNPSettings* stgs,
	SOLNPInfo* info // record running time
) {
	solnp_int ls_time = stgs->ls_time;
	solnp_int max_iter = stgs->max_iter;
	solnp_int min_iter = stgs->min_iter;
	solnp_float tol = stgs->tol;
	stgs->rescue = 1;
	stgs->ls_time = 0;
	stgs->max_iter = stgs->max_iter_rescue;
	stgs->min_iter = stgs->min_iter_rescue;
	if (stgs->max_iter == 0) {
		return w;
	}
	stgs->tol = MAX(stgs->tol* stgs->tol*10,1e-6);
	solnp_printf("SOLNP+--> Rescue Process begins. Finding a feasible solution...\n");
	SOLNPWork* w_rescue = SOLNP_RESCUE(w, stgs, info);

	w_rescue->n -= w_rescue->nec;

	SOLNPCost** ob = &w_rescue->ob;
	w->ob->cost(ob, w_rescue->bestp + w->nic, w->n, 1, 0);
	memcpy(w_rescue->best_fea_p, w_rescue->bestp + w->nic, w->n * sizeof(float));
	memcpy(w_rescue->best_fea_l, w_rescue->best_fea_l + w->nic, w->n * sizeof(float));
	copySOLNPCost(w_rescue->best_fea_ob, w_rescue->ob);

	memcpy(w_rescue->constraint, w_rescue->ob->ec, w_rescue->nec * sizeof(solnp_float));
	memcpy(&w_rescue->constraint[w->nec], w_rescue->ob->ic, w_rescue->nic * sizeof(solnp_float));

	if (w->nic > 0.5) {
		SOLNP(add_scaled_array)(&w_rescue->constraint[w->nec], w->p, w_rescue->nic, -1.0);
	}

	w_rescue->cons_nm2 = SOLNP(norm)(w_rescue->constraint, w_rescue->nc);
	w_rescue->j = w_rescue->ob->obj;
	w_rescue->cons_nm1 = w_rescue->cons_nm2;

	w_rescue->bestcon = w_rescue->cons_nm2;
	w_rescue->bestobj = w_rescue->ob->obj;
	stgs->rescue = 0;
	stgs->ls_time = ls_time;
	stgs->tol = tol;
	stgs->max_iter = max_iter;
	stgs->min_iter = min_iter;
	free_work(w);
	return w_rescue;
}
*/
SOLNPWork *SOLNP(init)(
	SOLNPIput *input,
	SOLNPCost *cost,
	solnp_float *ib0_p)
{
	SOLNPWork *w;
	w = init_work(input, cost);

	return w;
}

void restart(
	SOLNPWork *w,
	SOLNPSettings *stgs,
	solnp_int iter,
	solnp_float delta0)
{
	solnp_int n = stgs->rescue ? w->n - 2 * w->nec : w->n;
	w->restart++;
	if (w->exit != 0 || w->restart > 2)
	{
		if (w->bestobj != INFINITY)
		{
			memcpy(w->p, w->bestp, n * sizeof(solnp_float));
			memcpy(w->l, w->bestl, MAX(w->nc, 1) * sizeof(solnp_float));
			copySOLNPCost(w->ob, w->bestob);
			w->cons_nm1_orgin = w->bestcon;
		}
		else
		{
			memcpy(w->p, w->best_fea_p, n * sizeof(solnp_float));
			memcpy(w->l, w->best_fea_l, MAX(w->nc, 1) * sizeof(solnp_float));
			copySOLNPCost(w->ob, w->best_fea_ob);
			w->cons_nm1_orgin = w->best_fea_con;
			w->bestcon = w->best_fea_con;
		}
		w->jh[iter] = w->ob->obj;
		w->ch[iter] = w->cons_nm1_orgin;
		w->pen_l1 = 1;
		w->j = w->ob->obj + w->pen_l1 * 2 * w->nec;
		if (stgs->bfgs)
		{
			for (int diag = 0; diag < w->nic + w->n; diag++)
			{
				w->h[diag * (w->nic + w->n) + diag] = 1;
			}
		}
		if (stgs->rescue)
		{
			for (int k = 0; k < 2 * w->nec; k++)
			{
				w->p[w->nic + w->n - w->nec * 2 + k] = 1;
			}
		}
		w->exit = 0;
		w->mu = w->n;
	}
	if (stgs->noise)
	{
		stgs->delta = delta0;
	}
	if (stgs->drsom)
	{
		// solnp_free(w->p_old); // will be freed in free_work
		w->radius = 1;
	}
	// h=diag(diag(h));
	if (stgs->bfgs)
	{
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
}

SOLNPWork *SOLNP(solve)(
	SOLNPWork *w,
	SOLNPIput *input,
	solnp_float *ib0_p,
	SOLNPSol *sol,
	SOLNPInfo *info)
{
	solnp_int i = 0;
	solnp_int j;
	SOLNPSettings *stgs = input->stgs;
	solnp_int n = stgs->rescue ? w->n - 2 * w->nec : w->n;
	solnp_int nc = w->nc;
	solnp_int nec = w->nec;
	solnp_int nic = w->nic;
	solnp_int start_chg_l1_pen = w->nec > 3 * n ? 12 : 5;

	// solnp_int rescue_tag = 0;
	solnp_int max_iter = stgs->max_iter;
	solnp_int iter = 0;

	solnp_float delta0 = stgs->delta;
	update_work(w, input, ib0_p);

	for (i = 0; i < max_iter; i++)
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
			w->cons_nm2_orgin = stgs->rescue ? calculate_infeas_scaledob_l1(w->ob, w, w->p) : w->cons_nm2;
			// previous 10
			if (w->cons_nm2 < 1 * stgs->tol_con)
			{
				w->rho = 0;
				w->mu = MIN(w->mu, stgs->tol);
			}
			// previously 5
			// adjust the ALM penalty parameter
			if (w->cons_nm2 < w->cons_nm1 && w->cons_nm2 < 10 * stgs->tol_con)
			{
				w->rho /= 5;
			}
			else if (w->cons_nm2 > 10 * w->cons_nm1)
			{
				w->rho = 5 * MAX(w->rho, SQRTF(stgs->tol));
			}

			// adjust the l1 penalty parameter
			if (stgs->rescue && w->cons_nm2_orgin <= 3 * w->cons_nm1_orgin && w->cons_nm2_orgin < 10 * stgs->tol_con)
			{
				w->pen_l1 = w->pen_l1 / 2;
			}
			if (stgs->rescue && w->cons_nm2_orgin >= 5 * w->cons_nm1_orgin || (w->cons_nm2_orgin > 100 * stgs->tol_con && i >= MIN(stgs->max_iter / 4, start_chg_l1_pen)))
			{
				w->pen_l1 = MIN(w->pen_l1 * 3, 1e4);
			}
			if (w->exit == 0 && MAX(w->obj_dif, 3 * w->cons_nm1 - w->cons_nm2) <= 0 && w->restart < stgs->re_time && (stgs->delta <= MAX(3 * stgs->tol, stgs->delta_end) || stgs->grad))
			{
				SOLNP(scale_array)
				(w->l, 0, nc);
				w->restart++;
				// h=diag(diag(h));
				if (stgs->noise)
				{
					stgs->delta = delta0;
				}
				if (stgs->bfgs)
				{
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
			}

			w->cons_nm1 = w->cons_nm2;
			w->cons_nm1_orgin = w->cons_nm2_orgin;
		}

		w->count_h[1 + iter] = w->count_cost;
		w->jh[1 + iter] = w->j;
		w->ch[1 + iter] = w->cons_nm1_orgin;
		iter++;
		/*
		if (i >= max_iter / 2 && rescue_tag == 0 && w->cons_nm1 > 20 * stgs->tol_con) {

			rescue_tag = 1;
			stgs->delta = delta0;

			solnp_float* hisj = (solnp_float*)solnp_malloc(2 * (stgs->max_iter + 1) * sizeof(solnp_float));
			solnp_float* hisc = (solnp_float*)solnp_malloc(2 * (stgs->max_iter + 1) * sizeof(solnp_float));
			memcpy(hisj, w->jh, (iter + 1) * sizeof(solnp_float));
			memcpy(hisc, w->ch, (iter + 1) * sizeof(solnp_float));
			w = rescue(w, stgs, info);

			j = 0;
			while (w->jh[j] != INFINITY) {
				j++;
			}
			memcpy(hisj + (iter + 1), w->jh, j * sizeof(solnp_float));
			memcpy(hisc + (iter + 1), w->ch, j * sizeof(solnp_float));
			solnp_free(w->jh);
			solnp_free(w->ch);
			w->jh = hisj;
			w->ch = hisc;
			iter += j;
			continue;
		}
		*/
		// if (stgs->rs) {
		//	//FOR TEST RANDOM SAMPLING ONLY
		//	w->alm_crit = 1;
		// }
		if (w->nc && stgs->rescue &&
			w->const_time > 5 * (w->n - 2 * w->nec) && w->best_fea_con <= stgs->tol_con)
		{
			break;
		}
		if ((w->cons_nm1_orgin <= stgs->tol_con && w->obj_dif <= stgs->tol) || (stgs->drsom && w->radius <= stgs->tol / 10) || w->exit)
		{
			if (w->restart < stgs->re_time && (w->alm_crit > stgs->tol_restart || (w->exit != 0) || (stgs->drsom && w->radius > stgs->tol / 10)) && w->exit != 1)
			{
				// Restart from the best point or feasible point.
				restart(w, stgs, iter, delta0);
			}
			else
			{
				i++;
				break;
			}
		}
	}
	if (w->bestobj != INFINITY)
	{
		memcpy(w->p, w->bestp, n * sizeof(solnp_float));
		memcpy(w->l, w->bestl, MAX(nc, 1) * sizeof(solnp_float));
		copySOLNPCost(w->ob, w->bestob);
		w->cons_nm1_orgin = w->bestcon;
	}
	else
	{
		memcpy(w->p, w->best_fea_p, n * sizeof(solnp_float));
		memcpy(w->l, w->best_fea_l, MAX(nc, 1) * sizeof(solnp_float));
		copySOLNPCost(w->ob, w->best_fea_ob);
		w->cons_nm1_orgin = w->best_fea_con;
		w->bestcon = w->best_fea_con;
	}
	/*
	if (w->nec && ( isnan(w->cons_nm1) || w->cons_nm1 > 10 * stgs->tol_con) && stgs->max_iter_rescue) {
		stgs->delta = delta0;

		solnp_float* hisj = (solnp_float*)solnp_malloc(2 * (stgs->max_iter+1) * sizeof(solnp_float));
		solnp_float* hisc = (solnp_float*)solnp_malloc(2 * (stgs->max_iter + 1) * sizeof(solnp_float));
		memcpy(hisj, w->jh, (iter+1) * sizeof(solnp_float));
		memcpy(hisc, w->ch, (iter+1) * sizeof(solnp_float));
		w = rescue(w, stgs, info);
		j = 0;
		while (w->jh[j] != INFINITY) {
			j++;
		}
		memcpy(hisj+iter+1, w->jh, j * sizeof(solnp_float));
		memcpy(hisc+iter+1, w->ch, j * sizeof(solnp_float));
		solnp_free(w->jh);
		solnp_free(w->ch);
		w->jh = hisj;
		w->ch = hisc;
		iter += j;
	}*/
	sol->p = (solnp_float *)solnp_malloc(n * sizeof(solnp_float));
	sol->l = (solnp_float *)solnp_malloc(MAX(1, nc) * sizeof(solnp_float));
	if (w->bestobj != INFINITY)
	{
		memcpy(sol->p, w->bestp, n * sizeof(solnp_float));
		memcpy(sol->l, w->bestl, MAX(nc, 1) * sizeof(solnp_float));
	}
	else
	{
		memcpy(sol->p, &w->p[w->nic], n * sizeof(solnp_float));
		memcpy(sol->l, w->l, MAX(nc, 1) * sizeof(solnp_float));
		w->bestcon = w->best_fea_con;
		w->bestobj = w->ob->obj;
	}

	if (w->cons_nm1_orgin <= stgs->tol_con && ((w->obj_dif <= stgs->tol && stgs->delta <= MAX(stgs->tol, stgs->delta_end)) || w->const_time > 5 * (w->n - 2 * w->nec)) || stgs->grad == 1)
	{
		sol->status = 1; // Success
		printf("SOLNP+--> Success! Completed in %d iterations\n", iter);
		printf("         The infeasibility is %e.\n", w->bestcon);
	}
	else
	{
		if (w->exit == 1)
		{
			sol->status = 0;
			printf("SOLNP+--> Exiting after maximum number of function evaluation. Tolerance not achieved.\n");
			printf("         The infeasibility is %e.\n", w->bestcon);
			printf("         SOLNP has restarted %d times.\n", w->restart);
		}
		else if (w->exit == 2)
		{
			sol->status = -3;
			printf("SOLNP+--> Exiting because of unknown error. Tolerance not achieved.\n");
			printf("         The infeasibility is %e.\n", w->bestcon);
		}
		else if (w->cons_nm1_orgin > stgs->tol_con)
		{
			sol->status = -1; // Fail to find a feasible point.
			printf("SOLNP+--> Exiting after maximum number of iterations. Tolerance not achieved.\n");
			printf("         The infeasibility is %e. SOLNP has restarted %d times.\n", w->bestcon, w->restart);
		}
		else if (w->obj_dif > stgs->tol)
		{
			sol->status = -2; // Fail to converge
			printf("SOLNP+--> Exiting after maximum number of iterations. Tolerance of infeasibility achieved.\n");
			printf("         The infeasibility is %e. SOLNP has restarted %d times.\n", w->bestcon, w->restart);
			printf("         SOLNP fails to converge.\n");
		}
	}

	sol->iter = iter;

	sol->ic = (solnp_float *)solnp_malloc(MAX(nic, 1) * sizeof(solnp_float));
	if (nic == 0)
	{
		*sol->ic = 0;
	}
	else
	{
		memcpy(sol->ic, w->p, MAX(nic, 1) * sizeof(solnp_float));
	}

	sol->best_fea_p = (solnp_float *)solnp_malloc((n) * sizeof(solnp_float));
	memcpy(sol->best_fea_p, w->best_fea_p, (n) * sizeof(solnp_float));

	sol->count_h = (solnp_float *)solnp_malloc((iter + 1) * sizeof(solnp_float));
	memcpy(sol->count_h, w->count_h, (iter + 1) * sizeof(solnp_float));

	sol->jh = (solnp_float *)solnp_malloc((iter + 1) * sizeof(solnp_float));
	memcpy(sol->jh, w->jh, (iter + 1) * sizeof(solnp_float));

	sol->ch = (solnp_float *)solnp_malloc((iter + 1) * sizeof(solnp_float));
	memcpy(sol->ch, w->ch, (iter + 1) * sizeof(solnp_float));

	if (stgs->bfgs)
	{
		sol->h = (solnp_float *)solnp_malloc((nic + n) * (nic + n) * sizeof(solnp_float));
		memcpy(sol->h, w->h, (nic + n) * (nic + n) * sizeof(solnp_float));
	}
	else
	{
		sol->h = (solnp_float *)solnp_malloc(1 * sizeof(solnp_float));
		*sol->h = 1;
	}

	sol->obj = w->bestobj;

	sol->count_cost = w->count_cost;

	sol->count_grad = w->count_grad;

	sol->count_hess = w->count_hess;

	sol->constraint = w->bestcon;

	sol->restart_time = w->restart;

	return w;
}

solnp_int free_cost(SOLNPCost *cost)
{
	if (cost)
	// if (0)
	{
		if (cost->ec)
		{
			solnp_free(cost->ec);
		}
		if (cost->ic)
		{
			solnp_free(cost->ic);
		}
		if (cost->cost)
		{
			cost->cost = SOLNP_NULL;
		}
		if (cost->grad)
		{
			cost->grad = SOLNP_NULL;
		}
		if (cost->hess)
		{
			cost->hess = SOLNP_NULL;
		}
		solnp_free(cost);
	}

	return 0;
}

solnp_int free_constraint(SOLNPConstraint *cnstr)
{
	if (cnstr)
	{
		if (cnstr->il)
		{
			solnp_free(cnstr->il);
		}
		if (cnstr->iu)
		{
			solnp_free(cnstr->iu);
		}
		if (cnstr->pl)
		{
			solnp_free(cnstr->pl);
		}
		if (cnstr->pu)
		{
			solnp_free(cnstr->pu);
		}
		if (cnstr->Ipc)
		{
			solnp_free(cnstr->Ipc);
		}
		if (cnstr->Ipb)
		{
			solnp_free(cnstr->Ipb);
		}
		solnp_free(cnstr);
	}

	return 0;
}

solnp_int free_work(SOLNPWork *w)
{
	if (w)
	{
		if (w->ob)
		{
			free_cost(w->ob);
		}
		if (w->best_fea_ob)
		{
			free_cost(w->best_fea_ob);
		}
		if (w->bestob)
		{
			free_cost(w->bestob);
		}
		if (w->pb)
		{
			free_constraint(w->pb);
		}
		if (w->h)
		{
			solnp_free(w->h);
		}
		if (w->constraint)
		{
			solnp_free(w->constraint);
		}
		if (w->count_h)
		{
			solnp_free(w->count_h);
		}
		if (w->jh)
		{
			solnp_free(w->jh);
		}
		if (w->ch)
		{
			solnp_free(w->ch);
		}
		if (w->p)
		{
			solnp_free(w->p);
		}
		if (w->p_old)
		{
			solnp_free(w->p_old);
		}
		if (w->best_fea_p)
		{
			solnp_free(w->best_fea_p);
		}
		if (w->l)
		{
			solnp_free(w->l);
		}
		if (w->best_fea_l)
		{
			solnp_free(w->best_fea_l);
		}
		if (w->bestp)
		{
			solnp_free(w->bestp);
		}
		if (w->bestl)
		{
			solnp_free(w->bestl);
		}
		solnp_free(w);
	}

	return 0;
}

solnp_int free_input(SOLNPIput *input)
{
	if (input)
	{
		// already freed in SOLNPWork
		// if(input->cnstr){
		//    free_constraint(input->cnstr);
		//}
		if (input->stgs)
		{
			solnp_free(input->stgs);
		}
		if (input->h)
		{
			// free_hessian(input->h);
			solnp_free(input->h);
		}
		if (input->l)
		{
			solnp_free(input->l);
		}
		solnp_free(input);
	}

	return 0;
}

solnp_int SOLNP(finish)(SOLNPWork *w, SOLNPIput *input)
{
	if (w)
	{
		free_work(w);
	}
	if (input)
	{
		free_input(input);
	}

	return 0;
}

solnp_int SOLNP(main)(
	SOLNPIput *input,
	SOLNPCost *cost,
	solnp_float *ib0_p,
	SOLNPSol *sol,
	SOLNPInfo *info)
{
	solnp_int status;
	info->qpsolver_time = 0;

	SOLNPWork *w = SOLNP(init)(input, cost, ib0_p);

	if (w)
	{
		w = SOLNP(solve)(w, input, ib0_p, sol, info);
		SOLNP(finish)
		(w, input);
		cost = SOLNP_NULL;
	}

	return 0;
}

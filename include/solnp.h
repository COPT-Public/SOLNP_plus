#pragma once
#ifndef SOLNP_H_GUARD
#define SOLNP_H_GUARD

#ifdef __cplusplus
extern "C"
{
#endif

#include "solnp_glbopts.h"
#include <string.h>

      typedef struct SOLNP_WORK SOLNPWork;
      typedef struct SOLNP_COST SOLNPCost;
      typedef struct SOLNP_SETTINGS SOLNPSettings;
      typedef struct SOLNP_CONSTRAINT SOLNPConstraint;
      typedef struct SOLNP_INPUT SOLNPIput;
      typedef struct SOLNP_SOL SOLNPSol;
      typedef struct SOLNP_INFO SOLNPInfo;

      struct SOLNP_INFO
      {
            solnp_float total_time;
            solnp_float qpsolver_time;
            solnp_float cost_time;
      };

      struct SOLNP_INPUT
      {
            solnp_int n;
            SOLNPConstraint *cnstr;
            SOLNPSettings *stgs;
            solnp_float *l;
            solnp_float *h;
      };

      // define the lower and upper bounds of decision variable and inequality constrains
      struct SOLNP_CONSTRAINT
      {
            solnp_int n;      // dimension of decision variable
            solnp_int nic;    // number of inequality constrains
            solnp_float *il;  // lower bounds of inequality constrains, len of nic
            solnp_float *iu;  // upper bounds of inequality constrains, len of nic
            solnp_float *pl;  // lower bounds of decision variable, len of n
            solnp_float *pu;  // upper bounds of decision variable, len of n
            solnp_float *Ipc; // len of 2
            solnp_float *Ipb; // len of 2
      };

      // return of cost function
      struct SOLNP_COST
      {
            solnp_int nec;   // number of equality constrains
            solnp_int nic;   // number of inequality constrains
            solnp_float obj; // objective of cost function
            solnp_float *ec; // constraint function values of EC, len of nec
            solnp_float *ic; // constraint function values of IC, len of nic

            void (*hess)(solnp_float *h, solnp_float *p, solnp_int np, solnp_int nheval, solnp_int action); // function pointer to user-defined hess function
            void (*grad)(solnp_float *g, solnp_float *p, solnp_int np, solnp_int ngeval, solnp_int action); // function pointer to user-defined gradient function
            void (*cost)(SOLNPCost **c, solnp_float *p, solnp_int np, solnp_int nfeval, solnp_int action);  // function pointer to user-defined cost function
      };
      struct SOLNP_WORK
      {
            solnp_int n;   // dimension of decision variable
            solnp_int nec; // number of equality constrains
            solnp_int nic; // number of inequality constrains
            solnp_int nc;  // number constrains

            SOLNPCost *ob;           // observation of cost function
            solnp_float *constraint; // vector of constraint values
            solnp_float obj_dif;     // the difference of the objective values between two consecutive iterations
            solnp_float cons_nm1;    // NORM(CONSTRAINT) before a major iteration
            solnp_float cons_nm2;    // NORM(CONSTRAINT) after a major iteration

            solnp_float cons_nm1_orgin; // NORM(CONSTRAINT) before a major iteration
            solnp_float cons_nm2_orgin; // NORM(CONSTRAINT) after a major iteration

            SOLNPConstraint *pb;

            solnp_float j;
            solnp_float *count_h;
            solnp_float *jh;
            solnp_float *ch;
            solnp_float rho;
            solnp_float pen_l1; // l1 penalty paramenter
            solnp_float radius;

            solnp_float mu;
            solnp_float *p;
            solnp_float *p_old;
            // solnp_float *ib0;
            solnp_float *l;
            solnp_float *h; // Hessian matrix
            solnp_float *bestp;
            solnp_float *bestl;
            SOLNPCost *bestob;

            solnp_float bestcon;
            solnp_float bestobj;
            solnp_float alm_crit; // use alm stop criterion to decide whether to restart

            solnp_float *best_fea_p;
            solnp_float best_fea_con;
            solnp_float *best_fea_l;
            SOLNPCost *best_fea_ob;

            solnp_int const_time;
            solnp_int restart;

            solnp_int count_cost; // record cost times in solnp.c
            solnp_int count_grad; // record gradient times in solnp.c
            solnp_int count_hess;
            solnp_int exit;
      };

      struct SOLNP_SETTINGS
      {
            solnp_float pen_l1;
            solnp_float rho;
            solnp_int max_iter;
            solnp_int min_iter;
            solnp_int max_iter_rescue;
            solnp_int min_iter_rescue;
            solnp_float delta;
            solnp_float tol;
            solnp_float tol_con;
            solnp_int ls_time;
            solnp_float tol_restart;
            solnp_int re_time;
            solnp_float delta_end;
            solnp_int maxfev;
            solnp_int noise;
            solnp_int qpsolver;
            solnp_float k_r;
            solnp_float k_i;
            solnp_float c_r;
            solnp_float c_i;
            solnp_int batchsize;
            solnp_int hess;
            solnp_int grad;
            solnp_int rescue;
            solnp_int ls_way; // 1 means bisection, 2 means non-monotonic
            solnp_int bfgs;
            solnp_int rs;
            solnp_int cen_diff;
            solnp_float gd_step;
            solnp_int scale;
            solnp_int drsom;
      };

      struct SOLNP_SOL
      {
            solnp_int iter;
            solnp_int status;
            solnp_float *p;
            solnp_float *ic;
            solnp_float *count_h;
            solnp_float *jh;
            solnp_float *ch;
            solnp_float *best_fea_p;
            solnp_float *l;
            solnp_float *h;
            solnp_int count_cost;
            solnp_int count_grad; // record gradient times in solnp.c
            solnp_int count_hess;
            solnp_float constraint;
            solnp_float obj;
            solnp_int restart_time;
      };

      solnp_int SOLNP(main)(
          SOLNPIput *input,
          SOLNPCost *cost,
          solnp_float *ib0_p,
          SOLNPSol *sol,
          SOLNPInfo *info);

      solnp_int free_cost(SOLNPCost *cost);

      solnp_int subnp_qp(SOLNPWork *w, SOLNPSettings *stgs, SOLNPInfo *info);

#ifdef __cplusplus
}
#endif
#endif

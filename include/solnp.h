#ifndef SOLNP_H_GUARD
#define SOLNP_H_GUARD

#ifdef __cplusplus
extern "C" {
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
      solnp_float* h;
};


// define the lower and upper bounds of decision variable and inequality constrains
struct SOLNP_CONSTRAINT
{
      solnp_int n;   // dimension of decision variable
      solnp_int nic; // number of inequality constrains
      solnp_float *il; // lower bounds of inequality constrains, len of nic
      solnp_float *iu; // upper bounds of inequality constrains, len of nic
      solnp_float *pl; // lower bounds of decision variable, len of n
      solnp_float *pu; // upper bounds of decision variable, len of n
      solnp_float *Ipc; // len of 2
      solnp_float *Ipb; // len of 2

};

//return of cost function
struct SOLNP_COST
{     
      solnp_int nec; // number of equality constrains
      solnp_int nic; // number of inequality constrains
      solnp_float obj; // objective of cost function
      solnp_float *ec; // constraint function values of EC, len of nec
      solnp_float *ic; // constraint function values of IC, len of nic

      void (*cost)(SOLNPCost *c, solnp_float *p, solnp_int np); // function pointer to user-defined cost function
};
struct SOLNP_WORK 
{
      solnp_int n;   // dimension of decision variable
      solnp_int nec; // number of equality constrains
      solnp_int nic; // number of inequality constrains
      solnp_int nc; // number constrains

      SOLNPCost *ob; // observation of cost function
      solnp_float *constraint; // vector of constraint values
      solnp_float obj_dif; // the difference of the objective values between two consecutive iterations
      solnp_float cons_nm1; // NORM(CONSTRAINT) before a major iteration
      solnp_float cons_nm2; // NORM(CONSTRAINT) after a major iteration    
      SOLNPConstraint *pb;         

      solnp_float j;
      solnp_float *jh;
      solnp_float* ch;
      solnp_float rho;

      solnp_float mu;
      solnp_float *p;
      // solnp_float *ib0;
      solnp_float *l;
      solnp_float* h; // Hessian matrix
      solnp_float* bestp;
      solnp_float* bestl;
      solnp_float bestcon;
      solnp_float bestobj;
      solnp_float alm_crit; // use alm stop criterion to decide whether to restart

      solnp_int count_cost; // record cost times in solnp.c

      solnp_int exit;
};

struct SOLNP_SETTINGS
{
      solnp_float rho;
      solnp_int max_iter;
      solnp_int min_iter;
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
};

struct SOLNP_SOL
{
      solnp_int iter;
      solnp_int status;
      solnp_float *p;
      solnp_float *ic;
      solnp_float *jh;
      solnp_float* ch;
      solnp_float *l;
      solnp_float *h;
      solnp_int count_cost;
      solnp_float constraint;
      solnp_float obj;
};


solnp_int SOLNP(main)
(
    SOLNPIput *input,
    SOLNPCost *cost,
    solnp_float *ib0_p,
    SOLNPSol *sol,
    SOLNPInfo *info
);


solnp_int free_cost(SOLNPCost *cost);



#ifdef __cplusplus
}
#endif
#endif

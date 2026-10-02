#pragma once
#ifndef SUBNP_H_GUARD
#define SUBNP_H_GUARD
#ifdef __cplusplus
extern "C" {
#endif

#include "solnp.h"
typedef struct SUBNP_WORK SUBNPWork;
typedef struct SUBNP_JACOB SUBNPJacob;

struct SUBNP_JACOB
{
    solnp_int n;
    solnp_int nec;
    solnp_int nic;
    solnp_int nc;
    solnp_int npic;
    solnp_int a_row;
    solnp_float *g; // gradient of f(s, x)
    solnp_float *a; // Jacob of g(s, x) = 0
    solnp_float* anew;
};


struct SUBNP_WORK
{
    solnp_float *alp;
    solnp_int n_scale; // length of scale
    solnp_float *scale;
    solnp_float* p_cand;
    solnp_int ch;
    solnp_int mm;
    solnp_int nc;
    solnp_int npic;
    SOLNPCost* ob_cand;
    SUBNPJacob *J;
    solnp_float *b; // rhs if linear constraints exist
    solnp_int* constr_index;
};
solnp_int unscale_ob
(
    SOLNPCost* ob,
    const solnp_float* scale
);
solnp_int rescale_ob
(
    SOLNPCost* ob,
    const solnp_float* scale
);
solnp_int subnp_qp(SOLNPWork* w, SOLNPSettings* stgs, SOLNPInfo* info);
SOLNPCost* init_cost(solnp_int nec, solnp_int nic);
void copySOLNPCost(SOLNPCost* ob1, SOLNPCost* ob2);
// solnp_float qpsolver_time = 0;

#ifdef __cplusplus
}
#endif
#endif
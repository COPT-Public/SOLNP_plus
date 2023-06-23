#pragma once
#include "subnp.h"
#include "linalg.h"
#include "linsys.h"
#include <mkl.h>
#include "osqp.h"
#include <math.h>
#include "solnp_util.h" 

void proj
(
    solnp_float* p,
    SOLNPWork* w
);
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
);
solnp_int find_int_feas_sol_osqp
(
    SUBNPWork* w_sub,
    SOLNPWork* w,
    SOLNPSettings* stgs
);
solnp_int find_int_feas_sol_aff
(
    SUBNPWork* w_sub,
    SOLNPWork* w,
    SOLNPSettings* stgs
);
void solve_qp_osqp
(
    SOLNPWork* w,
    SUBNPWork* w_sub,
    solnp_float* p0,
    solnp_float* y,
    solnp_float* g
);
void solve_qp_aff
(
    SOLNPWork* w,
    SUBNPWork* w_sub,
    solnp_float* p0,
    solnp_float* y,
    solnp_float* g,
    solnp_int sys
);
void Gradient_descent(
    SOLNPWork* w,
    solnp_float* p0,
    solnp_float* g,
    solnp_float stepsize
);
solnp_float* sub_trustreg
(
    solnp_float* Q,
    solnp_float* c,
    solnp_float* G,
    solnp_float delta,
    solnp_float tol
);

solnp_float drsom(
    SOLNPWork* w,
    SUBNPWork* w_sub,
    SOLNPSettings* stgs,
    solnp_float* p0,
    solnp_float* g,
    solnp_float* m,
    solnp_float radius,
    solnp_float val
);
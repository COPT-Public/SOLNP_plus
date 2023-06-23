#include "subnp.h"

solnp_float calculate_infeas_scaledob
(
    SOLNPCost* ob,
    SOLNPWork* w,
    solnp_float* x,
    solnp_float* scale
);
solnp_float calculate_infeas_unscaleob
(
    SOLNPCost* ob,
    SOLNPWork* w,
    solnp_float* x,
    solnp_float* scale
);
solnp_float calculate_infeas_scaledob_l1
(
    SOLNPCost* ob,
    SOLNPWork* w,
    solnp_float* x
);
solnp_float calculate_almcrit_iq
(
    SOLNPWork* w,
    solnp_float* scale,
    solnp_float* p
);
void calculate_scaled_cost
(
    SOLNPCost** ob,
    solnp_float* p,
    const solnp_float* scale,
    SOLNPSettings* stgs,
    SOLNPWork* w,
    solnp_int nfeval
);
solnp_float line_search_merit(
    SOLNPCost* ob,
    SOLNPWork* w,
    SUBNPWork* w_sub,
    solnp_float* p
);
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
);*/
solnp_int calculate_scaled_grad
(
    solnp_float* g,
    solnp_float* p,
    const solnp_float* scale,
    SOLNPSettings* stgs,
    SOLNPWork* w
);
solnp_int calculate_scaled_grad_random
(
    solnp_float* g,
    solnp_float* p,
    SOLNPCost* ob_p,
    const solnp_float* scale,
    SOLNPSettings* stgs,
    SOLNPWork* w
);
solnp_int calculate_scaled_hess
(
    solnp_float* h,
    solnp_float* p,
    const solnp_float* scale,
    SOLNPSettings* stgs,
    SOLNPWork* w
);
void calculate_alm_criterion
(
    SOLNPWork* w,
    SUBNPWork* w_sub,
    solnp_float* grad
);
solnp_int calculate_Jacob_zero
(
    SUBNPWork* w_sub,
    SOLNPWork* w,
    SOLNPSettings* stgs
);
solnp_int calculate_Jacob_zero_rescue
(
    SUBNPWork* w_sub,
    SOLNPWork* w,
    SOLNPSettings* stgs
);
solnp_int calculate_ALMgradient_zero_rescue
(
    SUBNPWork* w_sub,
    SOLNPWork* w,
    SOLNPSettings* stgs,
    solnp_float* g,
    solnp_float j
);
solnp_int calculate_Jacob_first
(
    SUBNPWork* w_sub,
    SOLNPWork* w,
    SOLNPSettings* stgs
);
solnp_int calculate_Jacob_first_rescue
(
    SUBNPWork* w_sub,
    SOLNPWork* w,
    SOLNPSettings* stgs
);
solnp_float calculate_ALM
(
    SOLNPCost* ob,
    SOLNPSettings* stgs,
    solnp_float* p,
    const SOLNPWork* w,
    const SUBNPWork* w_sub
);
//solnp_float calculate_ALM_rescue
//(
//    SOLNPCost* ob,
//    SOLNPSettings* stgs,
//    solnp_float* p,
//    const SOLNPWork* w,
//    const SUBNPWork* w_sub
//);
solnp_int calculate_ALMgradient_zero
(
    SUBNPWork* w_sub,
    SOLNPWork* w,
    SOLNPSettings* stgs,
    solnp_float* g,
    solnp_float j
);
solnp_int calculate_ALMgradient_first
(
    SUBNPWork* w_sub,
    SOLNPWork* w,
    SOLNPSettings* stgs,
    solnp_float* g
);
solnp_int calculate_ALMgradient_first_rescue
(
    SUBNPWork* w_sub,
    SOLNPWork* w,
    SOLNPSettings* stgs,
    solnp_float* g,
    solnp_float j
);
solnp_int calculate_ALM_hess
(
    SUBNPWork* w_sub,
    SOLNPWork* w,
    SOLNPSettings* stgs,
    solnp_float* h
);
void BFGSudpate
(
    SOLNPWork* w,
    SUBNPWork* w_sub,
    SOLNPSettings* stgs,
    solnp_float* g,
    solnp_float* yg,
    solnp_float* sx
);
solnp_float fun_along_2d
(
    SOLNPWork* w,
    SOLNPSettings* stgs,
    SUBNPWork* w_sub,
    solnp_float* d1,
    solnp_float* d2,
    solnp_float* x,
    solnp_float* coeff
);
solnp_float* interpolate2d
(
    SOLNPWork* w,
    SOLNPSettings* stgs,
    SUBNPWork* w_sub,
    solnp_float radius,
    solnp_float* d1,
    solnp_float* d2,
    solnp_float* x,
    solnp_float val
);
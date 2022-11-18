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

    solnp_float *g; // gradient of f(s, x)
    solnp_float *a; // Jacob of g(s, x) = 0
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
};

solnp_int subnp(SOLNPWork *w, SOLNPSettings *stgs);

solnp_int subnp_qp(SOLNPWork *w, SOLNPSettings *stgs, SOLNPInfo *info);

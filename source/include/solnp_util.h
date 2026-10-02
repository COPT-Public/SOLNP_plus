#ifndef UTIL_H_GUARD
#define UTIL_H_GUARD

#ifdef __cplusplus
extern "C"
{
#endif

#include "solnp.h"
#include <stdlib.h>
#include <stdio.h>

#if (defined NOTIMER)
      typedef void *SOLNP(timer);

#elif (defined _WIN32 || defined _WIN64 || defined _WINDLL)

#include <windows.h>
typedef struct SOLNP(timer)
{
      LARGE_INTEGER tic;
      LARGE_INTEGER toc;
      LARGE_INTEGER freq;
} SOLNP(timer);

#elif (defined __APPLE__)

#include <mach/mach_time.h>
typedef struct SOLNP(timer)
{
      uint64_t tic;
      uint64_t toc;
      mach_timebase_info_data_t tinfo;
} SOLNP(timer);

#else

#include <time.h>
typedef struct SOLNP(timer)
{
      struct timespec tic;
      struct timespec toc;
} SOLNP(timer);

#endif

#if EXTRA_VERBOSE > 1
      extern SOLNP(timer) global_timer;
#endif

      void SOLNP(tic)(SOLNP(timer) * t);
      solnp_float SOLNP(toc)(SOLNP(timer) * t);
      solnp_float SOLNP(str_toc)(char *str, SOLNP(timer) * t);
      solnp_float SOLNP(tocq)(SOLNP(timer) * t);

      SOLNPConstraint *malloc_constriant(solnp_int n, solnp_int nic);
      SOLNPCost *malloc_cost(solnp_int nec, solnp_int nic, void *func, void *grad, void *hess);
      void SOLNP(set_default_settings)(SOLNPSettings *stgs, solnp_int np);
      void SOLNP(free_sol)(SOLNPSol *sol);
      void SOLNP(free_info)(SOLNPInfo *info);
      void SOLNP(free_stgs)(SOLNPSettings *stgs);
      SOLNPProb *SOLNP(init_prob)(solnp_int np, solnp_int nic, solnp_int nec);
      void SOLNP(free_prob)(SOLNPProb *prob);

      // C interface for solnp
      // typedef void cost_temple(double *p, double *result, int np, int nfeval);
      typedef void cost_temple(double *p, double *result, int np);
      typedef void g_temple(double *p, double *result);
      typedef void h_temple(double *p, double *result);
      void SOLNP_PLUS(
          SOLNPProb *prob,
          SOLNPSettings *stgs,
          SOLNPSol *sol,
          SOLNPInfo *info,
          cost_temple *cost_fun,
          g_temple *grad_fun,
          h_temple *hess_fun,
          solnp_float *l,
          solnp_float *h);
#ifdef __cplusplus
}
#endif
#endif

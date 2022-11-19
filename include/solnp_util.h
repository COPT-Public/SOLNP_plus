#ifndef UTIL_H_GUARD
#define UTIL_H_GUARD

#ifdef __cplusplus
extern "C" {
#endif

#include "solnp.h"
#include <stdlib.h>
#include <stdio.h>

#if (defined NOTIMER)
typedef void *SOLNP(timer);

#elif(defined _WIN32 || defined _WIN64 || defined _WINDLL)

#include <windows.h>
typedef struct SOLNP(timer) 
{
      LARGE_INTEGER tic;
      LARGE_INTEGER toc;
      LARGE_INTEGER freq;
} SOLNP(timer);

#elif(defined __APPLE__)

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

void SOLNP(tic)(SOLNP(timer) *t);
solnp_float SOLNP(toc)(SOLNP(timer) *t);
solnp_float SOLNP(str_toc)(char *str, SOLNP(timer) *t);
solnp_float SOLNP(tocq)(SOLNP(timer) *t);


SOLNPConstraint *malloc_constriant(solnp_int n, solnp_int nic);
void set_default_settings(SOLNPSettings *stgs);
SOLNPCost *malloc_cost(solnp_int nec, solnp_int nic, void *func);
void free_sol(SOLNPSol *sol);

#ifdef __cplusplus
}
#endif
#endif

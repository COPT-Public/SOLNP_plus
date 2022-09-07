#ifndef GLB_H_GUARD
#define GLB_H_GUARD

#ifdef __cplusplus
extern "C" {
#endif

#include <math.h>

//#define DLONG

#ifndef SOLNP
#define SOLNP(x) solnp_##x
#endif

#ifndef SUBNP
#define SUBNP(x) subnp_##x
#endif

/* SOLNP VERSION NUMBER ----------------------------------------------    */
#define SOLNP_VERSION                                                            \
  ("1.0.0") /* string literals automatically null-terminated */



#ifdef MATLAB_MEX_FILE
#include "mex.h"
#define solnp_printf mexPrintf
#define _solnp_free mxFree
#define _solnp_malloc mxMalloc
#define _solnp_calloc mxCalloc
#define _solnp_realloc mxRealloc
#elif defined PYTHON
#include <Python.h>
#include <stdlib.h>
#define solnp_printf(...)                                                            \
{                                                                                               \
      PyGILState_STATE gilstate = PyGILState_Ensure();       \
      PySys_WriteStdout(__VA_ARGS__);                                 \
      PyGILState_Release(gilstate);                                          \
}
#define _solnp_free free
#define _solnp_malloc malloc
#define _solnp_calloc calloc
#define _solnp_realloc realloc
#else
#include <stdio.h>
#include <stdlib.h>
#define solnp_printf printf
#define _solnp_free free
#define _solnp_malloc malloc
#define _solnp_calloc calloc
#define _solnp_realloc realloc
#endif

#define solnp_free(x)   \
      _solnp_free(x);          \
      x = SOLNP_NULL 
#define solnp_malloc(x) _solnp_malloc(x)
#define solnp_calloc(x, y) _solnp_calloc(x, y)
#define solnp_realloc(x, y) _solnp_realloc(x, y)

// //#ifdef DLONG
// //#ifdef _WIN64
// //typedef __int64 solnp_int; 
// //#else
// //typedef long solnp_int;
// //#endif
// //#else
// //typedef int solnp_int;
// //#endif
// typedef int solnp_int;


#ifdef DLONG
/*#ifdef _WIN64
#include <stdint.h>
typedef int64_t solnp_int;
#else
typedef long solnp_int;
#endif
*/
typedef long long solnp_int;
#else
typedef int solnp_int;
#endif


#ifndef SFLOAT
typedef double solnp_float;
#ifndef NAN
#define NAN ((solnp_float)0x7ff8000000000000)
#endif
#ifndef INFINITY
#define INFINITY NAN
#endif
#else
typedef float solnp_float;
#ifndef NAN
#define NAN ((float)0x7fc00000)
#endif
#ifndef INFINITY
#define INFINITY NAN
#endif
#endif

#define SOLNP_NULL 0

#ifndef MAX
#define MAX(a, b) (((a) > (b)) ? (a) : (b))
#endif

#ifndef MIN
#define MIN(a, b) (((a) < (b)) ? (a) : (b))
#endif

#ifndef ABS
#define ABS(x) (((x) < 0) ? -(x) : (x))
#endif

#ifndef POWF
#ifdef SFLOAT
#define POWF powf
#else
#define POWF pow
#endif
#endif

#ifndef SQRTF
#ifdef SFLOAT
#define SQRTF sqrtf
#else
#define SQRTF sqrt
#endif
#endif

#if EXTRA_VERBOSE > 1
#if (defined _WIN32 || defined _WIN64 || defined _WINDLL)
#define __func__ __FUNCTION__
#endif
#define DEBUG_FUNC solnp_printf("IN function: %s, time: %4f ms, file: %s, line: %i\n", __func__,  SOLNP(tocq)(&global_timer), __FILE__, __LINE__);
#define RETURN 
      solnp_printf("EXIT function: %s, time: %4f ms, file: %s, line: %i\n", __func__, SOLNP(tocq)(&global_timer), __FILE__, __LINE__);   \
      return
#else
#define DEBUG_FUNC
#define RETURN return
#endif

#define EPS_TOL (1E-18)
#define EPS (1E-8) // for condition number in subnp
#define SAFEDIV_POS(X, Y) ((Y) < EPS_TOL ? ((X) / EPS_TOL) : (X) / (Y))

#define CONVERGED_INTERVAL (1)
#define INDETERMINATE_TOL (1e-9)

#ifdef __cplusplus
}
#endif
#endif

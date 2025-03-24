#pragma once
#ifndef LINALG_H_GUARD
#define LINALG_H_GUARD
#define PI 3.14159265358979323846
#ifdef __cplusplus
extern "C"
{
#endif

#include "osqp.h"
#include "solnp.h"
#include <math.h>
#include <stdlib.h>
	// #include "osqp.h"

	solnp_float SOLNP(vec_mean)(
		solnp_float *x,
		solnp_int len);

	void SOLNP(set_as_scaled_array)(
		solnp_float *x,
		const solnp_float *a,
		const solnp_float b,
		solnp_int len);

	void SOLNP(set_as_sqrt)(
		solnp_float *x,
		const solnp_float *v,
		solnp_int len);

	void SOLNP(set_as_sq)(
		solnp_float *x,
		const solnp_float *v,
		solnp_int len);

	void SOLNP(scale_array)(
		solnp_float *a,
		const solnp_float b,
		solnp_int len);

	solnp_float SOLNP(dot)(
		const solnp_float *x,
		const solnp_float *y,
		solnp_int len);

	solnp_float SOLNP(norm_sq)(
		const solnp_float *v,
		solnp_int len);

	solnp_float SOLNP(norm_1)(
		const solnp_float *x,
		const solnp_int len);

	solnp_float SOLNP(cone_norm_1)(
		const solnp_float *x,
		const solnp_int len);

	solnp_float SOLNP(norm)(
		const solnp_float *v,
		solnp_int len);

	solnp_float SOLNP(norm_inf)(
		const solnp_float *a,
		solnp_int len);

	void SOLNP(add_array)(
		solnp_float *a,
		const solnp_float b,
		solnp_int len);

	void SOLNP(add_scaled_array)(
		solnp_float *a,
		const solnp_float *b,
		solnp_int n,
		const solnp_float sc);

	solnp_float SOLNP(norm_diff)(
		const solnp_float *a,
		const solnp_float *b,
		solnp_int len);

	solnp_float SOLNP(norm_inf_diff)(
		const solnp_float *a,
		const solnp_float *b,
		solnp_int len);

	void SOLNP(Ax)(
		solnp_float *ax,
		const solnp_float *a,
		const solnp_float *x,
		solnp_int m,
		solnp_int n);

	/* Rank 1 update of matrix :h  =  h + alpha * x x^T*/
	void SOLNP(rank1update)(
		solnp_int n,
		solnp_float *h,
		solnp_float alpha,
		solnp_float *x);

	void SOLNP(transpose)(
		const solnp_int m,
		const solnp_int n,
		const solnp_float *A,
		solnp_float *AT);

	void SOLNP(AB)(
		solnp_float *c,
		const solnp_float *a,
		const solnp_float *b,
		solnp_int m,
		solnp_int n,
		solnp_int p);

	solnp_float SOLNP(min)(solnp_float *a, solnp_int len);
	solnp_float SOLNP(max)(solnp_float *a, solnp_int len);

	solnp_int countA_sys(
		solnp_int m,
		solnp_int n,
		solnp_float *A);
	solnp_int countA(
		solnp_int m,
		solnp_int n,
		solnp_float *A);
	void calculate_csc_sys(
		c_int m,
		c_int n,
		c_float *A,
		c_float *A_x,
		c_int *A_i,
		c_int *A_p);
	void calculate_csc(
		c_int m,
		c_int n,
		c_float *A,
		c_float *A_x,
		c_int *A_i,
		c_int *A_p);
	void max_kelement(
		solnp_float *array,
		solnp_int len,
		solnp_int k,
		solnp_int *output);
	solnp_float Uniform_dis(solnp_float range);
	void Gaussian(
		solnp_float mean,
		solnp_float stddev,
		solnp_int length,
		solnp_float *x);
	void Uniform_sphere(
		solnp_float *x,
		solnp_int dim,
		solnp_float radius);
#ifdef __cplusplus
}
#endif
#endif

#include "linalg.h"
#include <math.h>
#include "mkl.h"

/* y = mean(x) */
solnp_float SOLNP(vec_mean)
(
      solnp_float *x,
      solnp_int len
)
{     
      if(x == SOLNP_NULL || len <= 0){
            return 0;
      }
      if(len<=0 || x == SOLNP_NULL){
            printf("invalid SOLNP(vec_mean) parameter");
            return -1;
      }
      solnp_float y=0;

      for(int i=0;i<len;i++){
            y += x[i];
      }

      return y/len;
}


/* x = b*a */
void SOLNP(set_as_scaled_array)
(
      solnp_float *x, 
      const solnp_float *a, 
      const solnp_float b,
      solnp_int len
) 
{     
      if(x == SOLNP_NULL || a == SOLNP_NULL || len <= 0){
            return;
      }

      solnp_int i;
      for (i = 0; i < len; ++i)
      {
            x[i] = b * a[i];
      }
}

/* x = sqrt(v) */
void SOLNP(set_as_sqrt)
(
      solnp_float *x, 
      const solnp_float *v, 
      solnp_int len
)
{
      if(x == SOLNP_NULL || v == SOLNP_NULL || len <= 0){
            return;
      }

      solnp_int i; 
      for (i = 0; i < len; ++i)
      {
            x[i] = SQRTF(v[i]); 
      }
}

/* x = v.^2 */
void SOLNP(set_as_sq)
(
      solnp_float *x, 
      const solnp_float *v, 
      solnp_int len
)
{
      if(x == SOLNP_NULL || v == SOLNP_NULL || len <= 0){
            return;
      }

      solnp_int i; 
      for (i = 0; i < len; ++i)
      {
            x[i] = v[i]*v[i]; 
      }
}  

/* a *= b */
void SOLNP(scale_array)
(
      solnp_float *a, 
      const solnp_float b, 
      solnp_int len
) 
{
      if(a == SOLNP_NULL || len <= 0){
            return;
      }

      solnp_int i;
      for (i = 0; i < len; ++i)
      {
            a[i] *= b;
      }
}

/* x'*y */
solnp_float SOLNP(dot)
(
      const solnp_float *x, 
      const solnp_float *y, 
      solnp_int len
) 
{     
      if(x == SOLNP_NULL || y == SOLNP_NULL || len <= 0){
            return 0;
      }

      solnp_int i;
      solnp_float ip = 0.0;
      for (i = 0; i < len; ++i) 
      {
            ip += x[i] * y[i];
      }
      return ip;
}

/* ||v||_2^2 */
solnp_float SOLNP(norm_sq)
(
      const solnp_float *v, 
      solnp_int len
) 
{    
      if(v == SOLNP_NULL || len <= 0){
            return 0;
      }

      solnp_int i;
      solnp_float nmsq = 0.0;
      for (i = 0; i < len; ++i) 
      {
            nmsq += v[i] * v[i];
      }
      return nmsq;
}

/* ||v||_2 */
solnp_float SOLNP(norm)
(
      const solnp_float *v, 
      solnp_int len
) 
{      
      if(v == SOLNP_NULL || len <= 0){
            return 0;
      }

      return SQRTF(SOLNP(norm_sq)(v, len));
}

/* ||x||_1 */
solnp_float SOLNP(norm_1)
(
	const solnp_float *x, 
	const solnp_int len
)
{     
      if(x == SOLNP_NULL || len <= 0){
            return 0;
      }

      solnp_float result = 0;
      for(int i=0; i<len; i++){
            result += ABS(x[i]);
      }
      return result;
}

/*the absolute value of the largest component of x*/
solnp_float SOLNP(cone_norm_1)
(
	const solnp_float *x,
	const solnp_int len
)
{
      if(x == SOLNP_NULL || len <= 0){
            return 0;
      }

      solnp_int i;
      solnp_float tmp; 
      solnp_float max = 0.0;
      for (i = 0; i < len; ++i) 
      {
            tmp = x[i];
            if (tmp > max) 
            {
                  max = tmp;
            }
      }
      return ABS(max);
}

/* max(|v|) */
solnp_float SOLNP(norm_inf)
(
      const solnp_float *a, 
      solnp_int len
) 
{
      if(a == SOLNP_NULL || len <= 0){
            return 0;
      }

      solnp_int i;
      solnp_float tmp; 
      solnp_float max = 0.0;
      for (i = 0; i < len; ++i) 
      {
            tmp = ABS(a[i]);
            if (tmp > max) 
            {
                  max = tmp;
            }
      }
      return max;
}

/* a .+= b */
void SOLNP(add_array)
(
      solnp_float *a, 
      const solnp_float b, 
      solnp_int len
)
{
      if(a == SOLNP_NULL || len <= 0){
            return;
      }

      solnp_int i; 
      for (i = 0; i <len; ++i) 
      {
            a[i] += b;
      }
} 

/* saxpy a += sc*b */
void SOLNP(add_scaled_array)
( 
      solnp_float *a, 
      const solnp_float *b, 
      solnp_int len,
      const solnp_float sc
) 
{
      if(a == SOLNP_NULL || b == SOLNP_NULL || len <= 0){
            return;
      }

      solnp_int i;
      for (i = 0; i < len; ++i) 
      {
            a[i] += sc * b[i];
      }
}

/* ||a-b||_2^2 */
solnp_float SOLNP(norm_diff)
(
      const solnp_float *a, 
      const solnp_float *b, 
      solnp_int len
) 
{
      if(a == SOLNP_NULL || b == SOLNP_NULL || len <= 0){
            return 0;
      }

      solnp_int i;
      solnp_float tmp;
      solnp_float nm_diff = 0.0; 
      for (i = 0; i < len; ++i) 
      {
            tmp = (a[i] - b[i]);
            nm_diff += tmp * tmp;
      }
      return SQRTF(nm_diff);
}

/* max(|a-b|) */
solnp_float SOLNP(norm_inf_diff)
(
      const solnp_float *a, 
      const solnp_float *b,
      solnp_int len
) 
{
      if(a == SOLNP_NULL || b == SOLNP_NULL || len <= 0){
            return 0;
      }

      solnp_int i;
      solnp_float tmp; 
      solnp_float max = 0.0;
      for (i = 0; i < len; ++i) 
      {
            tmp = ABS(a[i] - b[i]);
            if (tmp > max) 
            {
                  max = tmp;
            }
      }
      return max;
}

/* Ax where A \in R^m*n and x \in R^n */
/* ax = Ax */
void SOLNP(Ax)
(
      solnp_float *ax,
	const solnp_float *A,
	const solnp_float *x,
	solnp_int m,
	solnp_int n
)
{
      solnp_int i;
      SOLNP(set_as_scaled_array)(ax, A, x[0], m);
      for (i = 1; i < n; i++) {
            SOLNP(add_scaled_array)(ax, &(A[i*m]), m, x[i]);
      }

      return;
}


/* Rank 1 update of matrix :h  =  h + alpha * x x^T*/
void SOLNP(rank1update)
(
    solnp_int n,
    solnp_float* h,
    solnp_float alpha,
    solnp_float* x
) 
{
      solnp_int i,j;

      for(i=0; i<n; i++){//col index
            for(j=0; j<n; j++){//row index
            //element h(j,i)
                  h[i*n + j] += alpha * x[i]*x[j];
            }
      }
}
// /* Cholesky Decomposition: chol(h+ mu* diag(dx))*/
// void SOLNP(chol)
// (
//     solnp_int n,
//     SOLNPMatrix* h,
//     solnp_float* dx,
//     solnp_float** result,
//     solnp_float mu
// ) 
// {

// }


solnp_float SOLNP(max)(solnp_float* a, solnp_int len) {
    solnp_int i;
    solnp_float m = -INFINITY;
    for (i = 0; i < len; i++) {
        if (m < a[i]) {
            m = a[i];
        }
    }
    return m;
}
solnp_float SOLNP(min)(solnp_float* a, solnp_int len) {
    solnp_int i;
    solnp_float m = INFINITY;
    for (i = 0; i < len; i++) {
        if (m > a[i]) {
            m = a[i];
        }
    }
    return m;
}

// col major A in R^m*n, AT in R^n*m
void SOLNP(transpose)
(
      const solnp_int m,
      const solnp_int n,
      const solnp_float *A,
      solnp_float *AT
)
{
      solnp_int i;
      solnp_int j;

      for(i=0; i<n; i++){
            for(j=0; j<m; j++){
                  AT[i + j*n] = A[j + m*i];
            }
      }
}


/* C = AB where A \in R^m*n and B \in R^n*p */
void SOLNP(AB)
(
      solnp_float *c,
	const solnp_float *a,
	const solnp_float *b,
	solnp_int m,
	solnp_int n,
	solnp_int p
)
{
      solnp_int i;
      for(i=0; i<p; i++){
            SOLNP(Ax)(&c[i*m], a, &b[i*n], m, n);
      }
}

c_int countA_sys
(
    solnp_int m,
    solnp_int n,
    solnp_float* A
) {
    c_int i, j, count = 0;
    for (j = 0; j < n; j++) {
        for (i = 0; i < j + 1; i++) {
            if (A[i + j * m] != 0.) {
                count = count + 1;
            }
        }
    }return count;
}
c_int countA
(
    solnp_int m,
    solnp_int n,
    solnp_float* A
) {
    c_int i, j, count = 0;
    for (j = 0; j < n; j++) {
        for (i = 0; i < m; i++) {
            if (A[i + j * m] != 0.) {
                count = count + 1;
            }
        }count = count + 1;
    }return count;
}

void calculate_csc_sys
(
    solnp_int m,
    solnp_int n,
    solnp_float* A,
    c_float* A_x,
    c_int* A_i,
    c_int* A_p
) {
    c_int i, j, count = 0;
    A_p[0] = count;
    for (j = 0; j < n; j++) {
        for (i = 0; i < j + 1; i++) {
            if (A[i + j * m] != 0) {
                A_x[count] = A[i + j * m];
                A_i[count] = i;
                count = count + 1;
            }
        }A_p[j + 1] = count;
    }
}
void calculate_csc
(
    solnp_int m,
    solnp_int n,
    solnp_float* A,
    c_float* A_x,
    c_int* A_i,
    c_int* A_p
) {
    c_int i, j, count = 0;
    A_p[0] = count;
    for (j = 0; j < n; j++) {
        for (i = 0; i < m; i++) {
            if (A[i + j * m] != 0) {
                A_x[count] = A[i + j * m];
                A_i[count] = i;
                count = count + 1;
            }
        }A_x[count] = 1;
        A_i[count] = j + m;
        count = count + 1;
        A_p[j + 1] = count;
    }
}
/*
solnp_float SOLNP(cond)
(
    solnp_int m,
    solnp_int n,
    solnp_float* a
) 
{
}*/
/*

/*
 * Get number of nonzero elements in a dense matrix
 *//*
static solnp_int nonzero_elements(const solnp_float * in_matrix, const solnp_int nrows, const solnp_int ncols) {
	int i_row, i_col;
	solnp_int num_nonzero = 0;
	for (i_row = 0; i_row < nrows; i_row++) {
		for (i_col = 0; i_col < ncols; i_col++) {
			if (in_matrix[i_col + i_row * ncols] != 0.0) {
				num_nonzero++;
			}
		}
	}
	return num_nonzero;
}

static solnp_int arr_ind(const solnp_int i_col, const solnp_int i_row, const solnp_int nrows, const solnp_int ncols, const solnp_int format) {
	return (format == RowMajor) ? (i_col + i_row * ncols) : (i_row + i_col * nrows);
}

void *dense_to_csc_matrix(solnp_float * in_matrix, 
                        const solnp_int nrows, 
                        const solnp_int ncols, 
                        const solnp_int format, 
                        solnp_int **p, 
                        solnp_int **i, 
                        solnp_float **x) 
      {
	solnp_int i_row, i_col, ind_mat, ind_val = 0, num_col_nnz = 0;
	solnp_float *values;
	solnp_int *rows, *col_nnz;
	const solnp_int nnz = nonzero_elements(in_matrix, nrows, ncols);
	values = (solnp_float *)solnp_malloc(sizeof(solnp_float) * nnz);
	rows = (solnp_int *)solnp_malloc(sizeof(solnp_int) * nnz);
	col_nnz = (solnp_int *)solnp_malloc(sizeof(solnp_int) * (ncols + 1));

	// Fill values
	col_nnz[0] = (solnp_int)0;
	for (i_col = 0; i_col < ncols; i_col++) {
		num_col_nnz = 0;
		for (i_row = 0; i_row < nrows; i_row++) {
			ind_mat = arr_ind(i_col, i_row, nrows, ncols, format);
			if (in_matrix[ind_mat] != 0.0) {
				values[ind_val] = in_matrix[ind_mat];
				rows[ind_val] = i_row;
				ind_val++;
				num_col_nnz++;
			}
		}
		col_nnz[i_col + 1] = col_nnz[i_col] + num_col_nnz;
	}

	// Create CSR structure
      *p = col_nnz;
      *i = rows;
      *x = values;
}*/
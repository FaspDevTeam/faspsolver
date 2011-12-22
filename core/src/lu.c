/*! \file lu.c
 *  \brief LU decomposition and direct solve for dense matrix
 */

#include <math.h> 

#include "fasp.h"
#include "fasp_functs.h"

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

/**
 * \fn int fasp_smat_lu_decomp (double *A, int pivot[], int n) 
 * \brief LU decomposition of A usind Doolittle's method
 *
 * \param *A      pointer to the full matrix
 * \param pivot   pivoting positions
 * \param n       size of matrix A
 * \return        1 if succeed, 0 if fail
 *
 * Use Doolittle's method to decompose the n x n matrix A into a unit 
 * lower triangular matrix L and an upper triangular matrix U such that 
 * A = LU. The matrices L and U replace the matrix A. The diagonal elements 
 * of L are 1 and are not stored. 
 * 
 * The Doolittle method with partial pivoting is:  Determine the pivot 
 * row and interchange the current row with the pivot row, then assuming 
 * that row k is the current row, k = 0, ..., n - 1 evaluate in order the
 * following pair of expressions
 *       U[k][j] = A[k][j] - (L[k][0]*U[0][j] + ... + L[k][k-1]*U[k-1][j]) 
 *                                 for j = k, k+1, ... , n-1                
 *       L[i][k] = (A[i][k] - (L[i][0]*U[0][k] + . + L[i][k-1]*U[k-1][k]))    
 *                          / U[k][k]                                        
 *                                 for i = k+1, ... , n-1.              
 *
 * \author Xuehai Huang
 * \date 04/02/2009 
 */
int fasp_smat_lu_decomp (double *A, 
												 int pivot[], 
												 int n) 
{
	int i, j, k;
	double *p_k=NULL, *p_row=NULL, *p_col=NULL;
	double max;
	
	/* For each row and column, k = 0, ..., n-1, */
	for (k = 0, p_k = A; k < n; p_k += n, k++) {
		
		// find the pivot row
		pivot[k] = k;
		max = fabs( *(p_k + k) );
		for (j = k + 1, p_row = p_k + n; j < n; ++j, p_row += n) {
			if ( max < fabs(*(p_row + k)) ) {
				max = fabs(*(p_row + k));
				pivot[k] = j;
				p_col = p_row;
			}
		}
		
		// if the pivot row differs from the current row, interchange the two rows.		
		if (pivot[k] != k)
			for (j = 0; j < n; ++j) {
				max = *(p_k + j);
				*(p_k + j) = *(p_col + j);
				*(p_col + j) = max;
			}
		
		// if the matrix is singular, return error
		if ( fabs( *(p_k + k) ) < SMALLREAL ) return -1;
		
		// otherwise find the lower triangular matrix elements for column k. 		
		for (i = k+1, p_row = p_k + n; i < n; p_row += n, ++i) {
			*(p_row + k) /= *(p_k + k);
		}  
		
		// update remaining matrix		
		for (i = k+1, p_row = p_k + n; i < n; p_row += n, ++i)
			for (j = k+1; j < n; ++j)
				*(p_row + j) -= *(p_row + k) * *(p_k + j);
		
	}
	
	return 0;
}

/**
 * \fn int fasp_smat_lu_solve (double *A, double b[], int pivot[], double x[], int n)
 * \brief Solving Ax=b using LU decomposition
 *
 * \param *A    pointer to the full matrix
 * \param b     right hand side array 
 * \param pivot pivoting positions
 * \param x     solution array 
 * \param n     size of matrix A
 * \return      0 if succeed, else if failed
 *
 * This routine uses Doolittle's method to solve the linear equation Ax = b.
 * This routine is called after the matrix A has been decomposed into a product
 * of a unit lower triangular matrix L and an upper triangular matrix U with 
 * pivoting. The solution proceeds by solving the linear equation Ly = b for y 
 * and subsequently solving the linear equation Ux = y for x.
 *
 * \author Xuehai Huang
 * \date 04/02/2009
 */   
int fasp_smat_lu_solve (double *A, 
												double b[], 
												int pivot[], 
												double x[], 
												int n)
{
	int i, k;
	double *p_k;
	double dum;
	
	/* solve Ly = b	*/
	for (k = 0, p_k = A; k < n; p_k += n, k++) {
		if (pivot[k] != k) {dum = b[k]; b[k] = b[pivot[k]]; b[pivot[k]] = dum; }
		x[k] = b[k];
		for (i = 0; i < k; ++i) x[k] -= x[i] * *(p_k + i);
	}
	
	/* solve Ux = y */
	for (k = n-1, p_k = A + n*(n-1); k >= 0; k--, p_k -= n) {
		if (pivot[k] != k) {dum = b[k]; b[k] = b[pivot[k]]; b[pivot[k]] = dum; }
		for (i = k + 1; i < n; ++i) x[k] -= x[i] * *(p_k + i);
		if (*(p_k + k) == 0.0) return -1;
		x[k] /= *(p_k + k);
	}
  
	return 0;
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/

/** 
 //A simple test example can be written as the following
 
 int main (int argc, const char * argv[]) 
 {
 double A[3][3] = {{0.0, 1.0, 4.0},
 {4.0, 1.0, 0.0},
 {1.0, 4.0, 1.0}};  
 
 double b[3] = {1, 1, 1}, x[3];
 
 int pivot[3];          
 
 int ret, i, j;
 
 ret = lu_decomp(&A[0][0], pivot, 3); // LU decomposition
 
 ret = lu_solve(&A[0][0], b, pivot, x, 3); // Solve decomposed Ax=b
 
 return 1;
 }
 
 */

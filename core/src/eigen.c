/*! \file eigen.c
 *  \brief simple subroutines for compute the extreme eigenvalues.
 */

#include <math.h> 

#include "fasp.h"
#include "fasp_functs.h"

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

/**
 * \fn REAL fasp_dcsr_eig (dCSRmat *A, const REAL tol, const INT maxit)
 *
 * \brief Approximate the largest eigenvalue of A by the power method
 *
 * \param *A     pointer to the dCSRmat
 * \param tol    tolerance
 * \param maxit	 max number of iterations
 *
 * \return       the largest eigenvalue
 *
 * \author Xiaozhe Hu
 * \date 01/25/2011 
 */
REAL fasp_dcsr_eig (dCSRmat *A, 
                    const REAL tol, 
                    const INT maxit)
{
	REAL eigenvalue=0.0, temp=1.0;
	
	dvector x, y;
	fasp_dvec_alloc(A->row, &x); 
	fasp_dvec_rand(A->row,&x); 
	fasp_blas_array_ax(A->row, 1.0/fasp_blas_dvec_norm2(&x), x.val);
	fasp_dvec_alloc(A->row, &y);
	
	REAL L2_norm_y;
	unsigned INT i;
	
	for (i=0; i<maxit; i++)
	{
		// y = Ax;
		fasp_blas_dcsr_mxv(A, x.val, y.val);
		
		// y/||y||
		L2_norm_y = fasp_blas_dvec_norm2(&y);
		fasp_blas_array_ax(A->row, 1.0/L2_norm_y, y.val);
		
		// eigenvalue = y'Ay;
		eigenvalue = fasp_blas_dcsr_vmv(A, y.val, y.val);
		
		// convergence test
		if ((ABS(eigenvalue - temp)/ABS(temp))<tol) goto FINISHED;
		
		// 
		fasp_dvec_cp(&y, &x);
		temp = eigenvalue;
	}
    
FINISHED:
	fasp_dvec_free(&x);
	fasp_dvec_free(&y);
    
	return eigenvalue;
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/

/*! \file amg_solve.c
 *  \brief Algebraic multigrid iterations: SOLVE phase. 
 */

#include <time.h>

#include "fasp.h"
#include "fasp_functs.h"

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

/**
 * \fn INT amg_solve (AMG_data *mgl, AMG_param *param)
 *
 * \brief AMG -- SOLVE phase
 *
 * \param *mgl    pointer to AMG_data data
 * \param *param  pointer to AMG parameters
 *
 * \return        iteration number if succeed
 *
 * \note Solve Ax=b using multigrid method. This is SOLVE phase only and is 
 * independent of SETUP method used! Should be called after multigrid hierarchy 
 * has been generated!
 *
 * \author Xuehai Huang, Chensong Zhang
 * \date 04/02/2010
 */
INT fasp_amg_solve (AMG_data *mgl, 
                    AMG_param *param)
{
	dCSRmat  *ptrA=&mgl[0].A;
	dvector  *b=&mgl[0].b, *x=&mgl[0].x, *r=&mgl[0].w; 
    
	const INT   MaxIt = param->max_iter; 
	const INT   print_level = param->print_level;
	const INT   m=ptrA->row, n=ptrA->col, nnz=ptrA->nnz;	
	const REAL  tol = param->tol;
	const REAL  sumb = fasp_blas_dvec_norm2(b); // L2norm(b)	
	
    // local variables
	unsigned INT iter=0;
	REAL relres1=BIGREAL, absres0=BIGREAL, absres, factor;		
	clock_t solve_start=clock();
	
	if (print_level>=PRINT_MOST) printf("fasp_amg_solve: nr=%d, nc=%d, nnz=%d\n",m,n,nnz);	
	
	while ((++iter <= MaxIt) & (sumb > SMALLREAL)) // MG solver here
	{	
		// Call one multigrid cycle
		fasp_solver_mgcycle(mgl, param); 
        // This is the non-recursive MG cycle subroutine.
        // If you wish to call the recursive version instead, replace it with:
        //     fasp_solver_mgrecur(mgl, param, 0);         
		
		// Form residual r = b-A*x		
		fasp_dvec_cp(b,r); fasp_blas_dcsr_aAxpy(-1.0,ptrA,x->val,r->val);		
		
        // Compute norms of r and convergence factor
		absres  = fasp_blas_dvec_norm2(r); // residual ||r||
		relres1 = absres/sumb;             // relative residual ||r||/||b||
		factor  = absres/absres0;          // contraction factor
        absres0 = absres;                  // prepare for next iteration

		// Print iteration information if needed	
		print_itinfo(print_level, STOP_REL_RES, iter, relres1, absres, factor);
		
        // Check convergence
		if (relres1<tol) break;
    }
	
	if (print_level>PRINT_NONE) {
		if (iter>MaxIt)
			printf("Maximal iteration %d exceeded with relative residual %e.\n", 
                   MaxIt, relres1);
		else
			printf("Number of iterations = %d with relative residual %e.\n", 
                   iter, relres1);
		
		clock_t solve_end=clock();
		REAL solveduration = (REAL)(solve_end - solve_start)/(REAL)(CLOCKS_PER_SEC);
		printf("AMG solve costs %f seconds.\n", solveduration);
	}
	
	return iter;
}

/**
 * \fn INT fasp_amg_solve_amli(AMG_data *mgl, AMG_param *param)
 *
 * \brief AMLI -- SOLVE phase
 *
 * \param *mgl    pointer to AMG_data data
 * \param *param  pointer to AMG parameters
 *
 * \return        iteration number if succeed
 *
 * \note Solve Ax=b using multigrid method. Should be called after 
 * multigrid hierarchy has been generated!
 *
 * \author Xiaozhe Hu
 * \date 01/23/2011
 */
INT fasp_amg_solve_amli (AMG_data *mgl, 
                         AMG_param *param)
{
	dCSRmat     *ptrA=&mgl[0].A;
	dvector     *b=&mgl[0].b, *x=&mgl[0].x, *r=&mgl[0].w; 
    
	const INT   MaxIt = param->max_iter; 
	const INT   print_level = param->print_level;
	const INT   m=ptrA->row, n=ptrA->col, nnz=ptrA->nnz;	
	const REAL  tol = param->tol;
	const REAL  sumb = fasp_blas_dvec_norm2(b); // L2norm(b)	
	
    // local variables
	REAL relres1=BIGREAL, absres0=BIGREAL, absres, factor;		
	unsigned INT iter=0;
	clock_t solve_start=clock();
	
	if (print_level>=PRINT_MOST) printf("fasp_amli_solve: nr=%d, nc=%d, nnz=%d\n", m, n, nnz);
    
	while ((++iter <= MaxIt) & (sumb > SMALLREAL)) // MG solver here
	{	
		// Call one multigrid cycle
		fasp_solver_amli(mgl, param, 0); 
		
		// Form residual r = b-A*x		
		fasp_dvec_cp(b,r); 
		fasp_blas_dcsr_aAxpy(-1.0,ptrA,x->val,r->val);		
		
        // Compute norms of r and convergence factor
		absres  = fasp_blas_dvec_norm2(r); // residual ||r||
		relres1 = absres/sumb;             // relative residual ||r||/||b||
		factor  = absres/absres0;          // contraction factor
        absres0 = absres;                  // prepare for next iteration
		
		// Print iteration information if needed	
		print_itinfo(print_level, STOP_REL_RES, iter, relres1, absres, factor);
		
        // Check convergence
		if (relres1<tol) break;
    }
	
	if (print_level>PRINT_NONE) {
		if (iter>MaxIt)
			printf("Maximal iteration %d exceeded with relative residual %e.\n", MaxIt, relres1);
		else
			printf("Number of iterations = %d with relative residual %e.\n", iter, relres1);
		
		clock_t solve_end=clock();
		REAL solveduration = (REAL)(solve_end - solve_start)/(REAL)(CLOCKS_PER_SEC);
		printf("AMLI solve costs %f seconds.\n", solveduration);
	}
	
	return iter;
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/

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
	
	if (print_level>=PRINT_MOST) printf("fasp_amg_solve: nr=%d, nc=%d, nnz=%d\n", m, n, nnz);	
	
	while ((++iter <= MaxIt) & (sumb > SMALLREAL)) // MG solver here
	{	
		// Call one multigrid cycle
        // This is the non-recursive MG cycle subroutine.
        // If you wish to call the recursive version instead, replace it with:
        //     fasp_solver_mgrecur(mgl, param, 0);         
		fasp_solver_mgcycle(mgl, param); 
		
		// r = b-A*x		
		fasp_dvec_cp(b,r); fasp_blas_dcsr_aAxpy(-1.0,ptrA,x->val,r->val);		
		
        // compute norms of r and convergence factor
		absres  = fasp_blas_dvec_norm2(r); // residual ||r||
		relres1 = absres/sumb;             // relative residual ||r||/||b||
		factor  = absres/absres0;          // contraction factor
		
		// output iteration information if needed	
		print_itinfo(print_level, STOP_REL_RES, iter, relres1, absres, factor);
		
		if (relres1<tol) break; // early exit condition
		
		absres0 = absres;
	}
	
	if (print_level>PRINT_NONE) {
		if (iter>MaxIt)
			printf("Maximal iteration %d exceeded with relative residual %e.\n", MaxIt, relres1);
		else
			printf("Number of iterations = %d with relative residual %e.\n", iter, relres1);
		
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
		
		// r = b-A*x		
		fasp_dvec_cp(b,r); 
		fasp_blas_dcsr_aAxpy(-1.0,ptrA,x->val,r->val);		
		
        // compute norms of r and convergence factor
		absres  = fasp_blas_dvec_norm2(r); // residual ||r||
		relres1 = absres/sumb;             // relative residual ||r||/||b||
		factor  = absres/absres0;          // contraction factor
		
		// output iteration information if needed	
		print_itinfo(print_level, STOP_REL_RES, iter, relres1, absres, factor);
		
		if (relres1<tol) break; // early exit condition
		
		absres0 = absres;
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

/*---------------------------------omp----------------------------------------*/

/**
 * \fn int fasp_amg_solve_omp(AMG_data *mgl, AMG_param *param, int nthreads, int openmp_holds)
 * \brief AMG solve phase
 *
 * \param *mgl    pointer to AMG_data data
 * \param *param  pointer to AMG parameters
 * \param nthreads number of threads
 * \param openmp_holds threshold of parallelization
 *
 * Solve Ax=b using multigrid method.
 *
 * \note Should be called after multigrid hierarchy has been setup!
 *
 * \author Feng Chunsheng, Yue Xiaoqiang
 * \date 03/01/2011
 */
int fasp_amg_solve_omp (AMG_data *mgl, 
                        AMG_param *param, 
                        int nthreads, 
                        int openmp_holds)
{
	unsigned int iter=0;
	
#if FASP_USE_OPENMP
    param->max_iter= 100;
	const int     MaxIt = param->max_iter; 
	const int     print_level = param->print_level;
	const double  tol = param->tol;
	
	dCSRmat     *ptrA=&mgl[0].A;
	const int     m=ptrA->row, n=ptrA->col, nnz=ptrA->nnz;	
	dvector     *b=&mgl[0].b, *x=&mgl[0].x, *r=&mgl[0].w; 
	const double  sumb = fasp_blas_dvec_norm2_omp(b, nthreads, openmp_holds); // L2norm(b)	
	
	double relres1=BIGREAL, absres0=BIGREAL, absres, factor;		
	
	if (print_level>8) printf("fasp_amg_solve: nr=%d, nc=%d, nnz=%d  param->max_iter  %d MaxIt =%d \n", m, n, nnz, param->max_iter,MaxIt);
	
	double solve_start=omp_get_wtime();
	
	while ((++iter <= MaxIt) & (sumb > SMALLREAL)) // MG solver here
	{	
		// one multigrid cycle
		fasp_solver_mgcycle_omp(mgl, param, nthreads, openmp_holds); // fasp_solver_mgrecur(mgl, param, 0);
		
		// r = b-A*x		
		fasp_dvec_cp_omp(b,r,nthreads,openmp_holds);
		fasp_blas_dcsr_aAxpy_omp(-1.0,ptrA,x->val,r->val,nthreads,openmp_holds);
		
		absres  = fasp_blas_dvec_norm2_omp(r,nthreads,openmp_holds); // residual ||r||
		relres1 = absres/sumb;       // relative residual ||r||/||b||
		factor  = absres/absres0;    // contraction factor
		
		// output iteration information if needed	
		print_itinfo(print_level, STOP_REL_RES, iter, relres1, absres, factor);
		
		if (relres1<tol) break; // early exit condition
		
		absres0 = absres;
	}
	
	if (print_level>0) {
		if (iter>MaxIt)
			printf("Maximal iteration %d exceeded with relative residual %e.\n", MaxIt, relres1);
		else
			printf("Number of iterations = %d with relative residual %e.\n", iter, relres1);
		
		double solve_end=omp_get_wtime();
		double solveduration = solve_end - solve_start;
		printf("AMG solve costs %f seconds.\n", solveduration);
	}
	
#endif
	
	return iter;
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/

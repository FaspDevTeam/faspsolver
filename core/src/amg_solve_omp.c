/*! \file amg_solve_omp.c
 *  \brief Algebraic multigrid iterations: SOLVE phase. 
 */

#include <time.h>

#include "fasp.h"
#include "fasp_functs.h"

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

/**
 * \fn INT fasp_amg_solve_omp (AMG_data *mgl, AMG_param *param, INT nthreads, INT openmp_holds)
 * \brief AMG solve phase
 *
 * \param mgl    pointer to AMG_data data
 * \param param  pointer to AMG parameters
 * \param nthreads number of threads
 * \param openmp_holds threshold of parallelization
 *
 * Solve Ax=b using multigrid method.
 *
 * \note Should be called after multigrid hierarchy has been setup!
 *
 * \author FENG Chunsheng, Yue Xiaoqiang
 * \date 03/01/2011
 * \date Jan/11/2012 modified by FENG Chunsheng
 */
INT fasp_amg_solve_omp (AMG_data *mgl, 
                        AMG_param *param, 
                        INT nthreads, 
                        INT openmp_holds)
{
    dCSRmat  *ptrA=&mgl[0].A;
    dvector  *b=&mgl[0].b, *x=&mgl[0].x, *r=&mgl[0].w; 
    
    const INT   MaxIt = param->maxit; 
    const INT   print_level = param->print_level;
    const INT   m=ptrA->row, n=ptrA->col, nnz=ptrA->nnz;    
    const REAL  tol = param->tol;
    const REAL  sumb = fasp_blas_dvec_norm2_omp(b, nthreads, openmp_holds); // L2norm(b)    
    
    // local variables
    unsigned INT iter=0;
    REAL relres1=BIGREAL, absres0=BIGREAL, absres, factor;    
    REAL solve_start=omp_get_wtime();
    
    if (print_level>=PRINT_MOST) printf("fasp_amg_solve: nr=%d, nc=%d, nnz=%d\n",m,n,nnz);    
    
    while ((++iter <= MaxIt) & (sumb > SMALLREAL)) // MG solver here
        {    
            // Call one multigrid cycle
            fasp_solver_mgcycle_omp(mgl, param,nthreads,openmp_holds); 
            // This is the non-recursive MG cycle subroutine.
            // If you wish to call the recursive version instead, replace it with:
            //     fasp_solver_mgrecur(mgl, param, 0);         
    
            // Form residual r = b-A*x    
            fasp_dvec_cp_omp(b,r,nthreads,openmp_holds);
            fasp_blas_dcsr_aAxpy_omp(-1.0,ptrA,x->val,r->val,nthreads,openmp_holds);
    
            // Compute norms of r and convergence factor
            absres  = fasp_blas_dvec_norm2_omp(r,nthreads,openmp_holds); // residual ||r||
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
    
        REAL solve_end=omp_get_wtime();
        REAL solveduration = solve_end - solve_start;
        printf("AMG solve costs %f seconds.\n", solveduration);
    }
    
    return iter;
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/

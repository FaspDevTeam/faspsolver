/*! \file amg_solve.c
 *  \brief Algebraic multigrid iterations: SOLVE phase. 
 */

#include <time.h>
#include <omp.h>

#include "fasp.h"
#include "fasp_functs.h"

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

/**
 * \fn INT fasp_amg_solve (AMG_data *mgl, AMG_param *param)
 *
 * \brief AMG -- SOLVE phase
 *
 * \param mgl     Pointer to AMG_data data
 * \param param   Pointer to AMG parameters
 *
 * \return        Iteration number if succeed, ERROR otherwise
 *
 * \author Xuehai Huang, Chensong Zhang
 * \date 04/02/2010
 *
 * \note Solve Ax=b using multigrid method. This is SOLVE phase only and is 
 *       independent of SETUP method used! Should be called after multigrid 
 *       hierarchy has been generated!
 */

INT fasp_amg_solve (AMG_data *mgl, 
                    AMG_param *param)
{
    dCSRmat      *ptrA=&mgl[0].A;
    dvector      *b=&mgl[0].b, *x=&mgl[0].x, *r=&mgl[0].w; 
    
    const SHORT   print_level = param->print_level;
    const INT     MaxIt = param->maxit; 
    const REAL    tol = param->tol;
    const REAL    sumb = fasp_blas_dvec_norm2(b); // L2norm(b)    

	// local variables
	REAL solveduration;

#if FASP_USE_OPENMP
    REAL        solve_start=omp_get_wtime();
#else
    clock_t       solve_start=clock();
#endif

    REAL          relres1=BIGREAL, absres0=BIGREAL, absres, factor;    
    INT           iter=0;
    
#if DEBUG_MODE
    const INT     m=ptrA->row, n=ptrA->col, nnz=ptrA->nnz;    
    printf("### DEBUG: fasp_amg_solve ...... [Start]\n");
    printf("### DEBUG: nr=%d, nc=%d, nnz=%d\n", m, n, nnz);
#endif
    
    // MG solver here
    while ((++iter <= MaxIt) & (sumb > SMALLREAL)) {
        
#if TRUE
        // Call one multigrid cycle -- non recursive
        fasp_solver_mgcycle(mgl, param); 
#else
        // If you wish to call the recursive version instead, replace it with:
        fasp_solver_mgrecur(mgl, param, 0);         
#endif
        
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
            printf("Maximal iteration %d exceeded with relative residual %e.\n", 
                   MaxIt, relres1);
        else
            printf("Number of iterations = %d with relative residual %e.\n", 
                   iter, relres1);
#if FASP_USE_OPENMP
        REAL solve_end=omp_get_wtime();
        solveduration = (REAL)(solve_end - solve_start);
#else
        clock_t solve_end=clock();
        solveduration = (REAL)(solve_end - solve_start)/(REAL)(CLOCKS_PER_SEC);
#endif
        print_cputime("AMG solve",solveduration);
    }
    
#if DEBUG_MODE
    printf("### DEBUG: fasp_amg_solve ...... [Finish]\n");
#endif
    
    return iter;
}

/**
 * \fn INT fasp_amg_solve_amli (AMG_data *mgl, AMG_param *param)
 *
 * \brief AMLI -- SOLVE phase
 *
 * \param mgl     pointer to AMG_data data
 * \param param   pointer to AMG parameters
 *
 * \return        iteration number if succeed
 *
 * \note Solve Ax=b using multigrid method. Should be called after 
 *       multigrid hierarchy has been generated!
 *
 * \author Xiaozhe Hu
 * \date 01/23/2011
 */
INT fasp_amg_solve_amli (AMG_data *mgl, 
                         AMG_param *param)
{
    dCSRmat     *ptrA=&mgl[0].A;
    dvector     *b=&mgl[0].b, *x=&mgl[0].x, *r=&mgl[0].w; 
    
    const INT    MaxIt = param->maxit; 
    const SHORT  print_level = param->print_level;
    const REAL   tol = param->tol;
    const REAL   sumb = fasp_blas_dvec_norm2(b); // L2norm(b)    
    
    // local variables
    clock_t      solve_start=clock();
    REAL         relres1=BIGREAL, absres0=BIGREAL, absres, factor;    
    INT          iter=0;
    
#if DEBUG_MODE
    const INT    m=ptrA->row, n=ptrA->col, nnz=ptrA->nnz;    
    printf("### DEBUG: fasp_amg_solve_amli ...... [Start]\n");
    printf("### DEBUG: nr=%d, nc=%d, nnz=%d\n", m, n, nnz);
#endif
    
    // MG solver here
    while ((++iter <= MaxIt) & (sumb > SMALLREAL)) {

        // Call one AMLI cycle
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
            printf("Maximal iteration %d exceeded with relative residual %e.\n",
                   MaxIt, relres1);
        else
            printf("Number of iterations = %d with relative residual %e.\n",
                   iter, relres1);
    
        clock_t solve_end=clock();
        REAL solveduration = (REAL)(solve_end - solve_start)/(REAL)(CLOCKS_PER_SEC);
        print_cputime("AMLI solve",solveduration);
    }
    
#if DEBUG_MODE
    printf("### DEBUG: fasp_amg_solve_amli ...... [Finish]\n");
#endif
    
    return iter;
}

/**
 * \fn INT fasp_amg_solve_nl_amli (AMG_data *mgl, AMG_param *param)
 *
 * \brief Nonlinear AMLI --- SOLVE phase
 *
 * \param mgl     pointer to AMG_data data
 * \param param   pointer to AMG parameters
 *
 * \return        iteration number if succeed
 *
 * \note Solve Ax=b using nonlinear AMLI-cycle method. Should be called after 
 *       multigrid hierarchy has been setup!
 *
 * \author Xiaozhe Hu
 * \date 04/30/2011
 */
INT fasp_amg_solve_nl_amli (AMG_data *mgl, 
                            AMG_param *param)
{
    dCSRmat      *ptrA=&mgl[0].A;
    dvector      *b=&mgl[0].b, *x=&mgl[0].x, *r=&mgl[0].w; 
    
    const INT     MaxIt = param->maxit; 
    const SHORT   print_level = param->print_level;
    const REAL    tol = param->tol;
    const REAL    sumb = fasp_blas_dvec_norm2(b); // L2norm(b)    
    
    // local variables
    REAL          relres1=BIGREAL, absres0=BIGREAL, absres, factor;    
    INT           iter=0;
    
#if DEBUG_MODE
    const INT     m=ptrA->row, n=ptrA->col, nnz=ptrA->nnz;    
    printf("### DEBUG: fasp_amg_solve_nl_amli ...... [Start]\n");
    printf("### DEBUG: nr=%d, nc=%d, nnz=%d\n", m, n, nnz);
#endif
    
    clock_t solve_start=clock();
    
    while ((++iter <= MaxIt) & (sumb > SMALLREAL)) // MG solver here
        {    
            // one multigrid cycle
            fasp_solver_nl_amli(mgl, param, 0, mgl[0].num_levels); 
    
            // r = b-A*x    
            fasp_dvec_cp(b,r); 
            fasp_blas_dcsr_aAxpy(-1.0,ptrA,x->val,r->val);    
    
            absres  = fasp_blas_dvec_norm2(r); // residual ||r||
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
    
        clock_t solve_end=clock();
        REAL solveduration = (REAL)(solve_end - solve_start)/(REAL)(CLOCKS_PER_SEC);
        print_cputime("Nonlinear AMLI solve",solveduration);
    }
    
#if DEBUG_MODE
    printf("### DEBUG: fasp_amg_solve_nl_amli ...... [Finish]\n");
#endif
    
    return iter;
}

/**
 * \fn void fasp_famg_solve (AMG_data *mgl, AMG_param *param)
 *
 * \brief FMG -- SOLVE phase
 *
 * \param mgl     pointer to AMG_data data
 * \param param   pointer to AMG parameters
 *
 * \note Solve Ax=b using full multigrid method. This is SOLVE phase only and
 *       is independent of SETUP method used! Should be called after multigrid 
 *       hierarchy has been generated!
 *
 * \author Chensong Zhang
 * \date 01/10/2012
 */
void fasp_famg_solve (AMG_data *mgl, 
                      AMG_param *param)
{
    dCSRmat     *ptrA=&mgl[0].A;
    dvector     *b=&mgl[0].b, *x=&mgl[0].x, *r=&mgl[0].w; 
    
    const SHORT  print_level = param->print_level;
    const REAL   sumb = fasp_blas_dvec_norm2(b); // L2norm(b)    
    
    // local variables
    clock_t      solve_start=clock();
    REAL         relres1=BIGREAL, absres;    
    
#if DEBUG_MODE
    const INT    m=ptrA->row, n=ptrA->col, nnz=ptrA->nnz;    
    printf("### DEBUG: fasp_famg_solve ...... [Start]\n");
    printf("### DEBUG: nr=%d, nc=%d, nnz=%d\n", m, n, nnz);
#endif
    
    // Call one full multigrid cycle
    fasp_solver_fmgcycle(mgl, param); 
    
    // Form residual r = b-A*x    
    fasp_dvec_cp(b,r); fasp_blas_dcsr_aAxpy(-1.0,ptrA,x->val,r->val);    
    
    // Compute norms of r and convergence factor
    absres  = fasp_blas_dvec_norm2(r); // residual ||r||
    relres1 = absres/sumb;             // relative residual ||r||/||b||
    
    if (print_level>PRINT_NONE) {
        printf("FMG finishes with relative residual %e.\n", relres1);
    
        clock_t solve_end=clock();
        REAL solveduration = (REAL)(solve_end - solve_start)/(REAL)(CLOCKS_PER_SEC);
        print_cputime("FMG solve",solveduration);
    }
    
#if DEBUG_MODE
    printf("### DEBUG: fasp_famg_solve ...... [Finish]\n");
#endif
    
    return;
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/

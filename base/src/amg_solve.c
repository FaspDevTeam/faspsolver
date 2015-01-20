/*! \file amg_solve.c
 *
 *  \brief Algebraic multigrid iterations: SOLVE phase.
 *
 *  \note Solve Ax=b using multigrid method. This is SOLVE phase only and is
 *        independent of SETUP method used! Should be called after multigrid
 *        hierarchy has been generated!
 */

#include <time.h>

#include "fasp.h"
#include "fasp_functs.h"

#include "itsolver_util.inl"

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

/**
 * \fn INT fasp_amg_solve (AMG_data *mgl, AMG_param *param)
 *
 * \brief AMG -- SOLVE phase
 *
 * \param mgl    Pointer to AMG data: AMG_data
 * \param param  Pointer to AMG parameters: AMG_param
 *
 * \return       Iteration number if succeed, ERROR otherwise
 *
 * \author Xuehai Huang, Chensong Zhang
 * \date   04/02/2010
 *
 * Modified by Chensong 04/21/2013: Fix an output typo
 */
INT fasp_amg_solve (AMG_data *mgl,
                    AMG_param *param)
{
    dCSRmat      *ptrA = &mgl[0].A;
    dvector      *b = &mgl[0].b, *x = &mgl[0].x, *r = &mgl[0].w;
    
    const SHORT   print_level = param->print_level;
    const INT     MaxIt       = param->maxit;
    const REAL    tol         = param->tol;
    const REAL    sumb        = fasp_blas_dvec_norm2(b); // L2norm(b)
    
    // local variables
    REAL  solve_start, solve_end;
    REAL  relres1 = BIGREAL, absres0 = sumb, absres, factor;
    INT   iter = 0;

#if DEBUG_MODE
    printf("### DEBUG: %s ...... [Start]\n", __FUNCTION__);
    printf("### DEBUG: nrow = %d, ncol = %d, nnz = %d\n",
           mgl[0].A.row, mgl[0].A.col, mgl[0].A.nnz);
#endif

    fasp_gettime(&solve_start);
    
    // Print iteration information if needed
    print_itinfo(print_level, STOP_REL_RES, iter, 1.0, sumb, 0.0);
    
    // MG solver here
    while ( (++iter <= MaxIt) & (sumb > SMALLREAL) ) {
        
#if TRUE
        // Call one multigrid cycle -- non recursive version
        fasp_solver_mgcycle(mgl, param);
#else
        // Call one multigrid cycle -- recursive version
        fasp_solver_mgrecur(mgl, param, 0);
#endif
        
        // Form residual r = b - A*x
        fasp_dvec_cp(b, r);
        fasp_blas_dcsr_aAxpy(-1.0, ptrA, x->val, r->val);
        
        // Compute norms of r and convergence factor
        absres  = fasp_blas_dvec_norm2(r); // residual ||r||
        relres1 = absres/sumb;             // relative residual ||r||/||b||
        factor  = absres/absres0;          // contraction factor
        absres0 = absres;                  // prepare for next iteration
        
        // Print iteration information if needed
        print_itinfo(print_level, STOP_REL_RES, iter, relres1, absres, factor);
        
        // Check convergence
        if ( relres1 < tol ) break;
    }
    
    if ( print_level > PRINT_NONE ) {
        ITS_FINAL(iter, MaxIt, relres1);
        fasp_gettime(&solve_end);
        print_cputime("AMG solve",solve_end - solve_start);
    }
    
#if DEBUG_MODE
    printf("### DEBUG: %s ...... [Finish]\n", __FUNCTION__);
#endif
    
    return iter;
}

/**
 * \fn INT fasp_amg_solve_amli (AMG_data *mgl, AMG_param *param)
 *
 * \brief AMLI -- SOLVE phase
 *
 * \param mgl    Pointer to AMG data: AMG_data
 * \param param  Pointer to AMG parameters: AMG_param
 *
 * \return       Iteration number if succeed, ERROR otherwise
 *
 * \author Xiaozhe Hu
 * \date   01/23/2011
 *
 * \note AMLI polynomila computed by the best approximation of 1/x. Refer to Johannes K. Kraus, Panayot S. Vassilevski,
 *          Ludmil T. Zikatanov, "Polynomial of best uniform approximation to $x^{-1}$ and smoothing in two-leve
 *          methods", 2013.
 *
 *
 * Modified by Chensong 04/21/2013: Fix an output typo
 */
INT fasp_amg_solve_amli (AMG_data *mgl,
                         AMG_param *param)
{
    dCSRmat     *ptrA = &mgl[0].A;
    dvector     *b = &mgl[0].b, *x = &mgl[0].x, *r = &mgl[0].w;
    
    const INT    MaxIt       = param->maxit;
    const SHORT  print_level = param->print_level;
    const REAL   tol         = param->tol;
    const REAL   sumb        = fasp_blas_dvec_norm2(b); // L2norm(b)
    
    // local variables
    REAL         solve_start, solve_end, solve_time;
    REAL         relres1 = BIGREAL, absres0 = sumb, absres, factor;
    INT          iter = 0;
    
#if DEBUG_MODE
    printf("### DEBUG: %s ...... [Start]\n", __FUNCTION__);
    printf("### DEBUG: nrow = %d, ncol = %d, nnz = %d\n",
           mgl[0].A.row, mgl[0].A.col, mgl[0].A.nnz);
#endif
    
    fasp_gettime(&solve_start);

    // Print iteration information if needed
    print_itinfo(print_level, STOP_REL_RES, iter, 1.0, sumb, 0.0);
    
    // MG solver here
    while ( (++iter <= MaxIt) & (sumb > SMALLREAL) ) {
        
        // Call one AMLI cycle
        fasp_solver_amli(mgl, param, 0);
        
        // Form residual r = b-A*x
        fasp_dvec_cp(b, r);
        fasp_blas_dcsr_aAxpy(-1.0, ptrA, x->val, r->val);
        
        // Compute norms of r and convergence factor
        absres  = fasp_blas_dvec_norm2(r); // residual ||r||
        relres1 = absres/sumb;             // relative residual ||r||/||b||
        factor  = absres/absres0;          // contraction factor
        absres0 = absres;                  // prepare for next iteration
        
        // Print iteration information if needed
        print_itinfo(print_level, STOP_REL_RES, iter, relres1, absres, factor);
        
        // Check convergence
        if ( relres1 < tol ) break;
    }
    
    if ( print_level > PRINT_NONE ) {
        ITS_FINAL(iter, MaxIt, relres1);
        fasp_gettime(&solve_end);
        solve_time = solve_end - solve_start;
        print_cputime("AMLI solve", solve_time);
    }
    
#if DEBUG_MODE
    printf("### DEBUG: %s ...... [Finish]\n", __FUNCTION__);
#endif
    
    return iter;
}

/**
 * \fn INT fasp_amg_solve_nl_amli (AMG_data *mgl, AMG_param *param)
 *
 * \brief Nonlinear AMLI -- SOLVE phase
 *
 * \param mgl    Pointer to AMG data: AMG_data
 * \param param  Pointer to AMG parameters: AMG_param
 *
 * \return       Iteration number if succeed, ERROR otherwise
 *
 * \author Xiaozhe Hu
 * \date   04/30/2011
 *
 * Modified by Chensong 04/21/2013: Fix an output typo
 *
 * \note Nonlinar AMLI-cycle.  Refer to Xiazhe Hu, Panayot S. Vassilevski, Jinchao Xu
 *          "Comparative Convergence Analysis of Nonlinear AMLI-cycle Multigrid", 2013.
 *
 */
INT fasp_amg_solve_nl_amli (AMG_data *mgl,
                            AMG_param *param)
{
    dCSRmat      *ptrA = &mgl[0].A;
    dvector      *b = &mgl[0].b, *x = &mgl[0].x, *r = &mgl[0].w;
    
    const INT     MaxIt       = param->maxit;
    const SHORT   print_level = param->print_level;
    const REAL    tol         = param->tol;
    const REAL    sumb        = fasp_blas_dvec_norm2(b); // L2norm(b)
    
    // local variables
    REAL          solve_start, solve_end;
    REAL          relres1 = BIGREAL, absres0 = BIGREAL, absres, factor;
    INT           iter = 0;
    
#if DEBUG_MODE
    printf("### DEBUG: %s ...... [Start]\n", __FUNCTION__);
    printf("### DEBUG: nrow = %d, ncol = %d, nnz = %d\n",
           mgl[0].A.row, mgl[0].A.col, mgl[0].A.nnz);
#endif
    
    fasp_gettime(&solve_start);
    
    // Print iteration information if needed
    print_itinfo(print_level, STOP_REL_RES, iter, 1.0, sumb, 0.0);
    
    while ( (++iter <= MaxIt) & (sumb > SMALLREAL) ) // MG solver here
    {
        // one multigrid cycle
        fasp_solver_nl_amli(mgl, param, 0, mgl[0].num_levels);
        
        // r = b-A*x
        fasp_dvec_cp(b, r);
        fasp_blas_dcsr_aAxpy(-1.0, ptrA, x->val, r->val);
        
        absres  = fasp_blas_dvec_norm2(r); // residual ||r||
        relres1 = absres/sumb;       // relative residual ||r||/||b||
        factor  = absres/absres0;    // contraction factor
        
        // output iteration information if needed
        print_itinfo(print_level, STOP_REL_RES, iter, relres1, absres, factor);
        
        if ( relres1 < tol ) break; // early exit condition
        
        absres0 = absres;
    }
    
    if ( print_level > PRINT_NONE ) {
        ITS_FINAL(iter, MaxIt, relres1);
        fasp_gettime(&solve_end);
        print_cputime("Nonlinear AMLI solve", solve_end - solve_start);
    }
    
#if DEBUG_MODE
    printf("### DEBUG: %s ...... [Finish]\n", __FUNCTION__);
#endif
    
    return iter;
}

/**
 * \fn void fasp_famg_solve (AMG_data *mgl, AMG_param *param)
 *
 * \brief FMG -- SOLVE phase
 *
 * \param mgl    Pointer to AMG data: AMG_data
 * \param param  Pointer to AMG parameters: AMG_param
 *
 * \author Chensong Zhang
 * \date   01/10/2012
 */
void fasp_famg_solve (AMG_data *mgl,
                      AMG_param *param)
{
    dCSRmat     *ptrA = &mgl[0].A;
    dvector     *b = &mgl[0].b, *x = &mgl[0].x, *r = &mgl[0].w;
    
    const SHORT  print_level = param->print_level;
    const REAL   sumb        = fasp_blas_dvec_norm2(b); // L2norm(b)
    
    // local variables
    REAL         solve_start, solve_end;
    REAL         relres1 = BIGREAL, absres;
        
#if DEBUG_MODE
    printf("### DEBUG: %s ...... [Start]\n", __FUNCTION__);
    printf("### DEBUG: nrow = %d, ncol = %d, nnz = %d\n",
           mgl[0].A.row, mgl[0].A.col, mgl[0].A.nnz);
#endif
    
    fasp_gettime(&solve_start);

    // Call one full multigrid cycle
    fasp_solver_fmgcycle(mgl, param);
    
    // Form residual r = b-A*x
    fasp_dvec_cp(b, r);
    fasp_blas_dcsr_aAxpy(-1.0, ptrA, x->val, r->val);
    
    // Compute norms of r and convergence factor
    absres  = fasp_blas_dvec_norm2(r); // residual ||r||
    relres1 = absres/sumb;             // relative residual ||r||/||b||
    
    if ( print_level > PRINT_NONE ) {
        printf("FMG finishes with relative residual %e.\n", relres1);
        fasp_gettime(&solve_end);
        print_cputime("FMG solve",solve_end - solve_start);
    }
    
#if DEBUG_MODE
    printf("### DEBUG: %s ...... [Finish]\n", __FUNCTION__);
#endif
    
    return;
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
